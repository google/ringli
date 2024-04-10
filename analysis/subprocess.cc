// Copyright 2024 The Ringli Authors. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "analysis/subprocess.h"

#include <fcntl.h>
#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"

namespace ringli {

typedef struct pollfd pollfd;

SubprocessResult Execute(const std::string& binary_path,
                         const std::vector<std::string>& arguments,
                         const std::string& input) {
  int stdin[2] = {0};
  int stdout[2] = {0};
  int stderr[2] = {0};
  const int READ = 0;
  const int WRITE = 1;

  if (pipe(stdin) == -1) {
    return SubprocessResult{
        .stderr =
            absl::StrCat("failed to create stdin pipe: ", strerror(errno)),
        .status = -1};
  }
  if (pipe(stdout) == -1) {
    return SubprocessResult{
        .stderr =
            absl::StrCat("failed to create stdout pipe: ", strerror(errno)),
        .status = -1};
  }
  if (pipe(stderr) == -1) {
    return SubprocessResult{
        .stderr =
            absl::StrCat("failed to create stderr pipe: ", strerror(errno)),
        .status = -1};
  }

  const pid_t child_pid = fork();
  if (child_pid == -1) {
    return SubprocessResult{.stderr = "failed to fork", .status = -1};
  } else if (child_pid == 0) {
    close(stdin[WRITE]);
    close(stdout[READ]);
    close(stderr[READ]);

    if (dup2(stdin[READ], STDIN_FILENO) == -1) {
      std::cerr << "failed to dup stdin";
      exit(-1);
    }
    close(stdin[READ]);

    if (dup2(stdout[WRITE], STDOUT_FILENO) == -1) {
      std::cerr << "failed to dup stdout";
      exit(-1);
    }
    close(stdout[WRITE]);

    if (dup2(stderr[WRITE], STDERR_FILENO) == -1) {
      std::cerr << "failed to dup stderr";
      exit(-1);
    }
    close(stderr[WRITE]);

    std::vector<std::vector<char>> argument_vectors;
    argument_vectors.reserve(arguments.size());
    std::vector<char*> argument_pointers;
    argument_pointers.reserve(arguments.size() + 1);
    const auto push_argument = [&](const std::string& argument) {
      argument_vectors.emplace_back(argument.begin(), argument.end());
      argument_vectors.back().push_back('\0');
      argument_pointers.push_back(argument_vectors.back().data());
    };
    push_argument(binary_path);
    for (const auto& argument : arguments) {
      push_argument(argument);
    }
    argument_pointers.push_back(nullptr);
    execvp(binary_path.c_str(), argument_pointers.data());

    std::cerr << "failed to execute " << binary_path << std::endl;
    exit(-1);
  } else {
    close(stdin[READ]);
    close(stdout[WRITE]);
    close(stderr[WRITE]);
    int input_offset = 0;
    std::vector<char> buffer(1024);
    std::vector<char> out;
    std::vector<char> err;

    if (fcntl(stdin[WRITE], F_SETFL, O_NONBLOCK) < 0) {
      return SubprocessResult{
          .stderr = absl::StrCat("failed making child stdin nonblocking: ",
                                 strerror(errno)),
          .status = -1};
    }

    std::vector<pollfd> poll_fds;
    poll_fds.push_back({.fd = stdin[WRITE], .events = POLLOUT});
    poll_fds.push_back({.fd = stdout[READ], .events = POLLIN});
    poll_fds.push_back({.fd = stderr[READ], .events = POLLIN});

    if (input.empty()) {
      poll_fds[0].events = 0;
      poll_fds[0].fd = -1;
      close(stdin[WRITE]);
    }

    while (poll_fds[0].fd >= 0 || poll_fds[1].fd >= 0 || poll_fds[2].fd >= 0) {
      int ready = poll(poll_fds.data(), 3, -1);
      if (ready == -1) {
        return SubprocessResult{
            .stderr =
                absl::StrCat("failed polling child pipes: ", strerror(errno)),
            .status = -1};
      }
      if (poll_fds[0].revents & POLLOUT) {
        int n_bytes = write(stdin[WRITE], input.c_str() + input_offset,
                            input.size() - input_offset);
        if (n_bytes < 0) {
          return SubprocessResult{
              .stderr = absl::StrCat(
                  "failed writing to STDIN of child process: ", errno),
              .status = -1};
        }
        input_offset += n_bytes;
        if (input_offset >= input.length()) {
          poll_fds[0].events = 0;
          close(stdin[WRITE]);
        }
      } else if (poll_fds[0].revents & (POLLHUP | POLLERR | POLLNVAL)) {
        poll_fds[0].fd = -1;
        close(stdin[WRITE]);
      }

      if (poll_fds[1].revents & POLLIN) {
        int n_bytes = read(stdout[READ], buffer.data(), buffer.size());
        if (n_bytes < 0) {
          return SubprocessResult{
              .stderr =
                  absl::StrCat("failed reading from STDOUT of child process: ",
                               strerror(errno)),
              .status = -1};
        }
        out.insert(out.end(), buffer.begin(), buffer.begin() + n_bytes);
      } else if (poll_fds[1].revents & POLLHUP) {
        poll_fds[1].fd = -1;
      }

      if (poll_fds[2].revents & (POLLIN | POLLERR | POLLNVAL)) {
        int n_bytes = read(stderr[READ], buffer.data(), buffer.size());
        if (n_bytes < 0) {
          return SubprocessResult{
              .stderr =
                  absl::StrCat("failed reading from STDERR of child process: ",
                               strerror(errno)),
              .status = -1};
        }
        err.insert(err.end(), buffer.begin(), buffer.begin() + n_bytes);
      } else if (poll_fds[2].revents & (POLLHUP | POLLERR | POLLNVAL)) {
        poll_fds[2].fd = -1;
      }
    }

    int status;
    pid_t wait_pid = wait(&status);
    if (wait_pid != child_pid) {
      return SubprocessResult{.stderr = "wrong pid returned from wait",
                              .status = -1};
    }
    return SubprocessResult{
        .stdout = std::string(out.begin(), out.end()),
        .stderr = std::string(err.begin(), err.end()),
        .status = WEXITSTATUS(status),
    };
  }
}

}  // namespace ringli