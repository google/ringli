#!/bin/sh

# Copyright 2024 The Ringli Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

DST="${1}"

if test -z "${DST}"; then
    echo "Usage: download_war_q.sh DESTINATION_DIRECTORY"
    exit 1
fi

export PYENV_ROOT="${DST}/pyenv"
if test -d "${PYENV_ROOT}"; then
    echo "pyenv v2.3.36 already checkout out in ${PYENV_ROOT}."
else
    echo "Checking out pyenv v2.3.36 into ${PYENV_ROOT}/..."
    git -c advice.detachedHead=false clone --quiet --branch v2.3.36 --single-branch https://github.com/pyenv/pyenv.git "${PYENV_ROOT}"
fi

export PYTHON39_ROOT="${PYENV_ROOT}/versions/3.9.18"
if test -d "${PYTHON39_ROOT}"; then
    echo "Python 3.9 already installed in ${PYENV_ROOT}/versions/3.9.18."
else
    echo "Installing Python 3.9 into ${PYENV_ROOT}/versions/3.9.18/..."
    "${PYENV_ROOT}/bin/pyenv" install 3.9 2> "${PYENV_ROOT}/install_log.txt"
    if test "${?}" != "0"; then
        cat "${PYENV_ROOT}/install_log.txt"
        exit 2
    fi
fi

export WARP_Q_ROOT="${DST}/WARP-Q"
if test -d "${WARP_Q_ROOT}"; then
    echo "WARP-Q v1.0.0 already checkout in ${WARP_Q_ROOT}."
else
    echo "Checking out WARP-Q v1.0.0 into ${WARP_Q_ROOT}..."
    git -c advice.detachedHead=false clone --quiet --branch v1.0.0 --single-branch https://github.com/wjassim/WARP-Q.git "${WARP_Q_ROOT}"
fi

echo "Ensuring WARP-Q dependencies installed into ${PYENV_ROOT}/versions/3.9.18/lib/python3.9/site-packages/..."
"${PYTHON39_ROOT}/bin/pip" install -r "${WARP_Q_ROOT}/requirements.txt" > "${WARP_Q_ROOT}/install_log.txt" 2>&1
if test "${?}" != "0"; then
    cat "${WARP_Q_ROOT}/install_log.txt"
    exit 3
fi

WARP_Q_SCRIPT="${DST}/warp_q.sh"
echo "Dropping executable warp_q.sh script into ${WARP_Q_SCRIPT}..."
cat > "${WARP_Q_SCRIPT}" <<EOF
#!/bin/sh

USAGE="Usage: warp_q.sh REFERENCE_WAV_PATH DISTORTION_WAV_PATH"

if test ! -f "\${1}"; then
    echo "\${USAGE}"
    exit 1
fi

if test ! -f "\${2}"; then
    echo "\${USAGE}"
    exit 1
fi

"${PYTHON39_ROOT}/bin/python3.9" "${WARP_Q_ROOT}/warpq.py" --mode predict_file --org "\${1}" --deg "\${2}" --mapping_model "${WARP_Q_ROOT}/models/RandomForest/Genspeech.zip"
EOF

chmod +x "${WARP_Q_SCRIPT}"
