#!/bin/bash

python3 -m venv pyenv && . pyenv/bin/activate
pip3 install --upgrade pip
pip3 install coveralls 'cython<3.3.0' gcovr==8.3 matplotlib requests scikit-build-core tabulate numpy networkx ipycytoscape seaborn plotly anywidget

mkdir core_build && cd "$_"
export CPU_COUNT=$(python3 $GITHUB_WORKSPACE/.github/workflows/scripts/get_core_count.py)
cmake -DNETWORKIT_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug -DNETWORKIT_COVERAGE=ON ..
make -j$CPU_COUNT

cd ..
./core_build/networkit_tests -t

export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:+$CMAKE_LIBRARY_PATH:}$(pwd)/core_build
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}$(pwd)/core_build

CMAKE_BUILD_PARALLEL_LEVEL=$CPU_COUNT pip3 install --no-build-isolation -e .