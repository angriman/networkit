#!/bin/bash
set -e
set -o pipefail

python3 -m venv pyenv && . pyenv/bin/activate
pip3 install --upgrade pip

# cython is required because git does not contain _NetworKit.
pip3 install 'cython<3.3.0' setuptools

# Several modules are required to build the documentation. The ones that affect the warnings are
# pinned, since CI treats warnings as errors.
pip3 install sphinx==9.2.0
pip3 install Jinja2 sphinx_bootstrap_theme numpydoc==1.11.0 exhale==0.3.7 nbsphinx==0.9.8 breathe==5.1.0
pip3 install sphinx_copybutton sphinx_gallery sphinx_last_updated_by_git
pip3 install ipykernel ipython matplotlib nbconvert jupyter-client networkx tabulate
pip3 install ipycytoscape plotly seaborn powerlaw anywidget

# Build the C++ core library (no need for optimizations).
mkdir core_build && cd "$_"
cmake -DCMAKE_BUILD_TYPE=Debug -DNETWORKIT_CXX_STANDARD=$CXX_STANDARD -DNETWORKIT_WARNINGS=ON "-DNETWORKIT_SPHINX_OPTIONS=${SPHINX_OPTIONS:-}" ..
make -j2
cd ..

# Build the Python extension.
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:+$CMAKE_LIBRARY_PATH:}$(pwd)/core_build
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}$(pwd)/core_build
export CPU_COUNT=$(python3 $GITHUB_WORKSPACE/.github/workflows/scripts/get_core_count.py)
NETWORKIT_PARALLEL_JOBS=$CPU_COUNT python3 ./setup.py build_ext --inplace --networkit-external-core
NETWORKIT_PARALLEL_JOBS=$CPU_COUNT pip3 install -e .

# Ensure notebook sources are available under docs/, as expected by docs/notebooks.rst.
rm -rf docs/notebooks
cp -R notebooks docs/notebooks

# Build the documentation.
cd core_build
make docs
touch htmldocs/.nojekyll
rm -rf htmldocs/{.buildinfo,.doctrees}
