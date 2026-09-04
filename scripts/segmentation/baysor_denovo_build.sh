#!/usr/bin/env bash
# Creates the micromamba environment used by the baysor_denovo scripts and builds
# the pinned native Baysor C++ CLI into it. Required on the LRZ Linux Cluster, which
# has no enroot image. Thread count is set by OMP_NUM_THREADS. Run once.

set -euo pipefail

ROOT="/dss/dsshome1/0C/ra98gaq/micromamba"
ENV="baysor_denovo"
ENVDIR="${ROOT}/envs/${ENV}"
WORK="${HOME}/.cache/baysor-build"
TAG="cpp-0.8.3"

[[ -d "${ENVDIR}" ]] || micromamba create -y -r "${ROOT}" -n "${ENV}" -c conda-forge \
  python=3.12 pandas pyyaml cxx-compiler cmake ninja pkg-config \
  eigen spdlog cgal-cpp libarrow libparquet hdf5 nlohmann_json libtiff

rm -rf "${WORK}"
git clone --branch "${TAG}" --depth 1 https://github.com/kharchenkolab/Baysor.git "${WORK}/src"

# Make CMake use the environment's libtiff instead of the system one, which
# requires a libjbig that the conda toolchain does not provide
mkdir -p "${WORK}/cmake"
cat > "${WORK}/cmake/FindTIFF.cmake" <<EOF
add_library(TIFF::TIFF UNKNOWN IMPORTED GLOBAL)
set_target_properties(TIFF::TIFF PROPERTIES
    IMPORTED_LOCATION "${ENVDIR}/lib/libtiff.so"
    INTERFACE_INCLUDE_DIRECTORIES "${ENVDIR}/include")
set(TIFF_FOUND TRUE)
EOF

# RPATH lets the installed binary run without activating the environment
micromamba run -r "${ROOT}" -n "${ENV}" cmake -S "${WORK}/src" -B "${WORK}/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBAYSOR_WITH_TESTS=OFF \
  -DCMAKE_INSTALL_PREFIX="${ENVDIR}" -DCMAKE_PREFIX_PATH="${ENVDIR}" \
  -DCMAKE_MODULE_PATH="${WORK}/cmake" \
  -DCMAKE_INSTALL_RPATH='$ORIGIN/../lib' -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON
micromamba run -r "${ROOT}" -n "${ENV}" cmake --build "${WORK}/build" --target baysor
micromamba run -r "${ROOT}" -n "${ENV}" cmake --install "${WORK}/build"

rm -rf "${WORK}"
echo "installed Baysor ${TAG}"
