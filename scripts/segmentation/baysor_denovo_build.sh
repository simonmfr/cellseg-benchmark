#!/usr/bin/env bash
# Builds the pinned native Baysor C++ CLI in a micromamba environment and packs it
# into a copy of the enroot image, under /opt/baysor-cpp, so the segmentation
# environment and its older Baysor stay untouched. Submit with sbatch (needs enroot,
# >=64G, ~2h); raise JOBS only with more memory.
set -euo pipefail
MAMBA="${MAMBA_EXE:-micromamba}"
ROOT="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ENV="baysor_denovo"
ENVDIR="${ROOT}/envs/${ENV}"
WORK="${HOME}/.cache/baysor-build"
TAG="cpp-0.8.3"
RUN=("${MAMBA}" run -r "${ROOT}" -n "${ENV}")
IMG="$(realpath -m "$(dirname "${BASH_SOURCE[0]}")/../../data")/misc/enroot_images"
OUT="${IMG}/benchmark_baysor.sqsh"
# Scratch lives next to the images, so it needs no local disk; squashfs cannot
# store the GPFS ACL xattrs found there, so skip them instead of warning per file.
export ENROOT_DATA_PATH="${IMG}/enroot_data"
export ENROOT_SQUASH_OPTIONS="-no-xattrs -processors ${SLURM_CPUS_PER_TASK:-8}"

cleanup() {
  rm -rf "${WORK}"
  enroot remove -f baysor_image 2>/dev/null || true
  rm -rf "${ENROOT_DATA_PATH}" "${OUT}.tmp"
}
trap cleanup EXIT

[[ -d "${ENVDIR}" ]] || "${MAMBA}" create -y -r "${ROOT}" -n "${ENV}" -c conda-forge \
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
"${RUN[@]}" cmake -S "${WORK}/src" -B "${WORK}/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBAYSOR_WITH_TESTS=OFF \
  -DCMAKE_INSTALL_PREFIX="${ENVDIR}" -DCMAKE_PREFIX_PATH="${ENVDIR}" \
  -DCMAKE_MODULE_PATH="${WORK}/cmake" \
  -DCMAKE_INSTALL_RPATH='$ORIGIN/../lib' -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON
"${RUN[@]}" cmake --build "${WORK}/build" --target baysor --parallel "${JOBS:-4}"
"${RUN[@]}" cmake --install "${WORK}/build"

mkdir -p "${ENROOT_DATA_PATH}"
enroot create -n baysor_image "${IMG}/benchmark_new.sqsh"
PREFIX="${ENROOT_DATA_PATH}/baysor_image/opt/baysor-cpp"
mkdir -p "${PREFIX}/bin" "${PREFIX}/lib"
cp -L "${ENVDIR}/bin/baysor" "${PREFIX}/bin/"
ldd "${ENVDIR}/bin/baysor" | awk -v d="${ENVDIR}/" '$3 ~ "^" d {print $3}' | sort -u \
  | xargs -I{} cp -L {} "${PREFIX}/lib/"
enroot export -o "${OUT}.tmp" baysor_image
mv -f "${OUT}.tmp" "${OUT}"

# Only on success: the environment is a build artifact, kept on failure for retries
"${MAMBA}" env remove -y -r "${ROOT}" -n "${ENV}"
"${MAMBA}" clean -y -a
echo "installed Baysor ${TAG} in ${OUT}"
