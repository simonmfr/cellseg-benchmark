#!/usr/bin/env bash
# Builds the pinned native Baysor C++ CLI in a micromamba environment and appends it
# to the enroot image under /baysor-cpp. Submit with sbatch (>=64G, ~1h);
# raise JOBS only with more memory.
# The image is appended to, never unpacked: on GPFS, unsquashfs creates each file
# with its mode and the inherited ACL masks it, which silently strips the exec bit
# from every file in the image.
set -euo pipefail
MAMBA="${MAMBA_EXE:-micromamba}"
ROOT="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ENV="baysor_denovo"
ENVDIR="${ROOT}/envs/${ENV}"
WORK="${HOME}/.cache/baysor-build"
STAGE="${WORK}/stage"
TAG="cpp-0.8.3"
RUN=("${MAMBA}" run -r "${ROOT}" -n "${ENV}")
IMG="$(realpath -m "$(dirname "${BASH_SOURCE[0]}")/../../data")/misc/enroot_images"
OUT="${IMG}/benchmark_new.sqsh"
LL="${WORK}.ll"

cleanup() { rm -rf "${WORK}" "${LL}" "${OUT}.tmp"; }
trap cleanup EXIT

unsquashfs -ll "${OUT}" > "${LL}"
if grep -q 'squashfs-root/baysor-cpp' "${LL}"; then
  echo "${OUT} already carries /baysor-cpp; restore ${OUT}.bak first" >&2
  exit 1
fi

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

mkdir -p "${STAGE}/baysor-cpp/bin" "${STAGE}/baysor-cpp/lib"
cp -L "${ENVDIR}/bin/baysor" "${STAGE}/baysor-cpp/bin/"
ldd "${ENVDIR}/bin/baysor" | awk -v d="${ENVDIR}/" '$3 ~ "^" d {print $3}' | sort -u \
  | xargs -I{} cp -L {} "${STAGE}/baysor-cpp/lib/"
chmod -R a+rX,u+w "${STAGE}"
chmod 755 "${STAGE}/baysor-cpp/bin/baysor"

cp "${OUT}" "${OUT}.tmp"
mksquashfs "${STAGE}" "${OUT}.tmp" -all-root
unsquashfs -ll "${OUT}.tmp" > "${LL}"
if ! grep -qE '^-rwxr-xr-x .*squashfs-root/usr/bin/bash$' "${LL}"; then
  echo "/usr/bin/bash lost its exec bit, not promoting ${OUT}.tmp" >&2
  exit 1
fi
if ! grep -qE '^-rwxr-xr-x .*squashfs-root/baysor-cpp/bin/baysor$' "${LL}"; then
  echo "baysor missing or not executable, not promoting ${OUT}.tmp" >&2
  exit 1
fi
mv -f "${OUT}" "${OUT}.bak"
mv -f "${OUT}.tmp" "${OUT}"

"${MAMBA}" env remove -y -r "${ROOT}" -n "${ENV}"
"${MAMBA}" clean -y -a
echo "appended Baysor ${TAG} to ${OUT} (previous image at ${OUT}.bak)"
