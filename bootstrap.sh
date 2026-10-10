#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

VCPKG_ROOT="${SCRIPT_DIR}/vcpkg"

# 1. Clone vcpkg if not present
if [ ! -d "${VCPKG_ROOT}" ]; then
    echo "==> Cloning vcpkg..."
    git clone https://github.com/microsoft/vcpkg.git "${VCPKG_ROOT}"
fi

# 2. Bootstrap vcpkg executable
if [ ! -f "${VCPKG_ROOT}/vcpkg" ]; then
    echo "==> Bootstrapping vcpkg..."
    "${VCPKG_ROOT}/bootstrap-vcpkg.sh" -disableMetrics
fi

# 3. Configure CMake with vcpkg toolchain
echo "==> Configuring CMake with vcpkg toolchain..."
cmake -B build -S . \
    -DCMAKE_TOOLCHAIN_FILE="${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON

# Symlink compile_commands.json to root for VS Code / language servers
if [ -f "build/compile_commands.json" ]; then
    ln -sf build/compile_commands.json compile_commands.json
fi

echo "==> Building project and test suite..."
cmake --build build --config Release
ctest --test-dir build --output-on-failure --build-config Release