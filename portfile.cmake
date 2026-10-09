vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO Maxime2/tabulated-function
    REF "v${VERSION}"
    SHA512 0 # Replace with the SHA512 hash of the release archive
    HEAD_REF main
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME tabulated-function)
file(INSTALL "${SOURCE_PATH}/LICENSE" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}" RENAME copyright OPTIONAL)