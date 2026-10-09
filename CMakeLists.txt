cmake_minimum_required(VERSION 3.25)
project(tabulated-function
    VERSION 0.1.0
    LANGUAGES CXX
    DESCRIPTION "Tabulated function interpolation and numerical calculus library in C++23"
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

add_library(tabulated-function INTERFACE)
add_library(tabulatedfunction::tabulated-function ALIAS tabulated-function)

target_include_directories(tabulated-function
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(tabulated-function INTERFACE cxx_std_23)

install(
    TARGETS tabulated-function
    EXPORT tabulated-function-targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    FILES tabulated-function.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT tabulated-function-targets
    FILE tabulated-function-targets.cmake
    NAMESPACE tabulatedfunction::
    DESTINATION ${CMAKE_INSTALL_DATADIR}/tabulated-function
)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/tabulated-function-config.cmake.in" [=[
@PACKAGE_INIT@
include("${CMAKE_CURRENT_LIST_DIR}/tabulated-function-targets.cmake")
check_required_components(tabulated-function)
]=])

configure_package_config_file(
    "${CMAKE_CURRENT_BINARY_DIR}/tabulated-function-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/tabulated-function-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/tabulated-function
)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/tabulated-function-config.cmake"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/tabulated-function
)

if(BUILD_TESTING)
    add_executable(tabulated-function_test
        tabulated-function_test.cpp
    )
    target_link_libraries(tabulated-function_test PRIVATE tabulated-function)
    add_test(NAME TabulatedFunctionTests COMMAND tabulated-function_test)
endif()