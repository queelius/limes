cmake_minimum_required(VERSION 3.20)
project(algebraic_integrators VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_TESTS "Build test programs" ON)
option(ENABLE_AVX2 "Enable AVX2 SIMD instructions" ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Header-only library
add_library(algebraic_integrators INTERFACE)
add_library(ai::algebraic_integrators ALIAS algebraic_integrators)

target_include_directories(algebraic_integrators
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Compiler flags
target_compile_features(algebraic_integrators INTERFACE cxx_std_20)

# Platform-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(algebraic_integrators INTERFACE
        $<$<CONFIG:Release>:-O3 -march=native -ffast-math>
        $<$<CONFIG:Debug>:-O0 -g -fsanitize=address,undefined>
    )

    if(ENABLE_AVX2)
        target_compile_options(algebraic_integrators INTERFACE -mavx2 -mfma)
    endif()
endif()

if(MSVC)
    target_compile_options(algebraic_integrators INTERFACE
        $<$<CONFIG:Release>:/O2 /arch:AVX2>
        $<$<CONFIG:Debug>:/Od /MDd>
    )
endif()

# OpenMP support
if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(algebraic_integrators INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

# Threading support
find_package(Threads REQUIRED)
target_link_libraries(algebraic_integrators INTERFACE Threads::Threads)

# Examples
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS algebraic_integrators
    EXPORT algebraic_integratorsTargets
)

install(EXPORT algebraic_integratorsTargets
    FILE algebraic_integratorsTargets.cmake
    NAMESPACE ai::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/algebraic_integrators
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package configuration
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/algebraic_integratorsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/algebraic_integratorsConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/algebraic_integrators
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/algebraic_integratorsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/algebraic_integratorsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/algebraic_integratorsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/algebraic_integrators
)