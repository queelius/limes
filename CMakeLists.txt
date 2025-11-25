cmake_minimum_required(VERSION 3.20)
project(calckit VERSION 2.0.0 LANGUAGES CXX)

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
add_library(calckit INTERFACE)
add_library(ck::calckit ALIAS calckit)

target_include_directories(calckit
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Compiler flags
target_compile_features(calckit INTERFACE cxx_std_20)

# Platform-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(calckit INTERFACE
        $<$<CONFIG:Release>:-O3 -march=native -ffast-math>
        $<$<CONFIG:Debug>:-O0 -g -fsanitize=address,undefined>
    )

    if(ENABLE_AVX2)
        target_compile_options(calckit INTERFACE -mavx2 -mfma)
    endif()
endif()

if(MSVC)
    target_compile_options(calckit INTERFACE
        $<$<CONFIG:Release>:/O2 /arch:AVX2>
        $<$<CONFIG:Debug>:/Od /MDd>
    )
endif()

# OpenMP support
if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(calckit INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

# Threading support
find_package(Threads REQUIRED)
target_link_libraries(calckit INTERFACE Threads::Threads)

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

install(TARGETS calckit
    EXPORT calckitTargets
)

install(EXPORT calckitTargets
    FILE calckitTargets.cmake
    NAMESPACE ck::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/calckit
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package configuration
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/calckitConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/calckitConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/calckit
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/calckitConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/calckitConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/calckitConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/calckit
)