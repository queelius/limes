cmake_minimum_required(VERSION 3.20)
project(limes VERSION 3.1.0 LANGUAGES CXX)

option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_TESTS "Build test programs" ON)
option(BUILD_DOCS "Build API documentation with Doxygen" OFF)
option(ENABLE_AVX2 "Enable AVX2 SIMD instructions" ON)
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)

include(cmake/LimesCompilerOptions.cmake)

# Header-only library
add_library(limes INTERFACE)
add_library(limes::limes ALIAS limes)

target_include_directories(limes
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(limes INTERFACE cxx_std_20)
limes_apply_compiler_options(limes)

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(limes INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

find_package(Threads REQUIRED)
target_link_libraries(limes INTERFACE Threads::Threads)

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Documentation
if(BUILD_DOCS)
    find_package(Doxygen REQUIRED)

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        @ONLY
    )

    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )

    message(STATUS "Build docs with: cmake --build build --target docs")
endif()

# Installation
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS limes EXPORT limesTargets)

install(EXPORT limesTargets
    FILE limesTargets.cmake
    NAMESPACE limes::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/limes
)

install(DIRECTORY include/limes
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/limesConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/limesConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/limes
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/limesConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/limesConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/limesConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/limes
)
