# Find Google Test
find_package(GTest REQUIRED)

# Enable testing
enable_testing()
include(GoogleTest)

# Coverage setup for GCC/Clang
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Coverage enabled for tests")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

# Keep basic tests for simple validation
add_executable(basic_tests basic_tests.cpp)
target_link_libraries(basic_tests PRIVATE algebraic_integrators)
add_test(NAME BasicTests COMMAND basic_tests)
set_tests_properties(BasicTests PROPERTIES
    TIMEOUT 30
    LABELS "unit"
)

# Google Test-based test executables
add_executable(test_accumulators test_accumulators.cpp)
target_link_libraries(test_accumulators
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_accumulators)

add_executable(test_quadrature test_quadrature.cpp)
target_link_libraries(test_quadrature
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_quadrature)

add_executable(test_integrators test_integrators.cpp)
target_link_libraries(test_integrators
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_integrators)

add_executable(test_transforms test_transforms.cpp)
target_link_libraries(test_transforms
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_transforms)

add_executable(test_parallel test_parallel.cpp)
target_link_libraries(test_parallel
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_parallel)

add_executable(test_integration_result test_integration_result.cpp)
target_link_libraries(test_integration_result
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_integration_result)

add_executable(test_concepts test_concepts.cpp)
target_link_libraries(test_concepts
    PRIVATE
    algebraic_integrators
    GTest::gtest
    GTest::gtest_main
    pthread
)
gtest_discover_tests(test_concepts)

# Add custom target for running tests with coverage
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_custom_target(coverage
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/coverage
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/coverage
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        COMMAND lcov --capture --directory ${CMAKE_BINARY_DIR} --output-file ${CMAKE_BINARY_DIR}/coverage/coverage.info --no-external
        COMMAND lcov --remove ${CMAKE_BINARY_DIR}/coverage/coverage.info '*/tests/*' '*/build/*' '/usr/*' --output-file ${CMAKE_BINARY_DIR}/coverage/coverage.info
        COMMAND genhtml ${CMAKE_BINARY_DIR}/coverage/coverage.info --output-directory ${CMAKE_BINARY_DIR}/coverage/html
        COMMAND ${CMAKE_COMMAND} -E echo "Coverage report generated in ${CMAKE_BINARY_DIR}/coverage/html/index.html"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        DEPENDS test_accumulators test_quadrature test_integrators test_transforms test_parallel test_integration_result test_concepts
    )
endif()