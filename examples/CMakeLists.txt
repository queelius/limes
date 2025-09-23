add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE algebraic_integrators)

# Additional examples can be added here
# add_executable(benchmark benchmark.cpp)
# target_link_libraries(benchmark PRIVATE algebraic_integrators)