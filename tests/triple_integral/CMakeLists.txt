include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(triple_instable ../Random.cpp triple.cpp)
add_executable(ftriple ../Random.cpp ftriple.cpp)

target_compile_features(triple_instable PRIVATE cxx_range_for)
target_compile_features(ftriple PRIVATE cxx_range_for)

target_link_libraries(triple_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)
