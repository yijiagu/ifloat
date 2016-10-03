include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(crout_instable ../Random.cpp ALG067.cpp)
add_executable(fcrout ../Random.cpp fALG067.cpp)

target_compile_features(crout_instable PRIVATE cxx_range_for)
target_compile_features(fcrout PRIVATE cxx_range_for)

target_link_libraries(crout_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)
