include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(choleski_instable ../Random.cpp ALG066.cpp)
add_executable(fcholeski ../Random.cpp fALG066.cpp)

target_compile_features(choleski_instable PRIVATE cxx_range_for)
target_compile_features(fcholeski PRIVATE cxx_range_for)

target_link_libraries(choleski_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)
