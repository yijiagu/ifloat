include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(fft_instable ../Random.cpp FFT.cpp)
add_executable(ffft ../Random.cpp fFFT.cpp)

target_compile_features(fft_instable PRIVATE cxx_range_for)
target_compile_features(ffft PRIVATE cxx_range_for)

target_link_libraries(fft_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)
