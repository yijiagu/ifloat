include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(adam_instable ALG054.cpp)

add_executable(fadam fALG054.cpp)

target_compile_features(adam_instable PRIVATE cxx_range_for)

target_compile_features(fadam PRIVATE cxx_range_for)

target_link_libraries(adam_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})




