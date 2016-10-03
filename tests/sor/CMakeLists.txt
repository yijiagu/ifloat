include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(sor_instable ../Random.cpp SOR_Instable.cpp)

add_executable(fsor ../Random.cpp fSOR.cpp)

target_compile_features(sor_instable PRIVATE cxx_range_for)
target_compile_features(fsor PRIVATE cxx_range_for)

target_link_libraries(sor_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)
