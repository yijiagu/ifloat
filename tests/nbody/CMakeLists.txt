include_directories(${PROJECT_SOURCE_DIR}/lib/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

add_executable(nbody_instable nbody.cpp)
add_executable(fnbody fnbody.cpp)

target_compile_features(nbody_instable PRIVATE cxx_range_for)
target_compile_features(fnbody PRIVATE cxx_range_for)

target_link_libraries(nbody_instable ${PROJECT_BINARY_DIR}/libinstableFloat.a mpfr ppl gmpxx gmp)
