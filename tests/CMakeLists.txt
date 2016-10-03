include_directories(${PROJECT_SOURCE_DIR}/lib/include)
include_directories(${PROJECT_SOURCE_DIR}/lib/src)
include_directories(${PROJECT_SOURCE_DIR}/tests)

add_subdirectory(ray_tracing)
add_subdirectory(sor)
add_subdirectory(fft)
add_subdirectory(ldl)
add_subdirectory(nbody)
#add_subdirectory(reflection)
add_subdirectory(triple_integral)
add_subdirectory(adam_forth)
add_subdirectory(choleski)
add_subdirectory(crout)
add_subdirectory(stencil2d)

option(test "Build all tests." OFF)

if (test)
  add_subdirectory(gtest)

  enable_testing()

  # Include the gtest library. gtest_SOURCE_DIR is available due to
  # 'project(gtest)' above.
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

endif()
