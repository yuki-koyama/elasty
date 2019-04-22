cmake_minimum_required(VERSION 3.1)
project(elasty)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/elasty/*.hpp)

add_library(elasty STATIC ${sources} ${headers})
target_include_directories(elasty PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(elasty PUBLIC Eigen3::Eigen tinyobjloader)

option(ELASTY_EXAMPLES ON "")
if(ELASTY_EXAMPLES)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/bigger)

  add_executable(simple ${CMAKE_CURRENT_SOURCE_DIR}/examples/simple/main.cpp)
  target_link_libraries(simple bigger elasty)
  add_custom_command(TARGET simple POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${BIGGER_SHADER_OUTPUT_DIR} $<TARGET_FILE_DIR:simple>/shaders)
  add_custom_command(TARGET simple POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/models $<TARGET_FILE_DIR:simple>/models)
endif()

option(ELASTY_TESTS ON "")
if(ELASTY_TESTS)
    enable_testing()

    add_executable(test-bending-constraint ${CMAKE_CURRENT_SOURCE_DIR}/tests/test-bending-constraint.cpp)
    target_link_libraries(test-bending-constraint elasty)

    add_test(NAME test-bending-constraint COMMAND $<TARGET_FILE:test-bending-constraint>)
endif()
