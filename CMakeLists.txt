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
target_link_libraries(elasty PUBLIC Eigen3::Eigen)
