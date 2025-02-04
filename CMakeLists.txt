cmake_minimum_required(VERSION 3.10)

project(Spline)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3)

add_executable(${PROJECT_NAME} src/main.cpp src/spline.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE include)

#target_link_libraries(${PROJECT_NAME} spline)