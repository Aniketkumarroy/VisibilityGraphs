cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(main VERSION 1.0)

add_subdirectory(visgraph)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE ${VISGRAPH_LIBS} m)
target_include_directories(main PRIVATE ${VISGRAPH_INCLUDE_DIRS})