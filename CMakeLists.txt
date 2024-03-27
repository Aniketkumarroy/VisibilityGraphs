cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(main VERSION 1.0)

add_subdirectory(visgraph)

set(EXECUTABLE main.cpp) # your source file name

get_filename_component(TARGET ${EXECUTABLE} NAME_WE)
add_executable(${TARGET} ${EXECUTABLE})
target_link_libraries(${TARGET} PRIVATE ${VISGRAPH_LIBS} m)
target_include_directories(${TARGET} PRIVATE ${VISGRAPH_INCLUDE_DIRS})