add_library(visgraph STATIC graph.cpp ThreadPool.cpp vis_graph.cpp visible_vertices.cpp)
target_include_directories(visgraph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Threads REQUIRED)
target_link_libraries(visgraph PRIVATE Threads::Threads)

# Combine visgraph and THREAD_LIBS into VISGRAPH_LIBS
set(VISGRAPH_LIBS ${VISGRAPH_LIBS} visgraph ${CMAKE_THREAD_LIBS_INIT} CACHE INTERNAL "visibility graphs cpu libraries")
set(VISGRAPH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "visibility graphs cpu include directories")
