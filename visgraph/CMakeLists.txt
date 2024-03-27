add_library(visgraph STATIC graph.cpp ThreadPool.cpp vis_graph.cpp visgraph.cpp visible_vertices.cpp)
target_include_directories(visgraph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(VISGRAPH_LIBS visgraph CACHE INTERNAL "visibility graphs cpu libraries")
set(VISGRAPH_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "visibility graphs cpu include directories")