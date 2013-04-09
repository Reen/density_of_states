cmake_minimum_required(VERSION 2.6)

project( gengraph )

include_directories("include/")

file(GLOB src_files include/*.cpp)
add_executable(graph mr-connected.cpp ${src_files})
add_executable(distrib distrib.cpp ${src_files})

set_property(TARGET graph PROPERTY OUTPUT_NAME_DEBUG "graph_debug")
set_property(TARGET graph PROPERTY DEBUG_OUTPUT_NAME "graph_debug")
set_property(TARGET graph PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set_property(TARGET distrib PROPERTY OUTPUT_NAME_DEBUG "distrib_debug")
set_property(TARGET distrib PROPERTY DEBUG_OUTPUT_NAME "distrib_debug")
set_property(TARGET distrib PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
