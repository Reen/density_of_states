cmake_minimum_required(VERSION 2.6)

project( dos )

set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fomit-frame-pointer -g -mtune=native")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug ) 
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG" )
endif( NOT CMAKE_BUILD_TYPE )

add_subdirectory(src/)
add_subdirectory(gengraph/)
