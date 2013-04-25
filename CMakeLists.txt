cmake_minimum_required(VERSION 2.6)

project( dos )

set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -DDEBUG -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -fomit-frame-pointer -mtune=native")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug ) 
endif( NOT CMAKE_BUILD_TYPE )

add_subdirectory(src/)
add_subdirectory(gengraph/)
