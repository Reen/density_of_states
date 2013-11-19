cmake_minimum_required(VERSION 2.6)

project( dos )

set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -DDEBUG -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -g -DNDEBUG -DBOOST_UBLAS_NDEBUG -fomit-frame-pointer -mtune=native -march=native -ftree-vectorize")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug ) 
endif( NOT CMAKE_BUILD_TYPE )

add_subdirectory(src/)
add_subdirectory(gengraph/)
