cmake_minimum_required(VERSION 2.6)

project( dos )

# see http://voices.canonical.com/jussi.pakkanen/2013/03/26/a-list-of-common-cmake-antipatterns/
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -O1 -DDEBUG  -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -DBOOST_UBLAS_NDEBUG -fomit-frame-pointer -mtune=native -march=native -ftree-vectorize")

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Debug ) 
endif( NOT CMAKE_BUILD_TYPE )

add_subdirectory(src/)
add_subdirectory(gengraph/)
