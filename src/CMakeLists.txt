cmake_minimum_required(VERSION 2.6)

project( dostoy )
set (dostoy_VERSION_MAJOR 1 )
set (dostoy_VERSION_MINOR 0 )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-additional-modules/")

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in"
	"${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY
)
set(SOURCES_GENERAL "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" "${CMAKE_SOURCE_DIR}/include/GitSHA1.h")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	SET(HAVE_LINUX 1)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	SET(HAVE_DARWIN 1)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
)

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/include")

#Boost Library
set( Boost_USE_MULTITHREAD ON )

find_package( Boost COMPONENTS program_options REQUIRED )

if( NOT Boost_FOUND )
	message( fatal_error "Cannot build without boost. Please install the boost library first." )
endif( NOT Boost_FOUND )

find_package ( GSL REQUIRED )

include_directories( ${Boost_INCLUDE_DIRS} )

# TARGET: ising ###############
add_executable( ising ising.cpp simulation.cpp simulation_rng.cpp q_matrix_tools.cpp ising_system.cpp simulation_system.cpp toydos_system.cpp ${SOURCES_GENERAL} )
target_link_libraries( ising ${Boost_LIBRARIES} ${GSL_LIBRARIES})

set_property(TARGET ising PROPERTY OUTPUT_NAME_DEBUG "ising_debug")
set_property(TARGET ising PROPERTY DEBUG_OUTPUT_NAME "ising_debug")
set_property(TARGET ising PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

# TARGET: doxygen #############
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)
