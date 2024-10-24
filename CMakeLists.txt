# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project( pro12 )

# Specify the path to the Boost folder
if(DEFINED ENV{BOOST_ROOT})
    set(Boost_NO_SYSTEM_PATHS TRUE)
    set(BOOST_ROOT $ENV{BOOST_ROOT})
endif()

# Check for CGAL
if(DEFINED ENV{CGAL_ROOT})
    set(CGAL_DIR $ENV{CGAL_ROOT})
endif()

set(CMAKE_BUILD_TYPE Release)

# Add the include directory for header files
include_directories(${CMAKE_SOURCE_DIR}/includes/utils)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  Qt5)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# Boost and its components
find_package(Boost REQUIRED COMPONENTS json)  # Specify 'json' here

# Boost and its components
#find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

include_directories(${Boost_INCLUDE_DIRS})
# include for local package


# Creating entries for target: pro12
# ############################

add_executable(pro12 pro12.cpp includes/utils/functions.cpp)

add_to_cached_list( CGAL_EXECUTABLE_TARGETS pro12 )

# Link the executable to CGAL and third-party libraries
target_link_libraries(pro12 PUBLIC Qt5::Widgets Qt5::Gui Qt5::Core CGAL::CGAL Boost::boost Boost::json)

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER)
  target_link_libraries(pro12
 PRIVATE CGAL::CGAL_Qt5)
endif()

