# --- CMake Modules

cmake_minimum_required(VERSION 3.12)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include("Anaconda")

# --- PlantGL Project

project(treematching CXX)

# --- Build setup

set(CMAKE_INCLUDE_PATH "$ENV{CONDA_PREFIX}/include" ${CMAKE_INCLUDE_PATH})
set(CMAKE_LIBRARY_PATH "$ENV{CONDA_PREFIX}/lib" ${CMAKE_LIBRARY_PATH})
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)

if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")



# ---  CXX11 Compilation

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")


# --- (Win32) Multithreaded Compilation

if (MSVC)
    string(REGEX REPLACE "/W3" "/W0" ${CMAKE_CXX_FLAGS} "${${CMAKE_CXX_FLAGS}}") 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")
endif()



## ###################################################################
## Dependencies 
## ###################################################################

# --- Python
set(Python3_FIND_VIRTUALENV FIRST)
if (WIN32)
  # needed when we run cmake in a conda environment
  set(Python3_FIND_REGISTRY LAST)
endif()


find_package (Python3 COMPONENTS Interpreter Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})
message(STATUS "Python version : " ${Python3_VERSION})
message(STATUS "Python lib : " ${Python3_LIBRARIES})

# --- Libraries

set(Boost_NO_SYSTEM_PATHS ON)
set(Boost_USE_MULTITHREAD ON)
set(Boost_USE_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)

set(boost_python python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})
find_package(Boost 1.69 COMPONENTS ${boost_python} )

# --- Include Directories

include_directories("src/cpp")
include_directories(${Boost_INCLUDE_DIR})

# --- Library Directory

if (DEFINED CONDA_ENV)
    link_directories("${CONDA_ENV}/lib")
endif()

# --- Source Directories

add_subdirectory("src/cpp")
add_subdirectory("src/wrapper")

