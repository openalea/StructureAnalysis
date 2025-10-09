cmake_minimum_required(VERSION 3.20)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(wrapper)

project(stat_tool LANGUAGES CXX)

# Options
option(WITH_EFENCE "Build with efence library" OFF)
option(WITH_TEST "Build tests" OFF)

# Set C++ standard and compile options
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # For shared libs

if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -enable-stdcall-fixup -enable-auto-import -enable-runtime-pseudo-reloc -s")
endif()

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)


find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
set(boost_python python${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})
find_package(Boost REQUIRED COMPONENTS ${boost_python})

if(NOT Python3_EXTENSION_SUFFIX)
  set(Python3_EXTENSION_SUFFIX ".cpython-${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}-darwin.so")
endif()

message(STATUS "Python extension suffix: ${Python3_EXTENSION_SUFFIX}")

add_subdirectory("src/cpp/stat_tool")
add_subdirectory("src/wrapper/")

# Optionally handle tests
if(WITH_TEST)
  enable_testing()
  add_subdirectory(test/cpp)
endif()
