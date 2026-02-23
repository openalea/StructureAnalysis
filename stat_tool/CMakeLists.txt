# Initialize the CMake project
cmake_minimum_required(VERSION 3.20)

##################################################
# If we are running in a Conda environment, we automatically
# add the Conda env prefix to the CMAKE_PREFIX_PATH

if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
  #TODO: Windows Conda environments are structured differently,
  #      how unfortunate is this?
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}/Library")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include("Anaconda")

##################################################
# Define the version

if(SKBUILD)
  set(STAT_TOOL_VERSION ${SKBUILD_PROJECT_VERSION})
else()
  set(STAT_TOOL_VERSION "2.0.0")
endif()

project(openalea.stat_tool
  VERSION ${STAT_TOOL_VERSION}
  LANGUAGES CXX
  DESCRIPTION "OpenAlea Stat Tool library for structure analysis")

##################################################
# Set CMake policies for this project

# We allow <Package>_ROOT (env) variables for locating dependencies using find_paclage
cmake_policy(SET CMP0074 NEW)
# We allow target_sources to convert relative paths to absolute paths
cmake_policy(SET CMP0076 NEW)
# for Python*_FIND_STRATEGY=LOCATION
cmake_policy(SET CMP0094 NEW)
cmake_policy(SET CMP0167 NEW)
##################################################
# Initialize some default paths
# See https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include(GNUInstallDirs)

##################################################
# Set C++ standard and compile options

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS ON)

# This may be set in pyproject.toml
set(CMAKE_VERBOSE_MAKEFILE ON)

##################################################
# For Python bindings, we need to enable position-independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)  # For shared libs

##################################################
# TODO : REMOVE?
# if(WIN32)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -enable-stdcall-fixup -enable-auto-import -enable-runtime-pseudo-reloc -s")
# endif()


# RPath settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(
  FIND
  CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  isSystemDir
)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif("${isSystemDir}" STREQUAL "-1")

if (WIN32)
  string(REGEX REPLACE "/W3" "/W0" ${CMAKE_CXX_FLAGS} "${${CMAKE_CXX_FLAGS}}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")

  # To fix compilation error with vc14 and boost
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DHAVE_SNPRINTF")
endif()

# Default options (set in pyproject.toml or directly using "cmake -D")
# "option" does not override the cache value, contrary to cmake -D or pyproject setting
# These 2 lines are safety default values
option(WITH_EFENCE "Build with efence library" OFF)
option(WITH_TEST "Build tests" OFF)
message(STATUS "WITH_EFENCE = ${WITH_EFENCE}")
message(STATUS "WITH_TEST = ${WITH_TEST}")

##################################################
# Add a library for stat_tool and its Python wrappers
add_library(oastat_tool SHARED)

##################################################
# Find external dependencies
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(Boost REQUIRED COMPONENTS python)

##################################################
# Setup lib, inc, defines for C++ oastat_tool lib

target_include_directories(oastat_tool
  PUBLIC
    ${Boost_INCLUDE_DIRS}
)

# Export symbols on Windows
if (WIN32 OR MSVC)
  target_compile_definitions(oastat_tool PUBLIC STAT_TOOL_MAKEDLL)
endif()

# Optionally add efence
if(WITH_EFENCE)
  find_library(EFENCE_LIBRARY NAMES efence)
  if(EFENCE_LIBRARY)
    target_link_libraries(oastat_tool PUBLIC ${EFENCE_LIBRARY})
    target_compile_definitions(oastat_tool PUBLIC DEBUG)
  else()
    message(WARNING "Could not find efence library, disabling it")
  endif()
endif()

# Add files to build
add_subdirectory(src/cpp/stat_tool)


##################################################
# Add a Python binding library
python_add_library(_stat_tool MODULE)

target_include_directories(oastat_tool
  PUBLIC
    "src/cpp"
)

target_link_libraries(_stat_tool
  PRIVATE
    oastat_tool
    Boost::python
    Python::Module
)

target_compile_definitions(_stat_tool PRIVATE BOOST_ALL_NO_LIB)

# Control the output name of the produced shared library
set_target_properties(_stat_tool PROPERTIES PREFIX "")
set_target_properties(_stat_tool PROPERTIES OUTPUT_NAME "_stat_tool")
if(WIN32)
  set_target_properties(_stat_tool PROPERTIES SUFFIX ".pyd")
elseif(APPLE)
  set_target_properties(_stat_tool PROPERTIES SUFFIX ".so")
endif()

add_subdirectory("src/wrapper/")


if(APPLE)
  set_target_properties(_stat_tool PROPERTIES INSTALL_RPATH "@loader_path/.")
elseif(UNIX)
  set_target_properties(_stat_tool PROPERTIES INSTALL_RPATH "$ORIGIN/.")
endif()

##################################################
# Install targets and headers for stat_tool lib
install(TARGETS oastat_tool
            RUNTIME DESTINATION "${CONDA_ENV}bin/"
            LIBRARY DESTINATION "${CONDA_ENV}lib/"
            ARCHIVE DESTINATION "${CONDA_ENV}lib/"
)

install(
  TARGETS
    _stat_tool
  DESTINATION "${SKBUILD_PLATLIB_DIR}/openalea/stat_tool"
)

# Optionally handle tests
# TODO TO TEST
if(WITH_TEST)
  enable_testing()
  add_subdirectory(test/cpp)
endif()

