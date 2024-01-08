# --- Source Files

file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(tree_matching SHARED ${SRC_FILES})

# --- Linked Libraries


# --- Preprocessor

if (WIN32)
    # Export DLL on Windows
    target_compile_definitions(tree_matching PRIVATE TREEMATCH_MAKEDLL)
endif()

if (APPLE)
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()

# --- Output Library

install(TARGETS tree_matching LIBRARY DESTINATION "lib")

# --- Install Headers

install(DIRECTORY "." DESTINATION "include/tree_matching" FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
