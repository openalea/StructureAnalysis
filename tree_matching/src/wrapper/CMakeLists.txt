# --- Source Files

file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(__tree_matching__ SHARED ${SRC_FILES})

# --- Linked Library

target_link_libraries(__tree_matching__ tree_matching)
target_link_libraries(__tree_matching__ Boost::${boost_python} Boost::dynamic_linking Boost::disable_autolinking)
if (APPLE)
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
else()
    target_link_libraries(__tree_matching__ Python3::Python)
endif()


# --- Dependencies

add_dependencies(__tree_matching__ tree_matching)

# --- Output Library

set_target_properties(__tree_matching__ PROPERTIES PREFIX "")

if (WIN32)
    set_target_properties(__tree_matching__ PROPERTIES SUFFIX ".pyd")
endif()

if (APPLE)
    set_target_properties(__tree_matching__ PROPERTIES SUFFIX ".so")
endif()

install(TARGETS __tree_matching__ DESTINATION "${CMAKE_SOURCE_DIR}/src/openalea/tree_matching")
