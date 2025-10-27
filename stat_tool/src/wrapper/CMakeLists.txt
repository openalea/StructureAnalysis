file(GLOB SOURCES "*.cpp")

add_library(_stat_tool MODULE ${SOURCES})

target_include_directories(_stat_tool PRIVATE
  ${CMAKE_SOURCE_DIR}/src/cpp
  ${Boost_INCLUDE_DIRS}
  ${Python3_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# cf. https://scikit-build-core.readthedocs.io/en/latest/guide/dynamic_link.html#manual-patching
if(APPLE)
  set(origin_token "@loader_path")
else()
  set(origin_token "$ORIGIN")
endif()

# Link to vpstat_tool and Boost.Python
target_link_libraries(_stat_tool vpstat_tool)
target_compile_definitions(_stat_tool PRIVATE BOOST_ALL_NO_LIB)
target_link_libraries(_stat_tool ${Boost_LIBRARIES})
if(APPLE)
  target_link_libraries(_stat_tool "-undefined dynamic_lookup")
else()
  target_link_libraries(_stat_tool Python3::Python)
endif()

## wrapper_link_boost(_stat_tool)
## wrapper_link_python(_stat_tool)

set_target_properties(_stat_tool PROPERTIES
  INSTALL_RPATH "${origin_token}"
)

add_dependencies(_stat_tool vpstat_tool)
## wrapper_install(_stat_tool stat_tool)
message(STATUS "Installing wrappers to ${SKBUILD_PLATLIB_DIR}/openalea/stat_tool")
set_target_properties(_stat_tool PROPERTIES PREFIX "")
install(TARGETS _stat_tool
  LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/openalea/stat_tool"
)
