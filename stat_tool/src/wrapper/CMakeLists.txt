file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

add_library(_stat_tool MODULE ${SOURCES} ${HEADERS})

target_include_directories(_stat_tool PRIVATE
  ${CMAKE_SOURCE_DIR}/src/cpp
  ${Boost_INCLUDE_DIRS}
  ${Python3_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to vpstat_tool and Boost.Python
target_link_libraries(_stat_tool vpstat_tool)
wrapper_link_boost(_stat_tool)
wrapper_link_python(_stat_tool)

add_dependencies(_stat_tool vpstat_tool)

wrapper_install(_stat_tool)

