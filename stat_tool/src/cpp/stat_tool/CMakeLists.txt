# Glob includes and sources
file(GLOB HEADERS "*.h" "*.hpp")
file(GLOB SOURCES "*.cpp")

add_library(vpstat_tool SHARED ${SOURCES})

target_compile_definitions(vpstat_tool PRIVATE MESSAGE)

# Optionally add efence
if(WITH_EFENCE)
  target_compile_definitions(vpstat_tool PRIVATE DEBUG)
  target_link_libraries(vpstat_tool PRIVATE efence)
endif()

target_include_directories(vpstat_tool PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

set(STAT_TOOL_LIBDIR "${SKBUILD_PLATLIB_DIR}/openalea/stat_tool")

message(STATUS "Installing vpstat_tool to ${STAT_TOOL_LIBDIR}")
install(TARGETS vpstat_tool
  LIBRARY DESTINATION "${STAT_TOOL_LIBDIR}"
)

message(STATUS "Installing headers: ${HEADERS}")
message(STATUS "Installing headers to ${SKBUILD_HEADERS_DIR}/stat_tool")
install(FILES ${HEADERS}
  DESTINATION "${SKBUILD_HEADERS_DIR}/stat_tool"
)
