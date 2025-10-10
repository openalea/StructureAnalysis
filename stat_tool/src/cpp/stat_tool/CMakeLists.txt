# Glob includes and sources
file(GLOB HEADERS "*.h" "*.hpp")
file(GLOB SOURCES "*.cpp")

add_library(vpstat_tool STATIC ${SOURCES} ${HEADERS})

target_compile_definitions(vpstat_tool PRIVATE MESSAGE)

# Optionally add efence
if(WITH_EFENCE)
  target_compile_definitions(vpstat_tool PRIVATE DEBUG)
  target_link_libraries(vpstat_tool PRIVATE efence)
endif()

# Add include directories
target_include_directories(vpstat_tool PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
