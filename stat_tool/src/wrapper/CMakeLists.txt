file(GLOB SOURCES "*.cpp")


target_sources(
  _stat_tool
  PRIVATE
    ${SOURCES}
)


# # cf. https://scikit-build-core.readthedocs.io/en/latest/guide/dynamic_link.html#manual-patching
# if(APPLE)
#   set(origin_token "@loader_path")
# else()
#   set(origin_token "$ORIGIN")
# endif()


# set_target_properties(_stat_tool PROPERTIES
#  INSTALL_RPATH "${origin_token}"
#)

