function(wrapper_link_python libwrapname)
  if(Python3_FOUND)
    if(APPLE)
      target_link_libraries(${libwrapname} "-undefined dynamic_lookup")
    else()
      target_link_libraries(${libwrapname} Python3::Python)
    endif()
  elseif(Python2_FOUND)
    target_link_libraries(${libwrapname} Python2::Python)
  endif()
endfunction()

function(wrapper_link_boost libwrapname)
  target_compile_definitions(${libwrapname} PRIVATE BOOST_ALL_NO_LIB)
  target_link_libraries(${libwrapname} ${Boost_LIBRARIES})
endfunction()

function(wrapper_install libwrapname libname)
  set_target_properties(${libwrapname} PROPERTIES PREFIX "")

  ##   if(WIN32)
  ##     set_target_properties(${libwrapname} PROPERTIES SUFFIX ".pyd")
  ##   elseif(APPLE)
  ##     set_target_properties(${libwrapname} PROPERTIES SUFFIX ".so")
  ##   endif()

  install(TARGETS ${libwrapname}
    LIBRARY DESTINATION "${SKBUILD_PLATLIB_DIR}/openalea/${libname}"
  )

endfunction()

function(install_headers directory headers)
  message(STATUS "Installing headers to ${SKBUILD_HEADERS_DIR}/stat_tool")
  install(FILES ${headers}
    DESTINATION "${SKBUILD_HEADERS_DIR}/stat_tool"
)
endfunction()
