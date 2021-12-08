#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "iLQR::iLQR" for configuration ""
set_property(TARGET iLQR::iLQR APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(iLQR::iLQR PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/iLQR"
  )

list(APPEND _IMPORT_CHECK_TARGETS iLQR::iLQR )
list(APPEND _IMPORT_CHECK_FILES_FOR_iLQR::iLQR "${_IMPORT_PREFIX}/bin/iLQR" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
