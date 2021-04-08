#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "mlpack::mlpack" for configuration ""
set_property(TARGET mlpack::mlpack APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(mlpack::mlpack PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "/home/subodh/bin/mlpack/lib/libmlpack.so.3.4"
  IMPORTED_SONAME_NOCONFIG "libmlpack.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS mlpack::mlpack )
list(APPEND _IMPORT_CHECK_FILES_FOR_mlpack::mlpack "/home/subodh/bin/mlpack/lib/libmlpack.so.3.4" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
