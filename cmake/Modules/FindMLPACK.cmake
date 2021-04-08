# ===================================================================
# This is FindMLPACK.cmake file for the ParMooN Version 1.1
# written by Subodh Joshi, CMG, CDS, IISc 
# date: 14 Dec 2020
# Since MLPACK is provided with ParMooN in EXT_LIB, we do not need to 
# search it in the system path. We can directly include it here.
# if found, this will define
#  MLPACK_FOUND - System has MLPACK
#  MLPACK_INCLUDE_DIRS - The MLPACK include directories
#  MLPACK_LIBRARIES - The libraries needed to use MLPACK
# ===================================================================
if(MLPACK_INCLUDES AND MLPACK_LIBRARIES)
  set(MLPACK_FIND_QUIETLY TRUE)
endif(MLPACK_INCLUDES AND MLPACK_LIBRARIES)

find_path(MLPACK_INCLUDE_DIR 
                       NAMES mlpack/core.hpp mlpack/prereqs.hpp
                       PATHS "${PARMOON_EXTLIB_PATH}/mlpack/include/")
find_library(MLPACK_LIBRARY 
                       NAMES mlpack
                       PATHS "${PARMOON_EXTLIB_PATH}/mlpack/lib/")
  
if(MLPACK_LIBRARY)
  include(FindPackageHandleStandardArgs)

  set(MLPACK_LIBRARIES ${MLPACK_LIBRARY})
  set(MLPACK_INCLUDE_DIRS ${MLPACK_INCLUDE_DIR})

  # handle the QUIETLY and REQUIRED arguments and set MLPACK_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(MLPACK  DEFAULT_MSG
                                    MLPACK_LIBRARY MLPACK_INCLUDE_DIR)

  mark_as_advanced(MLPACK_INCLUDE_DIR MLPACK_LIBRARY)
endif(MLPACK_LIBRARY)  



