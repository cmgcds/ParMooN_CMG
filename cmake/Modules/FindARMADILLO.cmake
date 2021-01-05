# ===================================================================
# This is FindARMADILLO.cmake file for the ParMooN Version 1.1
# written by Subodh Joshi, CMG, CDS, IISc 
# date: 14 Dec 2020
# Since ARMADILLO is provided with ParMooN in EXT_LIB, we do not need to 
# search it in the system path. We can directly include it here.
# if found, this will define
#  ARMADILLO_FOUND - System has ARMADILLO
#  ARMADILLO_INCLUDE_DIRS - The ARMADILLO include directories
#  ARMADILLO_LIBRARIES - The libraries needed to use ARMADILLO
# ===================================================================
if(ARMADILLO_INCLUDES AND ARMADILLO_LIBRARIES)
  set(ARMADILLO_FIND_QUIETLY TRUE)
endif(ARMADILLO_INCLUDES AND ARMADILLO_LIBRARIES)

find_path(ARMADILLO_INCLUDE_DIR 
                       NAMES armadillo
                       PATHS "${PARMOON_EXTLIB_PATH}/Armadillo/include/")
find_library(ARMADILLO_LIBRARY 
                       NAMES armadillo
                       PATHS "${PARMOON_EXTLIB_PATH}/Armadillo/lib/")
  
if(ARMADILLO_LIBRARY)
  include(FindPackageHandleStandardArgs)

  set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARY})
  set(ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIR})

  # handle the QUIETLY and REQUIRED arguments and set ARMADILLO_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(ARMADILLO  DEFAULT_MSG
                                    ARMADILLO_LIBRARY ARMADILLO_INCLUDE_DIR)

  mark_as_advanced(ARMADILLO_INCLUDE_DIR ARMADILLO_LIBRARY)
endif(ARMADILLO_LIBRARY)  



