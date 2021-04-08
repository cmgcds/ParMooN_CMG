# ===================================================================
# This is FindENSMALLEN.cmake file for the ParMooN Version 1.1
# written by Subodh Joshi, CMG, CDS, IISc 
# date: 14 Dec 2020
# Since ENSMALLEN is provided with ParMooN in EXT_LIB, we do not need to 
# search it in the system path. We can directly include it here.
# if found, this will define
#  ENSMALLEN_FOUND - System has ENSMALLEN
#  ENSMALLEN_INCLUDE_DIRS - The ENSMALLEN include directories
#  ENSMALLEN_LIBRARIES - The libraries needed to use ENSMALLEN
# ===================================================================
if(ENSMALLEN_INCLUDES AND ENSMALLEN_LIBRARIES)
  set(ENSMALLEN_FIND_QUIETLY TRUE)
endif(ENSMALLEN_INCLUDES AND ENSMALLEN_LIBRARIES)

find_path(ENSMALLEN_INCLUDE_DIR 
                       NAMES ensmallen.hpp
                       PATHS "${PARMOON_EXTLIB_PATH}/Ensmallen/include/")
                     #find_library(ENSMALLEN_LIBRARY NAMES ensmallen PATHS "${PARMOON_EXTLIB_PATH}/Ensmallen/lib/")
set(ENSMALLEN_INCLUDE_DIRS ${ENSMALLEN_INCLUDE_DIR})


