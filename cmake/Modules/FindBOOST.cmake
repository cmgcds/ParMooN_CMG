# ===================================================================
# This is FindMLPACK.cmake file for the ParMooN Version 1.1
# written by Subodh Joshi, CMG, CDS, IISc 
# date: 14 Dec 2020
# Since BOOST is provided with ParMooN in EXT_LIB, we do not need to 
# search it in the system path. We can directly include it here.
# if found, this will define
#  BOOST_FOUND - System has BOOST
#  BOOST_INCLUDE_DIRS - The BOOST include directories
#  BOOST_LIBRARIES - The libraries needed to use BOOST
# ===================================================================
if(BOOST_INCLUDES AND BOOST_LIBRARIES)
  set(BOOST_FIND_QUIETLY TRUE)
endif(BOOST_INCLUDES AND BOOST_LIBRARIES)

find_path(BOOST_INCLUDE_DIR 
                       NAMES mlpack/core.hpp mlpack/prereqs.hpp
                       PATHS "${PARMOON_EXTLIB_PATH}/mlpack/include/")
find_library(BOOST_LIBRARY 
                       NAMES boost
                       PATHS "${PARMOON_EXTLIB_PATH}/mlpack/lib/")
  
if(BOOST_LIBRARY)
  include(FindPackageHandleStandardArgs)

  set(BOOST_LIBRARIES ${BOOST_LIBRARY})
  set(BOOST_INCLUDE_DIRS ${BOOST_INCLUDE_DIR})

  # handle the QUIETLY and REQUIRED arguments and set BOOST_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(BOOST  DEFAULT_MSG
                                    BOOST_LIBRARY BOOST_INCLUDE_DIR)

  mark_as_advanced(BOOST_INCLUDE_DIR BOOST_LIBRARY)
endif(BOOST_LIBRARY)  



