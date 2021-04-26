# ===================================================================
# This is FindMLPACK.cmake file for the ParMooN Version 1.1
# written by Thivin Anandh 
# date: 22 Apr 2020
# searching for a MLPACK lib in the system 
# if found, this will define
#  MLPACK_FOUND - System has MLPACK
#  MLPACK_INCLUDE_DIRS - The MLPACK include directories
#  MLPACK_LIBRARIES - The libraries needed to use MLPACK
# ===================================================================
if(MLPACK_INCLUDES AND MLPACK_LIBRARIES)
  set(MLPACK_FIND_QUIETLY TRUE)
endif(MLPACK_INCLUDES AND MLPACK_LIBRARIES)

if(NOT MLPACK_FOUND)

  if(NOT MLPACK_LIBRARY)
    message("MLPACK not found in the system, so checking the availability in ParMooN for the selected AParMooN_ARCH=${AParMooN_ARCH}")
    find_path(MLPACK_INCLUDE_DIR  NAMES mlpack/core.hpp mlpack/prereqs.hpp 
                                  PATHS "${PARMOON_EXTLIB_PATH}/mlpack/include" 
                                        "${PARMOON_EXTLIB_PATH}/mlpack/include/mlpack/ensmallen"
                                        "${PARMOON_EXTLIB_PATH}/mlpack/include/mlpack/armadillo"
                                        "${PARMOON_EXTLIB_PATH}/mlpack/include/mlpack/boost"
                                        "${PARMOON_EXTLIB_PATH}/mlpack/include/mlpack/hdf5")
    
    find_library(MLPACK_LIBRARY NAMES mlpack_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/mlpack/lib)
    #  find_library(MLPACK_LIBRARY_HDF NAMES hdf5_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/mlpack/lib)
    find_library(MLPACK_LIBRARY_BOOST_SEQ NAMES boost_math_c99_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/mlpack/lib)
    find_library(MLPACK_LIBRARY_BOOST_MATH NAMES boost_serialization_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/mlpack/lib)
    find_library(MLPACK_LIBRARY_BOOST_UNIT NAMES boost_unit_test_framework_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/mlpack/lib)
    
  endif(NOT MLPACK_LIBRARY)
  
  if(MLPACK_LIBRARY)      
    # set MLPACK
    if(MLPACK_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(MLPACK_LIBRARIES ${MLPACK_LIBRARY} ${MLPACK_LIBRARY_HDF} ${MLPACK_LIBRARY_BOOST_SEQ} ${MLPACK_LIBRARY_BOOST_MATH} ${MLPACK_LIBRARY_BOOST_UNIT})
      set(MLPACK_INCLUDE_DIRS ${MLPACK_INCLUDE_DIR})
      message(${MLPACK_LIBRARIES})
      # handle the QUIETLY and REQUIRED arguments and set MLPACK_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(MLPACK  DEFAULT_MSG
                                        MLPACK_LIBRARY MLPACK_INCLUDE_DIR)

      mark_as_advanced(MLPACK_INCLUDE_DIR MLPACK_LIBRARY)
    endif(MLPACK_LIBRARY)  
  endif()

endif()


