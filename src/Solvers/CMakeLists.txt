
include_directories("${CMAKE_SOURCE_DIR}/include/Solvers")
include_directories("${CMAKE_SOURCE_DIR}/include/Solvers/cuda_general")
include_directories("/usr/local/cuda/include")

if("${AParMooN_PARALLEL_TYPE}" STREQUAL "SCUDA" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "MCUDA" )
list(APPEND CUSPLIB_SOURCES "${PROJECT_SOURCE_DIR}/src/Solvers/cudaSparseLinearSolvers.cu")
add_library(cudaSpLinSolver STATIC ${CUSPLIB_SOURCES})
endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SCUDA"  OR  "${AParMooN_PARALLEL_TYPE}" STREQUAL "MCUDA")


if("${AParMooN_GEO}" STREQUAL "2D")
list(APPEND SYST_SOURCES "${PROJECT_SOURCE_DIR}/src/Solvers/DirectSparseLinearSolver.C")
elseif("${AParMooN_GEO}" STREQUAL "3D")

endif("${AParMooN_GEO}" STREQUAL "2D")
add_library(solver STATIC ${SYST_SOURCES})