
include_directories("${CMAKE_SOURCE_DIR}/include/Solvers")
include_directories("${CMAKE_SOURCE_DIR}/include/Solvers/cuda_general")
include_directories("/usr/local/cuda/include")



list(APPEND CUSPLIB_SOURCES "${PROJECT_SOURCE_DIR}/src/Solvers/cudaSparseLinearSolvers.cu")

add_library(cudaSpLinSolver STATIC ${CUSPLIB_SOURCES})
