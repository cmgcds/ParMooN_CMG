
include_directories("${CMAKE_SOURCE_DIR}/include/ANN")

list(APPEND ANN_SOURCES "${PROJECT_SOURCE_DIR}/src/ANN/ANN.C")

add_library(ann STATIC ${ANN_SOURCES})
