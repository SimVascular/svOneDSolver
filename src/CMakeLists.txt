
# COMPILE CODE IN SPARSE SUBDIRECTORY
IF(sparseSolverType STREQUAL "skyline")

  # JUST ADD THE SOURCE FILES
  ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_C} ${SRC_H})

ELSEIF(sparseSolverType STREQUAL "superlu")

  # FIND BLAS
  FIND_PACKAGE(BLAS REQUIRED)
  IF(BLAS_FOUND)
    INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIRS})
  ENDIF()

  # CHECK IF SUPERLU_DIR IS STILL UNDEFINED
  IF(${SUPERLU_DIR} STREQUAL " ")
    MESSAGE(FATAL_ERROR "The SUPERLU_DIR folder is undefined.\
                         Please enter the path to the superLU_MT library.")
  ENDIF()

  # CHECK IF SUPERLU_DIR EXISTS
  IF(EXISTS ${SUPERLU_DIR})

    SET(SUPERLU_INCLUDE_DIRS "${SUPERLU_DIR}/SRC")
    SET(SUPERLU_LIBRARY "${SUPERLU_DIR}/lib/libsuperlu_mt_PTHREAD.a")

    MESSAGE(STATUS "SuperLU Includes: ${SUPERLU_INCLUDE_DIRS}")
    MESSAGE(STATUS "SuperLU Library: ${SUPERLU_LIBRARY}")

  ELSE()
    MESSAGE(FATAL_ERROR "The SUPERLU_DIR folder is invalid.\
                         Please enter the path to the superLU_MT library.")
  ENDIF()

  # ADD FILES IN ADDITIONAL FOLDERS
  FILE(GLOB SUPERLU_SRC_C "${CMAKE_SOURCE_DIR}/src/sparse/*.c"
                          "${CMAKE_SOURCE_DIR}/src/sparse/*.cpp"
                          "${CMAKE_SOURCE_DIR}/src/sparse/superlu/*.c"
                          "${CMAKE_SOURCE_DIR}/src/sparse/superlu/*.cpp")
  FILE(GLOB SUPERLU_SRC_H "${CMAKE_SOURCE_DIR}/src/sparse/*.h"
                          "${CMAKE_SOURCE_DIR}/src/sparse/superlu/*.h"
                          "${SUPERLU_INCLUDE_DIRS}/slu_mt_ddefs.h")

  # INCLUDE DIRECTORIES
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src
                      ${CMAKE_SOURCE_DIR}/src/sparse 
                      ${CMAKE_SOURCE_DIR}/src/sparse/superlu
                      ${SUPERLU_INCLUDE_DIRS})

  # LINK DIRECTORIES
  LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/src
                   ${CMAKE_SOURCE_DIR}/src/sparse 
                   ${CMAKE_SOURCE_DIR}/src/sparse/superlu)

  ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_C} ${SRC_H} ${SUPERLU_SRC_C} ${SUPERLU_SRC_H})
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${BLAS_LIBRARIES} ${SUPERLU_LIBRARY} pthread)  

ELSEIF(sparseSolverType STREQUAL "csparse")

  # ADD FILES IN ADDITIONAL FOLDERS
  FILE(GLOB CSPARSE_SRC_C "${CMAKE_SOURCE_DIR}/src/sparse/*.c"
                          "${CMAKE_SOURCE_DIR}/src/sparse/*.cpp"
                          "${CMAKE_SOURCE_DIR}/src/sparse/csparse/*.c"
                          "${CMAKE_SOURCE_DIR}/src/sparse/csparse/*.cpp")
  FILE(GLOB CSPARSE_SRC_H "${CMAKE_SOURCE_DIR}/src/sparse/*.h"
                          "${CMAKE_SOURCE_DIR}/src/sparse/csparse/*.h")

  # INCLUDE DIRECTORIES
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src
                      ${CMAKE_SOURCE_DIR}/src/sparse 
                      ${CMAKE_SOURCE_DIR}/src/sparse/csparse)

  # LINK DIRECTORIES
  LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/src
                   ${CMAKE_SOURCE_DIR}/src/sparse 
                   ${CMAKE_SOURCE_DIR}/src/sparse/csparse)

  # ADD SOURCE ON SPARSE FOLDER
  ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_C} ${SRC_H} ${CSPARSE_SRC_C} ${CSPARSE_SRC_H})

ENDIF()

install( TARGETS ${PROJECT_NAME}
         RUNTIME 
         DESTINATION bin 
         COMPONENT executable )


