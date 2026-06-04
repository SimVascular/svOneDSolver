# Set the library name
set(lib svoned_interface)

set(CXXSRCS interface.cpp)
set(HDRS interface.h)


if(${sparseSolverType} STREQUAL "skyline")
  file(GLOB MAIN_SRC_C "${CMAKE_SOURCE_DIR}/Code/Source/*.c"
                       "${CMAKE_SOURCE_DIR}/Code/Source/*.cxx")
  file(GLOB MAIN_SRC_H "${CMAKE_SOURCE_DIR}/Code/Source/*.h")

elseif(${sparseSolverType} STREQUAL "superlu")
  file(GLOB MAIN_SRC_C "${CMAKE_SOURCE_DIR}/Code/Source/*.c"
                       "${CMAKE_SOURCE_DIR}/Code/Source/*.cxx"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/*.cxx"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/superlu/*.c"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/superlu/*.cxx")
  file(GLOB MAIN_SRC_H "${CMAKE_SOURCE_DIR}/Code/Source/*.h"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/*.h"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/superlu/*.h")

elseif(${sparseSolverType} STREQUAL "csparse")
  file(GLOB MAIN_SRC_C "${CMAKE_SOURCE_DIR}/Code/Source/*.c"
                       "${CMAKE_SOURCE_DIR}/Code/Source/*.cxx"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/*.cxx"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/csparse/*.c"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/csparse/*.cxx")
  file(GLOB MAIN_SRC_H "${CMAKE_SOURCE_DIR}/Code/Source/*.h"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/*.h"
                       "${CMAKE_SOURCE_DIR}/Code/Source/sparse/csparse/*.h")
endif()

list(FILTER MAIN_SRC_C EXCLUDE REGEX "main.cxx")
list(FILTER MAIN_SRC_C EXCLUDE REGEX "interface.cpp")

list(APPEND MAIN_SRC_C "${CMAKE_SOURCE_DIR}/Code/Source/cvOneDSubdomain.cxx")

# Create shared library
add_library(${lib} SHARED ${CXXSRCS} ${MAIN_SRC_C} ${MAIN_SRC_H})

# Include directories
target_include_directories(${lib} PUBLIC 
  ${CMAKE_SOURCE_DIR}/Code/Source
)

# Link against required libraries
target_link_libraries(${lib} PUBLIC 
  nlohmann_json::nlohmann_json
)

# Link additional libraries for specific sparse solvers
if(${sparseSolverType} STREQUAL "superlu")
  find_package(BLAS REQUIRED)
  set(SUPERLU_INCLUDE_DIRS "${SUPERLU_DIR}/SRC")
  set(SUPERLU_LIBRARY "${SUPERLU_DIR}/lib/libsuperlu_mt_PTHREAD.a")
  
  target_include_directories(${lib} PRIVATE ${SUPERLU_INCLUDE_DIRS})
  target_link_libraries(${lib} PRIVATE ${BLAS_LIBRARIES} ${SUPERLU_LIBRARY} pthread)
endif()

# Set output directory
set(_iface_outdir "${CMAKE_BINARY_DIR}/lib")
set_target_properties(${lib} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY  "${_iface_outdir}"
  LIBRARY_OUTPUT_DIRECTORY  "${_iface_outdir}"
  ARCHIVE_OUTPUT_DIRECTORY  "${_iface_outdir}"
  OUTPUT_NAME               "svoned_interface"
)

# Add lib prefix on Windows
if(WIN32)
  set_target_properties(${lib} PROPERTIES
    PREFIX "lib"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
  )
endif()

message(STATUS "Configuring 1D solver interface library: ${lib}")
message(STATUS "Interface library output directory: ${_iface_outdir}")
message(STATUS "Sparse solver type: ${sparseSolverType}")