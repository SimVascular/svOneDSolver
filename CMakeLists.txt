
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# COMMON SETTINGS
PROJECT(OneDSolver)
ENABLE_LANGUAGE(C CXX)

# Get the SimVascular variables for libraries, includes, etc.
if (BUILD_SV_INSTALLER) 
    find_package(SimVascular REQUIRED)
    # Set the path to find the SimVascular CMake modules.
    set(CMAKE_MODULE_PATH "${SimVascular_CMAKE_DIR}" "${CMAKE_MODULE_PATH}")
    include(SimVascularMacros)
    include(SimVascularSystemSetup)
ENDIF()

# SET OPTIONS
OPTION(BUILD_SV_INSTALLER "Build SimVascular Installer" OFF)
OPTION(buildPy "Build Python Interface" OFF)
OPTION(buildDocs "Build Documentation" OFF)
SET(sparseSolverType "skyline" CACHE STRING "Use Sparse Solver")
SET_PROPERTY(CACHE sparseSolverType PROPERTY STRINGS skyline superlu csparse)

# FIX ERROR WITH BOOST
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

# ADD DEFINITION SO THE C++ CODE CAN SEE IT
IF(sparseSolverType STREQUAL "skyline")
  ADD_DEFINITIONS("-DUSE_SKYLINE")
ENDIF()
IF(sparseSolverType STREQUAL "superlu")  
  ADD_DEFINITIONS("-DUSE_SUPERLU")
ENDIF()
IF(sparseSolverType STREQUAL "csparse")
  ADD_DEFINITIONS("-DUSE_CSPARSE")
ENDIF()

# ASK THE USER TO ENTER THE SUPERLU FOLDER
IF(sparseSolverType STREQUAL "superlu")
  SET(SUPERLU_DIR " " CACHE PATH "Enter SuperLU library folder")
ENDIF()

# SET FOLDERS WITH CMAKE MODULES
SET (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# DEFINE FOLDERS VARIABLES
SET(SRCS_DIR "./Code/Source")
# ADD SUBDIRECTORY FOR DOCUMENTATION
SET(DOCS_DIR "./docs")
# ADD SUBDIRECTORY FOR PYTHON LIBRARIES
SET(PY_DIR "./py")

# INCLUDE SOURCE AND HEADER FILES
IF(sparseSolverType STREQUAL "skyline")
	FILE(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.c"
        	        "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.cxx")
	FILE(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.h")

ELSEIF(sparseSolverType STREQUAL "superlu")
	FILE(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.c"
	                "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.cxx"
     		        "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.cxx"
	                "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/superlu/*.c"
	                "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/superlu/*.cxx")
	FILE(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.h"
        	        "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.h"
                	"${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/superlu/*.h")

ELSEIF(sparseSolverType STREQUAL "csparse")
	FILE(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.c"
        	        "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.cxx"
			"${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.cxx"
			"${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/csparse/*.c"
			"${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/csparse/*.cxx")
	FILE(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/*.h"
		        "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/*.h"
			"${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/sparse/csparse/*.h")
ENDIF()

# COPY DATASET FOLDER
FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# COMPILER FLAGS
SET(CMAKE_CXX_FLAGS "-g -m64 -O3 -std=c++0x -fPIC")

# PLACE EXECUTABLE IN BIN FOLDER
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# BUILD C++ CODE
ADD_SUBDIRECTORY(${SRCS_DIR})

# BUILD DOCUMENTATION
IF(buildDocs)
  ADD_SUBDIRECTORY(${DOCS_DIR})
ENDIF()

# INSTALL PYTHON INTERFACE
IF(buildPy)

  # SWIG INTERFACE FILES
  FILE(GLOB ONED_SRC_I "${CMAKE_CURRENT_SOURCE_DIR}/Code/Source/oneDSolver.i")

  OPTION(USE_SYSTEM_SWIG "Use a pre-compiled version of SWIG 2.0 previously configured for your system" ON)

  # INSTALL SWIG IF REQUIRED
  IF( NOT USE_SYSTEM_SWIG )
    INCLUDE(installSwig)
    LIST(APPEND ${CMAKE_PROJECT_NAME}_DEPENDENCIES Swig)
  ENDIF()

  # USE SWIG TO GENERATE THE PYTHON LIBRARIES
  ADD_SUBDIRECTORY(${PY_DIR})

ENDIF()

if (BUILD_SV_INSTALLER) 
    add_subdirectory("${CMAKE_SOURCE_DIR}/Distribution")
ENDIF()
