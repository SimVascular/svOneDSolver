
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

# COMMON SETTINGS
PROJECT(OneDSolver)
ENABLE_LANGUAGE(C CXX)

# SET OPTIONS
OPTION(BUILD_SV_INSTALLER "Build SimVascular Installer" OFF)
OPTION(buildPy "Build Python Interface" OFF)
OPTION(buildDocs "Build Documentation" OFF)
OPTION(ENABLE_UNIT_TEST "Enable unit tests" ON)
SET(sparseSolverType "skyline" CACHE STRING "Use Sparse Solver")
SET_PROPERTY(CACHE sparseSolverType PROPERTY STRINGS skyline superlu csparse)

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

if (BUILD_SV_INSTALLER) 
    include(SimVascularMacros)
    include(SimVascularSystemSetup)
ENDIF()

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
# FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# COMPILER FLAGS
SET(CMAKE_CXX_FLAGS "-g -m64 -O3 -std=c++20 -fPIC")

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

# Unit tests and Google Test
if(ENABLE_UNIT_TEST)

  # Link pthread on ubuntu20
  find_package(Threads REQUIRED)

  # Install Google Test
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/heads/main.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  include(GoogleTest)

  set(TESTBIN_DIRECTORY ${CMAKE_BINARY_DIR}/testbin)

  # Copy test files into the build directory so tests can access them
  file(GLOB TEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Tests/UnitTests/TestFiles/*")
  foreach(TEST_FILE ${TEST_FILES})
    # Print the file being copied for debugging
    file(COPY ${TEST_FILE} DESTINATION ${TESTBIN_DIRECTORY}/TestFiles/)
  endforeach()

  # For now, we're filtering out "main" so we can include all the 
  # other source files in the same way as the main executable.
  #
  # It would probably be better if we refactored to avoid including
  # main in the overall list of source files. Instead, we could just 
  # include main only when we're building the executable. That could 
  # be a good first refactoring to avoid this hacky post-processing 
  # step.
  set(SRC_C_FOR_TESTS ${SRC_C})
  list(FILTER SRC_C_FOR_TESTS EXCLUDE REGEX "Source/main.cxx")

  set(SRC_H_FOR_TESTS ${SRC_H})
  list(FILTER SRC_H_FOR_TESTS EXCLUDE REGEX "Source/main.h")

  # Find all .cpp files in the UnitTests folder and subfolders
  file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Tests/UnitTests/*.cxx")

  # Combine all test files into a single list
  add_executable(UnitTestsExecutable ${TEST_SOURCES} ${SRC_C_FOR_TESTS} ${SRC_H_FOR_TESTS})

  # For some solvers, we may need to link against additional libraries here.
  # When we bump up against unit tests that require that, we should consider
  # changing the strategy for how we generate executables in cmakelists to
  # handle both use cases.
  target_link_libraries(UnitTestsExecutable gtest gtest_main pthread)

  # It's likely we'll need to include additional directories for some
  # versions of the unit tests. That can be updated when/if we update
  # the general strategy.
  target_include_directories(UnitTestsExecutable PRIVATE ${SRCS_DIR})

  # Set the output directory for the executable
  set_target_properties(UnitTestsExecutable PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${TESTBIN_DIRECTORY}
  )

  # Add the test to CTest, we're setting the working directory to testbin
  # because we're copying the input files there for now
  add_test(NAME UnitTests COMMAND UnitTestsExecutable WORKING_DIRECTORY ${TESTBIN_DIRECTORY})

endif()
