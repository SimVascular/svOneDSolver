# Set cpack variables to generate the solver installer.
#
# The installers for the different platforms are
#
#   Linux - The solver is packaged in a DEB file using the 'deb' cpack generator.
#   MacOS - The solver is packaged in a PKG file using the 'productbuild' cpack generator.
#
# The solver install directory structure is defined by the variables
#
#    SV_SOLVER_INSTALL_DIR/sv1dsolver/SOLVER_INSTALL_VERSION_DIR
#
# These variables are derived from variables set in SimVascular CMake files
#
#     SV_SOLVER_INSTALL_DIR - The plugins install directory, set in SimVascular/Code/CMake/SimVascularSystemSetup.cmake.
#                             This includes the main install directory and SimVascular version.
#
include(CPackComponent)

## Set general package information.
#
set(SOLVER_NAME "oneDSolver")
set(CPACK_PACKAGE_NAME "oneD-solver")
set(CPACK_PACKAGE_VENDOR "SimVascular.org")
set(CPACK_PACKAGE_CONTACT "SimVascular<sv@simvascular.org>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "1D Solver")
set(SOLVER_INSTALL_SCRIPT_DIR ".")

## Set the solver installation directory.
if (UNIX)
  set(SV_SOLVER_INSTALL_DIR "/usr/local/sv/")
endif()

string(REPLACE "." "-" SV_SOLVER_INSTALL_VERSION_DIR ${SimVascular_VERSION})

if (UNIX)
  set(CPACK_PACKAGING_INSTALL_PREFIX "${SV_SOLVER_INSTALL_DIR}/${SOLVER_NAME}/${SV_SOLVER_INSTALL_VERSION_DIR}")
endif()

message(STATUS "SimVascular_VERSION: ${SimVascular_VERSION}")
message(STATUS "SV_SOLVER_INSTALL_DIR: ${SV_SOLVER_INSTALL_DIR}")
message(STATUS "SV_SOLVER_INSTALL_VERSION_DIR: ${SV_SOLVER_INSTALL_VERSION_DIR}")
message(STATUS "CPACK_PACKAGING_INSTALL_PREFIX: ${CPACK_PACKAGING_INSTALL_PREFIX}")

## Set the generators used to create the install file and create a setup.sh script
#  that sets environment variables to locate the solver libraries.
#
if(APPLE)
  set(CPACK_GENERATOR "productbuild")
  #set(CPACK_COMPONENTS_ALL libraries)
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-darwin-x64-${SV_FULL_VERSION}${INSTALL_TYPE_STR}")
  set(SOLVER_SETUP_SCRIPT "setup.sh") 
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup-macos.in ${SOLVER_SETUP_SCRIPT})

elseif(LINUX)
  set(CPACK_GENERATOR DEB)
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-linux-x64-${SV_FULL_VERSION}${INSTALL_TYPE_STR}")
  set(SOLVER_SETUP_SCRIPT "setup.sh") 
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup-linux.in ${SOLVER_SETUP_SCRIPT})

# [TODO] This needs to be finished.
elseif(WIN32)
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-windows-win-x64-${SV_FULL_VERSION}${INSTALL_TYPE_STR}")
  set(CPACK_MONOLITHIC_INSTALL 1)
  set(CPACK_COMPONENTS_IGNORE_GROUPS 1)
endif()

## Install setup script into plugin directory.
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${SOLVER_SETUP_SCRIPT} DESTINATION ${SOLVER_INSTALL_SCRIPT_DIR})

## Include CPack to introduce the appropriate targets
include(CPack)

