# Include CMake External Project Module
INCLUDE(ExternalProject)

SET ( WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external )

# Make sure this file is included only once
#GET_FILENAME_COMPONENT( CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE )
#if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
#  RETURN()
#ENDIF()

#SET(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
IF(DEFINED SWIG_DIR AND NOT EXISTS ${SWIG_DIR})
  MESSAGE(FATAL_ERROR "Swig_DIR variable is defined but corresponds to non-existing directory")
ENDIF()

IF(NOT SWIG_DIR)

  SET(SWIG_TARGET_VERSION 2.0.12-1)
  SET(SWIG_DOWNLOAD_SOURCE_HASH "44af22bffb53d1795b0f5cb3bff5eb1a")
  SET(SWIG_DOWNLOAD_WIN_HASH "439bc49355dc76490b3fe0dffac2774d")


  IF(WIN32)
    # binary SWIG for windows
    set(swig_source_dir ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${SWIG_TARGET_VERSION})

    # swig.exe available as pre-built binary on Windows:
    ExternalProject_Add(SWIG
      URL http://midas3.kitware.com/midas/api/rest?method=midas.bitstream.download&checksum=${SWIG_DOWNLOAD_WIN_HASH}&name=swigwin-${SWIG_TARGET_VERSION}.zip
      URL_MD5 ${SWIG_DOWNLOAD_WIN_HASH}
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${SWIG_TARGET_VERSION}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      )

    set(SWIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${SWIG_TARGET_VERSION}) # path specified as source in ep
    set(SWIG_EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/swigwin-${SWIG_TARGET_VERSION}/swig.exe)

  else()
    # compiled SWIG for others

    # Set dependency list
    SET(SWIG_DEPENDENCIES "PCRE")

    # Install PCRE if needed 
    # PCRE (Perl Compatible Regular Expressions)
    INCLUDE(External_PCRE)

    # Install SWIG

    # swig uses bison find it by cmake and pass it down
    FIND_PACKAGE(BISON)
    SET( BISON_FLAGS "" CACHE STRING "Flags used by bison" )
    mark_as_advanced( BISON_FLAGS )

    # follow the standard EP_PREFIX locations
    SET(SWIG_PREFIX ${CMAKE_BINARY_DIR}/external)
    SET(SWIG_BINARY_DIR ${CMAKE_BINARY_DIR}/external/Swig-prefix/src/Swig-build)
    SET(SWIG_SOURCE_DIR ${CMAKE_BINARY_DIR}/external/Swig-prefix/src/Swig)
    SET(SWIG_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/Swig)

    # configure step
    configure_file(
      ${CMAKE_MODULE_PATH}/swig_configure_step.cmake.in
      ${CMAKE_BINARY_DIR}/external/swig_configure_step.cmake
      @ONLY)
    SET(SWIG_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/external/swig_configure_step.cmake)


    ExternalProject_add(SWIG
      URL http://midas3.kitware.com/midas/api/rest?method=midas.bitstream.download&checksum=${SWIG_DOWNLOAD_SOURCE_HASH}&name=swig-${SWIG_TARGET_VERSION}.tar.gz
      URL_MD5 ${SWIG_DOWNLOAD_SOURCE_HASH}
      CONFIGURE_COMMAND ${SWIG_CONFIGURE_COMMAND}
      DEPENDS "${Swig_DEPENDENCIES}"
      PREFIX ${SWIG_PREFIX}
      BINARY_DIR ${SWIG_BINARY_DIR}
      SOURCE_DIR ${SWIG_SOURCE_DIR}
      INSTALL_DIR ${SWIG_INSTALL_DIR}
    )

    set(SWIG_DIR ${SWIG_INSTALL_DIR}/share/swig/${SWIG_TARGET_VERSION})
    set(SWIG_EXECUTABLE ${SWIG_INSTALL_DIR}/bin/swig)

  endif()
endif(NOT SWIG_DIR)