# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED PCRE_DIR AND NOT EXISTS ${PCRE_DIR})
  message(FATAL_ERROR "PCRE_DIR variable is defined but corresponds to non-existing directory")
endif()

# Set dependency list
set(PCRE_DEPENDENCIES "")

if(NOT PCRE_DIR)

  # PCRE (Perl Compatible Regular Expressions)
  SET(PCRE_TARGET_VERSION 8.12)
  SET(PCRE_DOWNLOAD_SOURCE_HASH "fa69e4c5d8971544acd71d1f10d59193")

  # Follow the standard EP_PREFIX locations
  SET(PCRE_PREFIX ${CMAKE_BINARY_DIR}/external)
  SET(PCRE_DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/external/PCRE-prefix)
  SET(PCRE_BINARY_DIR ${CMAKE_BINARY_DIR}/external/PCRE-prefix/src/PCRE-build)
  SET(PCRE_SOURCE_DIR ${CMAKE_BINARY_DIR}/external/PCRE-prefix/src/PCRE)
  SET(PCRE_INSTALL_DIR ${CMAKE_BINARY_DIR}/external/PCRE)

  CONFIGURE_FILE(
    ${CMAKE_MODULE_PATH}/pcre_configure_step.cmake.in
    ${CMAKE_BINARY_DIR}/external/pcre_configure_step.cmake
    @ONLY)

  SET( pcre_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/external/pcre_configure_step.cmake )

  ExternalProject_add(PCRE
    URL http://midas3.kitware.com/midas/api/rest?method=midas.bitstream.download&checksum=${PCRE_DOWNLOAD_SOURCE_HASH}&name=pcre-${PCRE_TARGET_VERSION}.tar.gz
    URL_MD5 "${PCRE_DOWNLOAD_SOURCE_HASH}"
    CONFIGURE_COMMAND ${pcre_CONFIGURE_COMMAND}
    DEPENDS "${PCRE_DEPENDENCIES}"
    PREFIX ${PCRE_PREFIX}
    DOWNLOAD_DIR ${PCRE_DOWNLOAD_DIR}
    BINARY_DIR ${PCRE_BINARY_DIR}
    SOURCE_DIR ${PCRE_SOURCE_DIR}
    INSTALL_DIR ${PCRE_INSTALL_DIR}
  )
endif(NOT PCRE_DIR)