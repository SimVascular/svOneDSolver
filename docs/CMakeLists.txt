FIND_PACKAGE(Sphinx REQUIRED)
 
IF(NOT DEFINED SPHINX_THEME)
  SET(SPHINX_THEME default)
ENDIF()
 
IF(NOT DEFINED SPHINX_THEME_DIR)
  SET(SPHINX_THEME_DIR)
ENDIF()
 
# configured documentation tools and intermediate build results
SET(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")
 
# Sphinx cache with pickled ReST documents
SET(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/.doctrees")
 
# HTML output directory
SET(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_static")
 
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/source/conf.py.in"
    "${BINARY_BUILD_DIR}/conf.py"
    @ONLY)
 
ADD_CUSTOM_TARGET(docs
        ${SPHINX_EXECUTABLE}
        -v -b html
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/source"
        "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx")
