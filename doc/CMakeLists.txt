########################################################
# cmake file for building LCFI documentation
# @author Jan Engels, DESY
########################################################


FIND_PACKAGE( Doxygen )
FIND_PACKAGE( LATEX )

# global doc target for calling doc_html and doc_latex
ADD_CUSTOM_TARGET( doc )

#doxygen
IF( DOXYGEN_FOUND )
    
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/html"
        COMMAND "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Building API Documentation..."
    )
    
    # add doc_html target
    ADD_CUSTOM_TARGET( doc_html DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/html" )
    # tell doc target to call doc_html
    ADD_DEPENDENCIES( doc doc_html )
    
ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

#latex
IF( DOXYGEN_FOUND AND LATEX_COMPILER AND
    MAKEINDEX_COMPILER AND PDFLATEX_COMPILER )
     
     # generate latex with doxygen
     ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/latex"
        COMMAND "${DOXYGEN_EXECUTABLE}" ARGS Doxyfile.latex
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Generating Latex Documentation with Doxygen..."
        VERBATIM )
    
    # run latex + pdflatex
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${PROJECT_SOURCE_DIR}/doc/latex/refman.pdf"
        COMMAND "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/doc/VtxPackageDirStruct.jpg" "${PROJECT_SOURCE_DIR}/doc/latex"
				COMMAND "${CMAKE_COMMAND}" -E copy "${PROJECT_SOURCE_DIR}/doc/VertexPackageExeFlow.jpg" "${PROJECT_SOURCE_DIR}/doc/latex"
        COMMAND "${CMAKE_BUILD_TOOL}"
        # copy pdf to doc directory
        COMMAND "${CMAKE_COMMAND}" -E copy "refman.pdf" "${PROJECT_SOURCE_DIR}/doc"
        DEPENDS "${PROJECT_SOURCE_DIR}/doc/latex"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc/latex"
        COMMENT "Building Latex Documentation..."
        VERBATIM )
    
    # add doc_latex target
    ADD_CUSTOM_TARGET( doc_latex DEPENDS
        "${PROJECT_SOURCE_DIR}/doc/latex/refman.pdf" )
    # tell doc target to call doc_latex
    ADD_DEPENDENCIES( doc doc_latex )
    
ELSE()
    MESSAGE( STATUS "Latex not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# install documentation
IF( INSTALL_DOC )
    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXECUTE_PROCESS( COMMAND ${CMAKE_BUILD_TOOL} doc)" )

    # install documentation
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
            DESTINATION .
            PATTERN "*CVS*" EXCLUDE
            PATTERN "*/.svn*" EXCLUDE
    )
ENDIF()

