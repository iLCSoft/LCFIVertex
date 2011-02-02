########################################################
# cmake file for building LCFI documentation
# @author Jan Engels, DESY
########################################################


FIND_PACKAGE( Doxygen )
FIND_PACKAGE( LATEX )

SET( DOC_SRC_DIR "${PROJECT_SOURCE_DIR}/doc" )
SET( DOC_BIN_DIR "${PROJECT_BINARY_DIR}/docbuild" )

# global doc target for calling doc_html and doc_latex
ADD_CUSTOM_TARGET( doc )

#doxygen
IF( DOXYGEN_FOUND )

    #SET( DOX_INPUT ../include ../src ) #TODO
    
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_SRC_DIR}/html/index.html"
        COMMAND DOX_PROJECT_NAME=${PROJECT_NAME}
                DOX_PROJECT_NUMBER="${${PROJECT_NAME}_VERSION}"
                #DOX_OUTPUT_DIRECTORY="${DOC_BIN_DIR}" # TODO
                #DOX_INPUT="${DOX_INPUT}" #TODO
                "${DOXYGEN_EXECUTABLE}"
        WORKING_DIRECTORY "${DOC_SRC_DIR}"
        COMMENT "Building API Documentation..."
        DEPENDS Doxyfile CMakeLists.txt
    )
    
    # add doc_html target
    ADD_CUSTOM_TARGET( doc_html DEPENDS "${DOC_SRC_DIR}/html/index.html" )
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
        OUTPUT  "${DOC_SRC_DIR}/latex/index.tex"
        COMMAND "${DOXYGEN_EXECUTABLE}" ARGS Doxyfile.latex
        WORKING_DIRECTORY "${DOC_SRC_DIR}"
        COMMENT "Generating Latex Documentation with Doxygen..."
    )
    
    # run latex + pdflatex
    ADD_CUSTOM_COMMAND(
        OUTPUT  "${DOC_SRC_DIR}/latex/refman.pdf"
        COMMAND "${CMAKE_COMMAND}" -E copy "${DOC_SRC_DIR}/VtxPackageDirStruct.jpg" "${DOC_SRC_DIR}/latex"
		COMMAND "${CMAKE_COMMAND}" -E copy "${DOC_SRC_DIR}/VertexPackageExeFlow.jpg" "${DOC_SRC_DIR}/latex"
        COMMAND "${CMAKE_BUILD_TOOL}"
        # copy pdf to doc directory
        COMMAND "${CMAKE_COMMAND}" -E copy "refman.pdf" "${DOC_SRC_DIR}"
        DEPENDS "${DOC_SRC_DIR}/latex/index.tex"
        WORKING_DIRECTORY "${DOC_SRC_DIR}/latex"
        COMMENT "Building Latex Documentation..."
    )
    
    # add doc_latex target
    ADD_CUSTOM_TARGET( doc_latex DEPENDS "${DOC_SRC_DIR}/latex/refman.pdf" )
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

    # in-source vs out-of-source installations
    IF( CMAKE_INSTALL_PREFIX STREQUAL "${PROJECT_SOURCE_DIR}" )

        # --- in-source installations ---
        SET( DOC_INSTALL_PREFIX doc )

        ## keep source tree clean
        ADD_CUSTOM_TARGET( uninstall_doc
            COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_INSTALL_PREFIX}/doc/html"
            COMMAND "${CMAKE_COMMAND}" -E remove_directory "${CMAKE_INSTALL_PREFIX}/doc/latex"
            COMMAND "${CMAKE_COMMAND}" -E remove "${CMAKE_INSTALL_PREFIX}/doc/refman.pdf"
        )
        ADD_DEPENDENCIES( uninstall uninstall_doc )

    ELSE()

        # --- out-of-source installations ---
        SET( DOC_INSTALL_PREFIX doc/${PROJECT_NAME} )

        INSTALL( FILES refman.pdf zvtop.pdf DESTINATION ${DOC_INSTALL_PREFIX} )

    ENDIF( CMAKE_INSTALL_PREFIX STREQUAL "${PROJECT_SOURCE_DIR}" )

    # install html documentation
    INSTALL_DIRECTORY( "${DOC_SRC_DIR}/html" DESTINATION ${DOC_INSTALL_PREFIX} )

ENDIF()

