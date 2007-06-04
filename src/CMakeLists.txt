##################################################
# cmake file for building Marlin Example Package
# @author Jan Engels, DESY
##################################################


####################################################################################################
### INCLUDES
####################################################################################################

#include directories
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/include" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi" )
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi/nnet/inc" )

# boost
INCLUDE_DIRECTORIES( "${${PROJECT_NAME}_SOURCE_DIR}/boost" )

# install include files
INSTALL( DIRECTORY
        "${${PROJECT_NAME}_SOURCE_DIR}/include"
        "${${PROJECT_NAME}_SOURCE_DIR}/boost"
		DESTINATION .
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )

INSTALL( DIRECTORY
        "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi/inc"
        DESTINATION ./vertex_lcfi
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )

INSTALL( DIRECTORY
        "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi/algo/inc"
        DESTINATION ./vertex_lcfi/algo
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )

INSTALL( DIRECTORY
        "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi/nnet/inc"
        DESTINATION ./vertex_lcfi/nnet
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )

INSTALL( DIRECTORY
        "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi/util/inc"
        DESTINATION ./vertex_lcfi/util
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )

INSTALL( DIRECTORY
        "${${PROJECT_NAME}_SOURCE_DIR}/vertex_lcfi/zvtop/include"
        DESTINATION ./vertex_lcfi/zvtop
		PATTERN "*~" EXCLUDE
		PATTERN "*CVS*" EXCLUDE )

####################################################################################################
### SOURCES
####################################################################################################

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# get list of all source files
AUX_SOURCE_DIRECTORY( . currDir )
AUX_SOURCE_DIRECTORY( ../vertex_lcfi/src vtx_src )
AUX_SOURCE_DIRECTORY( ../vertex_lcfi/algo/src vtx_algo )
AUX_SOURCE_DIRECTORY( ../vertex_lcfi/nnet/src vtx_nnet )
AUX_SOURCE_DIRECTORY( ../vertex_lcfi/util/src vtx_util )
AUX_SOURCE_DIRECTORY( ../vertex_lcfi/zvtop/src vtx_zvtop )

SET( lib_srcs ${currDir}
        ${vtx_src} ${vtx_algo} ${vtx_nnet} ${vtx_util} ${vtx_zvtop} )


####################################################################################################
### DEPENDENCIES
####################################################################################################

# load macro
IF( DEFINED ${PROJECT_NAME}_DEPENDS OR DEFINED BUILD_WITH OR DEFINED LINK_WITH )
    # load macro
    IF( NOT EXISTS "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
        MESSAGE( FATAL_ERROR "\nSorry, could not find MacroCheckDeps.cmake...\n"
            "Please set CMAKE_MODULE_PATH correctly with: "
            "cmake -DCMAKE_MODULE_PATH=<path_to_cmake_modules>" )
    ENDIF()
    INCLUDE( "${CMAKE_MODULE_PATH}/MacroCheckDeps.cmake" )
    CHECK_DEPS()
ENDIF()


####################################################################################################
### LIBRARY
####################################################################################################

# create library
ADD_LIBRARY( lib_${PROJECT_NAME} ${lib_srcs} )

# create symbolic lib target for calling target lib_XXX
ADD_CUSTOM_TARGET( lib DEPENDS lib_${PROJECT_NAME} )

# change lib_target properties
SET_TARGET_PROPERTIES( lib_${PROJECT_NAME} PROPERTIES
	# create *nix style library versions + symbolic links
	VERSION ${${PROJECT_NAME}_VERSION}
	SOVERSION ${${PROJECT_NAME}_SOVERSION}
	# allow creating static and shared libs without conflicts
	CLEAN_DIRECT_OUTPUT 1
	# avoid conflicts between library and binary target names
	OUTPUT_NAME ${PROJECT_NAME}	)

# install library
INSTALL( TARGETS lib_${PROJECT_NAME} DESTINATION lib PERMISSIONS
		OWNER_READ OWNER_WRITE OWNER_EXECUTE
		GROUP_READ GROUP_EXECUTE
		WORLD_READ WORLD_EXECUTE )

