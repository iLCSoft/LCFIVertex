########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( LCFIVertex )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 5 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )


FIND_PACKAGE( MarlinUtil REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${MarlinUtil_INCLUDE_DIRS} )
LINK_LIBRARIES( ${MarlinUtil_LIBRARIES} )


# optional package
FIND_PACKAGE( AIDA )
IF( AIDA_FOUND )
    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
    MESSAGE( STATUS "AIDA: ${AIDA_DIR}" )
ELSE()
    MESSAGE( STATUS "AIDA -- not found" )
ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include ./vertex_lcfi ./vertex_lcfi/nnet/inc ./diagnostics/include )
INCLUDE_DIRECTORIES( ./boost )


# get list of all source files
AUX_SOURCE_DIRECTORY( ./src processor_srcs )
AUX_SOURCE_DIRECTORY( ./diagnostics/src diagnostics_srcs )
FILE( GLOB_RECURSE vertex_lcfi_srcs "vertex_lcfi/*.cpp" )
SET( library_sources ${processor_srcs} ${diagnostics_srcs} ${vertex_lcfi_srcs} )

# add library
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

