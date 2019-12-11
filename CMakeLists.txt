########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( LCFIVertex )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 8 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM ${LCIO_INCLUDE_DIRS} )
ADD_DEFINITIONS( ${LCIO_DEFINITIONS} )

FIND_PACKAGE( MarlinUtil REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( SYSTEM ${MarlinUtil_INCLUDE_DIRS} )


# optional package
FIND_PACKAGE( AIDA )
IF( AIDA_FOUND )
    INCLUDE_DIRECTORIES( SYSTEM ${AIDA_INCLUDE_DIRS} )
    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
    MESSAGE( STATUS "AIDA: ${AIDA_DIR}" )
ELSE()
    MESSAGE( STATUS "AIDA -- not found" )
ENDIF()



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################


# -------------------- boost ------------------------------------------------------
find_package( Boost REQUIRED )
include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )
SET( LCFI_USE_EXTERNAL_BOOST 1)
# ---------------------------------------------------------------------------------


# definitions to pass to the compiler
#ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include ./vertex_lcfi ./vertex_lcfi/nnet/inc ./diagnostics/include )

INSTALL_DIRECTORY( ./vertex_lcfi DESTINATION include FILES_MATCHING PATTERN "*.h*" )

# get list of all source files
AUX_SOURCE_DIRECTORY( ./src processor_srcs )
AUX_SOURCE_DIRECTORY( ./diagnostics/src diagnostics_srcs )
FILE( GLOB_RECURSE vertex_lcfi_srcs "vertex_lcfi/*.cpp" )


SET( library_sources ${vertex_lcfi_srcs} )
# add library
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )

ADD_SHARED_LIBRARY( ${PROJECT_NAME}Processors ${processor_srcs} ${diagnostics_srcs} )

TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${MarlinUtil_LIBRARIES} ${LCIO_LIBRARIES} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}Processors ${Marlin_LIBRARIES} ${AIDA_LIBRARIES} ${PROJECT_NAME} )

INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME}Processors DESTINATION lib )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake  )
