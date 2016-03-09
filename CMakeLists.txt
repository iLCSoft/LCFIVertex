########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( LCFIVertex )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 7 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )



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
ENDIF()



### LIBRARY #################################################################


# -------------------- boost ------------------------------------------------------
if( BOOST_ROOT )
  find_package( Boost REQUIRED ) 
  include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )
  SET( LCFI_USE_EXTERNAL_BOOST 1)
else()  
  # custom command to unpack boost (if needed)
  ADD_CUSTOM_COMMAND(
    OUTPUT "${PROJECT_SOURCE_DIR}/include/boost/config.hpp"
    COMMAND tar xzvf boost.tgz
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    )
  ADD_CUSTOM_TARGET( unpack_boost DEPENDS "${PROJECT_SOURCE_DIR}/include/boost/config.hpp" )
endif()
# ---------------------------------------------------------------------------------


# definitions to pass to the compiler
#ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
#ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include ./vertex_lcfi ./vertex_lcfi/nnet/inc ./diagnostics/include )

INSTALL_DIRECTORY( ./vertex_lcfi DESTINATION include FILES_MATCHING PATTERN "*.h*" )

if(NOT BOOST_ROOT )
  INSTALL_DIRECTORY( ./include/boost DESTINATION include FILES_MATCHING PATTERN "*.h*" )
endif()

# get list of all source files
AUX_SOURCE_DIRECTORY( ./src processor_srcs )
AUX_SOURCE_DIRECTORY( ./diagnostics/src diagnostics_srcs )
FILE( GLOB_RECURSE vertex_lcfi_srcs "vertex_lcfi/*.cpp" )
SET( library_sources ${processor_srcs} ${diagnostics_srcs} ${vertex_lcfi_srcs} )


# add library
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
if(NOT BOOST_ROOT )
  ADD_DEPENDENCIES( ${PROJECT_NAME} unpack_boost ) # unpack boost (if needed)
endif()

INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake )
