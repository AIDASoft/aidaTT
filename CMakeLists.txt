# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 2.8.2 )

# declare the project name
PROJECT(aidaTT)

# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 2 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

include(MacroAddSharedLibrary)
include(library_versioning)

# require proper c++
ADD_DEFINITIONS( "-Wall -g -ansi -pedantic -Wno-long-long" )

#------------- set the default installation directory to be the source directory

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET( CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} CACHE PATH
    "install prefix path  - overwrite with -D CMAKE_INSTALL_PREFIX = ..."
    FORCE )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} - overwrite with -D CMAKE_INSTALL_PREFIX" )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#-------------------------------------------------------------------------------

# include directories
INCLUDE_DIRECTORIES( BEFORE ./core/include  ./util/include ./geometry/include ./fields/include ./persistency/include ./propagations/include ./fittingAlgorithms/include )

#INSTALL( DIRECTORY ./include DESTINATION . PATTERN ".svn" EXCLUDE )

AUX_SOURCE_DIRECTORY( ./core/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fields/src  library_sources )
AUX_SOURCE_DIRECTORY( ./persistency/src  library_sources )
AUX_SOURCE_DIRECTORY( ./propagations/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fittingAlgorithms/src  library_sources )
AUX_SOURCE_DIRECTORY( ./util/src  library_sources )
AUX_SOURCE_DIRECTORY( ./geometry/src  library_sources )

ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )

FIND_PACKAGE( GSL REQUIRED)
FIND_PACKAGE( GBL )
FIND_PACKAGE( DD4hep REQUIRED COMPONENTS DDRec )
FIND_PACKAGE( LCIO )

FOREACH( pkg  GSL GBL DD4hep LCIO)  
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${DD4hep_COMPONENT_LIBRARIES} )

#INSTALL( TARGETS ${PROJECT_NAME} DESTINATION lib )

# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)

    INCLUDE(GlobalSVNRevision)
    
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
ELSE()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generated!" )
ENDIF(DOXYGEN_FOUND)

IF( GBL_FOUND )
        ADD_DEFINITIONS( "-DUSE_GBL" )
ELSE()
    MESSAGE( STATUS "GBL not found, track fitting with GBL will not be available." )
ENDIF()

IF( DD4hep_FOUND )
        ADD_DEFINITIONS( "-DAIDATT_USE_DD4HEP" )
ELSE()
    MESSAGE( STATUS "DD4hep not found, the geometry decription from DD4hep will not be available." )
ENDIF()

IF( LCIO_FOUND )
        ADD_DEFINITIONS( "-DUSE_LCIO" )
ELSE()
    MESSAGE( STATUS "LCIO not found, LCIO persistency will not be available." )
ENDIF()

# add the examples directory; only contains executables
# build examples with `make examples`
add_subdirectory(examples EXCLUDE_FROM_ALL)

# add the test directory that contains all the tests
# build tests with `make tests`
add_subdirectory(test EXCLUDE_FROM_ALL)




#--- install target-------------------------------------
FILE(GLOB hfiles "*/include/*.hh")

INSTALL(FILES ${hfiles}
  DESTINATION include/${PROJECT_NAME}
  )

install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  )

# generate and install following configuration files
#GENERATE_PACKAGE_CONFIGURATION_FILES( aidaTTConfig.cmake aidaTTConfigVersion.cmake )
GENERATE_PACKAGE_CONFIGURATION_FILES( aidaTTConfig.cmake )

