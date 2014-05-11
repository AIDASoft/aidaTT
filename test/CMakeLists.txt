ADD_CUSTOM_TARGET(tests)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/core/include)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/fields/include)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/helpers/include)

INCLUDE_DIRECTORIES( BEFORE ./unitTests)

AUX_SOURCE_DIRECTORY( ./unitTests _sources )

ADD_LIBRARY (unitTests ${_sources}) 

ADD_EXECUTABLE (aidaTT-unitTests aidaTT-unitTests.cpp ${_sources})

TARGET_LINK_LIBRARIES (aidaTT-unitTests ${PROJECT_NAME}) 

ADD_DEPENDENCIES(tests aidaTT-unitTests ${_sources})
