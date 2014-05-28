ADD_CUSTOM_TARGET(examples)

MACRO( ADD_AIDATT_EXAMPLE _example )
   ADD_EXECUTABLE( ${_example} EXCLUDE_FROM_ALL ${ARGN} )
   ADD_DEPENDENCIES( examples  ${_example} )
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/core/include)
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/fields/include)
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/geometryAPIs/include)
   TARGET_LINK_LIBRARIES( ${_example} ${PROJECT_NAME}  )
ENDMACRO( ADD_AIDATT_EXAMPLE )

MACRO( ADD_AIDATT_EXAMPLE_WITH_ROOT _example )
   ADD_EXECUTABLE( ${_example} EXCLUDE_FROM_ALL ${ARGN} )
   ADD_DEPENDENCIES( examples  ${_example} )
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/core/include)
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/fields/include)
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/geometryAPIs/include)
   TARGET_LINK_LIBRARIES( ${_example} ${PROJECT_NAME} ${ROOT_LIBS} )
ENDMACRO( ADD_AIDATT_EXAMPLE_WITH_ROOT )

IF(DD4HEP_FOUND)
ADD_AIDATT_EXAMPLE_WITH_ROOT( geometry_test geometry_test/geometry_test.cpp)
ADD_AIDATT_EXAMPLE_WITH_ROOT( intersection_calculator intersection_calculator/intersection_calculator.cpp)
    IF(LCIO_FOUND)
    ADD_AIDATT_EXAMPLE_WITH_ROOT ( lcio_read_example  lcio_read_example/lcio_read_example.cpp )
ENDIF(LCIO_FOUND)

ENDIF(DD4HEP_FOUND)


IF( GBL_FOUND)
ADD_AIDATT_EXAMPLE_WITH_ROOT( gbl_example gbl_example/gbl_example.cpp)
ENDIF( GBL_FOUND )
