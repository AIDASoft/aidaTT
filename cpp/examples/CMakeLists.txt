ADD_CUSTOM_TARGET(examples)

MACRO( ADD_AIDATT_EXAMPLE _example )
   ADD_EXECUTABLE( ${_example} EXCLUDE_FROM_ALL ${ARGN} )
   ADD_DEPENDENCIES( examples  ${_example} )
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/core/include)
   INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/fields/include)
   TARGET_LINK_LIBRARIES( ${_example} ${PROJECT_NAME}  )
ENDMACRO( ADD_AIDATT_EXAMPLE )

ADD_AIDATT_EXAMPLE( simpleton simple/simpleton.cpp)

IF( GBL_FOUND)
ADD_AIDATT_EXAMPLE( gbl_example gbl_example/gbl_example.cpp)
ENDIF( GBL_FOUND )
