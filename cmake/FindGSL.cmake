#########################################################
# cmake module for finding GSL
#########################################################

# find gsl-config
SET( GSL_CONFIG_EXECUTABLE GSL_CONFIG_EXECUTABLE-NOTFOUND )
MARK_AS_ADVANCED( GSL_CONFIG_EXECUTABLE )

# use the specified/pre-defined directory
FIND_PROGRAM( GSL_CONFIG_EXECUTABLE gsl-config PATHS ${GSL_DIR}/bin NO_DEFAULT_PATH )

# alternative: look for the executable
IF( NOT GSL_DIR )
    FIND_PROGRAM( GSL_CONFIG_EXECUTABLE gsl-config )
ENDIF()

# set the variables: GSL_PREFIX, GSL_INCLUDE_DIRS, GSL_VERSION and GSL_LIBRARIES
IF( GSL_CONFIG_EXECUTABLE )

    # set basic path for include
    EXECUTE_PROCESS( COMMAND "${GSL_CONFIG_EXECUTABLE}" --prefix
        OUTPUT_VARIABLE GSL_PREFIX
        RESULT_VARIABLE _exit_code
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    FIND_PATH( GSL_INCLUDE_DIRS NAME gsl/gsl_matrix.h PATHS ${GSL_PREFIX}/include NO_DEFAULT_PATH )

     # set version number
    EXECUTE_PROCESS( COMMAND "${GSL_CONFIG_EXECUTABLE}" --version
        OUTPUT_VARIABLE GSL_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # get the libraries
    EXEC_PROGRAM(${GSL_CONFIG_EXECUTABLE}
        ARGS --libs
        OUTPUT_VARIABLE GSL_LIBRARIES 
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
SET( GSL_FOUND TRUE)
ELSE()
    MESSAGE (STATUS "No gsl-config found. Is the gsl installed, and gsl-config in $PATH?")

ENDIF( GSL_CONFIG_EXECUTABLE )
