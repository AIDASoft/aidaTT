#!/bin/bash


if [ "$( cat /etc/*-release | grep Scientific )" ]; then
    OS=slc6
elif [ "$( cat /etc/*-release | grep CentOS )" ]; then
    OS=centos7
else
    echo "UNKNOWN OS"
    exit 1
fi


if [ -z ${GCC_VERSION} ]; then
    GCC_VERSION=4.9.3
fi
  
if [ -z ${BUILD_TYPE} ]; then
    BUILD_TYPE=opt
fi

GCC_VER=`echo ${GCC_VERSION} | sed -e 's/\.//g' | cut -c 1-2`

# General variables
CLICREPO=/cvmfs/clicdp.cern.ch
BUILD_FLAVOUR=x86_64-${OS}-gcc${GCC_VER}-${BUILD_TYPE}

#--------------------------------------------------------------------------------
#     Ninja
#--------------------------------------------------------------------------------

export Ninja_HOME=${CLICREPO}/software/Ninja/1.7.1/${BUILD_FLAVOUR}
export PATH="$Ninja_HOME:$PATH"

#--------------------------------------------------------------------------------
#     Eigen
#--------------------------------------------------------------------------------
#export EIGEN_DIR="/cvmfs/clicdp.cern.ch/software/Eigen/3.2.9/x86_64-slc6-gcc48-opt/"

#--------------------------------------------------------------------------------
#     ILCSOFT
#--------------------------------------------------------------------------------

ILCSOFT=/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/x86_64-slc6-gcc48-opt

source $ILCSOFT/init_ilcsoft.sh

export EIGEN_DIR=$EIGEN_INCLUDE_DIR
