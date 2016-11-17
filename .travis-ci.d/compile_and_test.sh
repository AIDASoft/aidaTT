#!/bin/bash


ILCSOFT=/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_gcc
source $ILCSOFT/init_ilcsoft.sh

cd /aidaTT
mkdir build
cd build
cmake -GNinja -C $ILCSOFT/ILCSoft.cmake ..
ninja

