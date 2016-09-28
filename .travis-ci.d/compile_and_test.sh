#!/bin/bash

export GCC_VERSION="4.8.5"
export BUILD_TYPE="opt"

source /aidaTT/.travis-ci.d/init_x86_64.sh

cd /aidaTT
mkdir build
cd build
cmake -GNinja -C $ILCSOFT/ILCSoft.cmake ..
ninja

