#!/bin/bash

cd ..
rm -rf build
mkdir build
cd build
export LD_LIBRARY_PATH=/usr/local/cuda-9.0/targets/aarch64-linux/lib:$LD_LIBRARY_PATH
cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_ECLIPSE_VERSION=4.4 -DCMAKE_BUILD_TYPE=Debug -DJETSON_BUILD=true ..
cd ../scripts
