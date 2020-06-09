#!/bin/bash

wget --no-check-certificate http://bitbucket.org/eigen/eigen/get/3.3.5.tar.gz
tar zxvf 3.3.5.tar.gz
cd eigen-eigen-*
mkdir build
cd build
cmake -DCMAKE_C_COMPILER:PATH="/usr/bin/aarch64-linux-gnu-gcc" -DCMAKE_CXX_COMPILER:PATH="/usr/bin/aarch64-linux-gnu-g++" ..
sudo make install
cd ../..
rm -rf eigen-eigen-*
rm -f 3.3.5.tar.gz
#cd /usr/local/include
cd /usr/include
sudo ln -s eigen3/Eigen Eigen
# Edit Eigen so it plays nicely with CUDA
#sudo sed -e'410i#ifndef EIGEN_NO_HALF' /usr/local/include/eigen3/Eigen/Core
#sudo sed -i '414i#endif' /usr/local/include/eigen3/Eigen/Core
sudo sed -i -e '410i//CUDARecv hack to permit cross-compilation for NVIDIA Jetson architecture\n#ifndef EIGEN_NO_HALF' -e '414i#endif' -e '411,413s/^/  /' /usr/local/include/eigen3/Eigen/Core
#export EIGEN3_INCLUDE_DIR=/usr/local/include/eigen3
#echo "\n\nPlease restart your computer for the system variables to include the path to the Eigen libraries.\n"
