
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)


message([STATUS] "Running cudarecv/CMakeLists.txt makefile")
message([STATUS] "Build for AARCH64 architecture for NVIDIA Jetson. Project settings are forced for AARCH64 \
in this cmake script.")
message([STATUS] "Note: This cmake script uses static filepaths since GGRG Jetsons are configured \
to specification per the setup instructions in the README.")
message([STATUS] "The authors welcome pull requests with improvements to this script for \
compiling to different platforms and dynamic filepaths.")



# Set broad specifications of the architecture/OS
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)
set(CMAKE_C_COMPILER "aarch64-linux-gnu-gcc" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "aarch64-linux-gnu-g++" CACHE FILEPATH "" FORCE)

set(CUDA_TOOLKIT_TARGET_NAME "aarch64-linux-gnu")	
set(CUDA_HOST_COMPILER "aarch64-linux-gnu-gcc")

set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/usr/lib/gcc/aarch64-linux-gnu/5;/usr/lib/aarch64-linux-gnu;/usr/lib;/lib/aarch64-linux-gnu;/lib")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "aarch64-linux-gnu")
	
	
# Specify project
project(cudarecv)


# Generate build commands for clang autocomplete
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Get CUDA Packages
unset(CUDA_rt_LIBRARY:FILEPATH CACHE)
set(CUDA_rt_LIBRARY "/usr/aarch64-linux-gnu/lib/librt.so" CACHE FILEPATH "" FORCE)

find_package(CUDA REQUIRED)
get_filename_component(CUDA_LIBRARY_DIR ${CUDA_CUDART_LIBRARY} DIRECTORY)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-L${CUDA_LIBRARY_DIR}")

find_package(CUDA REQUIRED)
unset(CUDA_rt_LIBRARY:FILEPATH CACHE)
set(CUDA_rt_LIBRARY "/usr/aarch64-linux-gnu/lib/librt.so" CACHE FILEPATH "" FORCE)

get_filename_component(CUDA_LIBRARY_DIR ${CUDA_CUDART_LIBRARY} DIRECTORY)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-L${CUDA_LIBRARY_DIR}")


# Get Eigen library
find_package(Eigen3 REQUIRED)
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3/" CACHE FILEPATH "" FORCE)
include_directories(${EIGEN3_INCLUDE_DIR})




# Specify compilation details
set ( TORCH_NVCC_FLAGS="-D__CUDA_NO_HALF_OPERATORS__" )
set (
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -D__CUDA_NO_HALF_OPERATORS__
    -G -g -O3 -ccbin aarch64-linux-gnu-g++
	--cudart static
	-lcublas
	
    -std=c++11
    -gencode arch=compute_62,code=[sm_62,compute_62]
    -lrt
    --expt-relaxed-constexpr
    )



# Link to subdirectories
subdirs(auxil console cudarecv dsp modules utils)
if (UHD)
    subdirs(guhd)
endif()




# Strong-arm project configuration to the aarch64 installation folders
set(CUDA_TOOLKIT_TARGET_NAME "aarch64-linux-gnu")
set(CUDA_TOOLKIT_TARGET_DIR "${CUDA_TOOLKIT_ROOT}/targets/${CUDA_TOOLKIT_TARGET_NAME}" CACHE PATH "CUDA Toolkit target location.")
SET(CUDA_TOOLKIT_ROOT_DIR ${CUDA_TOOLKIT_ROOT})
mark_as_advanced(CUDA_TOOLKIT_TARGET_DIR)
set(CUDA_TOOLKIT_TARGET_NAME "aarch64-linux-gnu")

unset(CUDA_rt_LIBRARY:FILEPATH CACHE)
set(CUDA_rt_LIBRARY "/usr/aarch64-linux-gnu/lib/librt.so" CACHE FILEPATH "" FORCE)	
set(CUDA_cudart_static_LIBRARY "/usr/local/cuda-9.0/targets/aarch64-linux/lib/libcudart_static.a" CACHE FILEPATH "" FORCE)
	
get_filename_component(CUDA_LIBRARY_DIR ${CUDA_CUDART_LIBRARY} DIRECTORY)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-L${CUDA_LIBRARY_DIR}")
		
unset(CUDA_cublas_LIBRARY:FILEPATH CACHE)
unset(CUDA_cufft_LIBRARY:FILEPATH CACHE)
set(CUDA_cublas_LIBRARY "/usr/local/cuda-9.0/targets/aarch64-linux/lib/stubs/libcublas.so" CACHE FILEPATH "" FORCE)
set(CUDA_cufft_LIBRARY "/usr/local/cuda-9.0/targets/aarch64-linux/lib/stubs/libcufft.so" CACHE FILEPATH "" FORCE)
