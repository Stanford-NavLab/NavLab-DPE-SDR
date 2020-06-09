include_directories(${cudarecv_SOURCE_DIR}/auxil/inc)
include_directories(${cudarecv_SOURCE_DIR}/auxil/src)
include_directories(${cudarecv_SOURCE_DIR}/console/inc)
include_directories(${cudarecv_SOURCE_DIR}/dsp/inc)
include_directories(${cudarecv_SOURCE_DIR}/modules/inc)
include_directories(${cudarecv_SOURCE_DIR}/utils/inc)

set(DSP_MODULES_SRC
    module.cpp
    sampleblock.cu
    batchcorrscores.cu
    batchcorrmanifold.cu
    dpinit.cpp
    cuekf.cu
    cuchanmgr.cu
    satpos.cu
    datalogger.cu
    )

cuda_add_library(modules ${DSP_MODULES_SRC})

if (NAV_BENCH)
    add_definitions(-DNAVBENCH)
endif()

target_link_libraries(modules console auxil utils ${CUDA_cufft_LIBRARY} ${CUDA_cublas_LIBRARY})
target_compile_features(modules PRIVATE cxx_lambdas cxx_variable_templates)

if (JETSON_BUILD)
	message([STATUS] "[Module] Building for Jetson")
endif (JETSON_BUILD)
