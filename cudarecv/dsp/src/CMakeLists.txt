include_directories(${cudarecv_SOURCE_DIR}/dsp/inc)
include_directories(${cudarecv_SOURCE_DIR}/dsp/src/utils)
include_directories(${cudarecv_SOURCE_DIR}/auxil/inc)
include_directories(${cudarecv_SOURCE_DIR}/utils/inc)

include_directories(${cudarecv_SOURCE_DIR}/modules/inc)
link_directories(${cudarecv_SOURCE_DIR}/modules/src)

file(GLOB dsp_internal_headers
    "*.h"
)

cuda_add_library(dsp
    flow.cu
    dpeflow.cpp
    flowmgr.cpp
    ${dsp_internal_headers})

target_link_libraries(dsp modules)

target_compile_features(dsp PRIVATE cxx_lambdas cxx_binary_literals)

