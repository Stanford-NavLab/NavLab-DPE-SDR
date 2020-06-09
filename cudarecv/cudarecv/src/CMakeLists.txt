include_directories(${cudarecv_SOURCE_DIR}/cudarecv/inc)

# Modules
include_directories(${cudarecv_SOURCE_DIR}/naveng/inc)
link_directories(${cudarecv_SOURCE_DIR}/naveng/src)
include_directories(${cudarecv_SOURCE_DIR}/console/inc)
link_directories(${cudarecv_SOURCE_DIR}/console/src)
include_directories(${cudarecv_SOURCE_DIR}/modules/inc)
link_directories(${cudarecv_SOURCE_DIR}/modules/src)
include_directories(${cudarecv_SOURCE_DIR}/dsp/inc)
include_directories(${cudarecv_SOURCE_DIR}/dsp/src)
link_directories(${cudarecv_SOURCE_DIR}/dsp/src)
include_directories(${cudarecv_SOURCE_DIR}/auxil/inc)
link_directories(${cudarecv_SOURCE_DIR}/auxil/src)

cuda_add_executable(cudarecv main.cu startoptions.cpp cmdCommon.cpp cmdFlow.cpp)
target_link_libraries(cudarecv
    console
    modules
    dsp
    pthread
    ${CUDA_cufft_LIBRARY}
)

target_compile_features(cudarecv PRIVATE cxx_lambdas)
target_compile_features(cudarecv PRIVATE cxx_binary_literals)
target_compile_features(cudarecv PRIVATE cxx_range_for)

set_target_properties(cudarecv
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    INTERPROCEDURAL_OPTIMIZATION True
)

