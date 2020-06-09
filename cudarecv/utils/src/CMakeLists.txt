include_directories(${cudarecv_SOURCE_DIR}/utils/inc)
include_directories(${cudarecv_SOURCE_DIR}/auxil/inc)

cuda_add_library(utils converters.cpp rinexparse.cpp)

target_link_libraries(utils auxil)
target_compile_features(utils PRIVATE cxx_variadic_templates)
