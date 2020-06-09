
include_directories(${cudarecv_SOURCE_DIR}/auxil/inc)
cuda_add_library(auxil auxil.cpp eigenwrapper.cpp tm_usage.cpp maybe.cpp)

target_compile_features(auxil PRIVATE cxx_variadic_templates)
