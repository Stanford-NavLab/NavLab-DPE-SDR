
include_directories(${cudarecv_SOURCE_DIR}/console/inc)

include_directories(${cudarecv_SOURCE_DIR}/auxil/inc)
#link_directories(${cudarecv_SOURCE_DIR}/auxil/src)

add_library(console console.cpp cmdDofile.cpp cmdExec.cpp cmdHistory.cpp cmdParser.cpp cmdReader.cpp sysmsg.cpp)

target_link_libraries(console auxil)

target_compile_features(console PRIVATE cxx_variadic_templates)
