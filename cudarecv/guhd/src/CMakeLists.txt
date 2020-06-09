
### Configure Compiler ########################################################
include_directories(
    ${Boost_INCLUDE_DIRS}
    ${UHD_INCLUDE_DIRS}
)
link_directories(
    ${Boost_LIBRARY_DIRS}
)

### Make the executable #######################################################
add_executable(guhd_record main.cpp guhd.cpp buffer.cpp)

target_compile_features(guhd_record PRIVATE cxx_range_for)
target_compile_features(guhd_record PRIVATE cxx_trailing_return_types)


message(STATUS "Linking against shared UHD library.")
target_link_libraries(guhd_record pthread ${UHD_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(guhd_record
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${guhd_SOURCE_DIR}"
    INTERPROCEDURAL_OPTIMIZATION True
)
