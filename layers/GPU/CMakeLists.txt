file(GLOB_RECURSE GPU_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/*.cu"
)
add_library(GPU STATIC ${GPU_SRC})
target_include_directories(GPU PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
set_property(TARGET GPU PROPERTY CUDA_ARCHITECTURES native)
target_include_directories(GPU PRIVATE "/usr/local/cuda-13/include")
target_link_libraries(GPU PRIVATE utils ${FLINT_LIB} ${GMP_LIB})
