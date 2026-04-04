file(GLOB_RECURSE COMPLEX_MATRIX_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp"
)
add_library(complex_matrix STATIC ${COMPLEX_MATRIX_SRC})
target_link_libraries(complex_matrix PRIVATE utils polynomial GPU)
target_include_directories(complex_matrix PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
# target_include_directories(complex_matrix SYSTEM PRIVATE 
#     ${CMAKE_SOURCE_DIR}/layers/GPU/include
#     ${CMAKE_SOURCE_DIR}/layers/operators/include
# )
set_property(TARGET complex_matrix PROPERTY CUDA_ARCHITECTURES native)
