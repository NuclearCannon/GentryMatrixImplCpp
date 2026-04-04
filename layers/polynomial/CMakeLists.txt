file(GLOB_RECURSE POLYNOMIAL_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/GentryPoly/*.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/GentryPoly/*.cu"
)
add_library(polynomial STATIC ${POLYNOMIAL_SRC})
target_link_libraries(polynomial PRIVATE utils operators complex_matrix ${FLINT_LIB} ${GMP_LIB})
target_include_directories(polynomial PRIVATE 
    ${CMAKE_CURRENT_LIST_DIR}/src
    ${CMAKE_CURRENT_LIST_DIR}/include
)
target_include_directories(polynomial SYSTEM PRIVATE 
    ${CMAKE_SOURCE_DIR}/layers/utils/include
    ${CMAKE_SOURCE_DIR}/layers/GPU/include
    ${CMAKE_SOURCE_DIR}/layers/operators/include
    ${CMAKE_SOURCE_DIR}/layers/complex_matrix/include
)
set_property(TARGET polynomial PROPERTY CUDA_ARCHITECTURES native)
