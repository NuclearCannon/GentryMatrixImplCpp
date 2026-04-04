file(GLOB OPERATORS_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/ntt/*.cpp"
)
add_library(operators STATIC ${OPERATORS_SRC})
target_link_libraries(operators PRIVATE utils GPU ${FLINT_LIB} ${GMP_LIB})
target_include_directories(operators PRIVATE 
    ${CMAKE_CURRENT_LIST_DIR}/src
    ${CMAKE_CURRENT_LIST_DIR}/include
)
target_include_directories(operators SYSTEM PRIVATE 
    ${CMAKE_SOURCE_DIR}/layers/utils/include
    ${CMAKE_SOURCE_DIR}/layers/GPU/include
)
