file(GLOB_RECURSE UTILS_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp"
)
add_library(utils STATIC ${UTILS_SRC})
target_link_libraries(utils PUBLIC ${FLINT_LIB} ${GMP_LIB})
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
