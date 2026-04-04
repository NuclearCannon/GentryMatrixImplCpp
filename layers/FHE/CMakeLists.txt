file(GLOB_RECURSE FHE_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/FHE/*.cpp"
)
add_library(FHE STATIC ${FHE_SRC})
target_link_libraries(FHE PRIVATE polynomial utils)
target_include_directories(FHE PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_include_directories(FHE SYSTEM PRIVATE 
    ${CMAKE_SOURCE_DIR}/layers/GPU/include
    ${CMAKE_SOURCE_DIR}/layers/polynomial/include
    ${CMAKE_SOURCE_DIR}/layers/operators/include
)
