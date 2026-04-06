file(GLOB_RECURSE APPLICATION_SRC 
    "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp"
)
add_library(application STATIC ${APPLICATION_SRC})
target_link_libraries(application PRIVATE polynomial utils)
target_include_directories(application PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_include_directories(application SYSTEM PRIVATE 
    ${CMAKE_SOURCE_DIR}/layers/complex_matrix/include
    ${CMAKE_SOURCE_DIR}/layers/GPU/include      
    ${CMAKE_SOURCE_DIR}/layers/FHE/include      
    ${CMAKE_SOURCE_DIR}/layers/polynomial/include
    ${CMAKE_SOURCE_DIR}/layers/operators/include
)
# TODO: 我们不应该直接依赖GPU，而是应该依赖于GentryPoly，而它再依赖于GPU