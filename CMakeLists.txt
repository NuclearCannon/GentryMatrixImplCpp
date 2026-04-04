cmake_minimum_required(VERSION 3.10)
project(poly_impl_flint LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)

find_library(FLINT_LIB flint REQUIRED)
find_library(GMP_LIB gmp REQUIRED)

add_subdirectory(layers/utils)
add_subdirectory(layers/GPU)
add_subdirectory(layers/operators)
add_subdirectory(layers/polynomial)
add_subdirectory(layers/FHE)
add_subdirectory(layers/complex_matrix)

add_library(poly_impl_flint INTERFACE)
target_link_libraries(poly_impl_flint INTERFACE 
    utils GPU operators polynomial FHE complex_matrix 
    ${FLINT_LIB} ${GMP_LIB} profiler tcmalloc
)

file(GLOB TEST_SUBDIRS RELATIVE "${CMAKE_SOURCE_DIR}/tests" "${CMAKE_SOURCE_DIR}/tests/*")
foreach(subdir ${TEST_SUBDIRS})
    if(IS_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/${subdir}")
        file(GLOB TEST_SRC_FILES "${CMAKE_SOURCE_DIR}/tests/${subdir}/*.cpp")
        if(TEST_SRC_FILES)
            add_executable(${subdir} ${TEST_SRC_FILES})
            set_property(TARGET ${subdir} PROPERTY CUDA_ARCHITECTURES native)
            target_link_libraries(${subdir} PRIVATE poly_impl_flint ${FLINT_LIB} ${GMP_LIB} profiler tcmalloc)
            target_include_directories(${subdir} PRIVATE 
                "${CMAKE_SOURCE_DIR}/layers/utils/include"
                "${CMAKE_SOURCE_DIR}/layers/GPU/include"
                "${CMAKE_SOURCE_DIR}/layers/operators/include"
                "${CMAKE_SOURCE_DIR}/layers/polynomial/include"
                "${CMAKE_SOURCE_DIR}/layers/FHE/include"
                "${CMAKE_SOURCE_DIR}/layers/complex_matrix/include"
            )
        endif()
    endif()
endforeach()
