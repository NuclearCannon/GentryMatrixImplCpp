cmake_minimum_required(VERSION 3.10)
project(poly_impl_flint LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/include)

# 查找src下所有cpp文件
file(GLOB_RECURSE SRC_FILES 
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*.cu"
    
)

# 查找FLINT库（假定为静态库）
find_library(FLINT_LIB flint REQUIRED)
find_library(GMP_LIB gmp REQUIRED)  # 新增：查找 GMP

# 生成静态库
add_library(poly_impl_flint STATIC ${SRC_FILES})
target_compile_options(poly_impl_flint PRIVATE -O3)
target_link_libraries(poly_impl_flint PRIVATE ${FLINT_LIB} ${GMP_LIB} profiler tcmalloc)
set_property(TARGET poly_impl_flint PROPERTY CUDA_ARCHITECTURES native)

# 编译tests下每个次级文件夹为一个可执行文件
file(GLOB TEST_SUBDIRS RELATIVE "${CMAKE_SOURCE_DIR}/tests" "${CMAKE_SOURCE_DIR}/tests/*")
foreach(subdir ${TEST_SUBDIRS})
    if(IS_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/${subdir}")
        file(GLOB TEST_SRC_FILES "${CMAKE_SOURCE_DIR}/tests/${subdir}/*.cpp")
        if(TEST_SRC_FILES) # 仅当有cpp文件时才生成可执行文件
            add_executable(${subdir} ${TEST_SRC_FILES})
            # 设置cuda架构
            set_property(TARGET ${subdir} PROPERTY CUDA_ARCHITECTURES native)
            target_include_directories(${subdir} PRIVATE ${CMAKE_SOURCE_DIR}/include)
            target_link_libraries(${subdir} PRIVATE poly_impl_flint ${FLINT_LIB} ${GMP_LIB} profiler tcmalloc)
        endif()
    endif()
endforeach()