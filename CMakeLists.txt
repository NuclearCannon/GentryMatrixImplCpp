cmake_minimum_required(VERSION 3.10)
project(poly_impl_flint LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

# 查找src下所有cpp文件
file(GLOB_RECURSE SRC_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# 查找FLINT库（假定为静态库）
find_library(FLINT_LIB flint REQUIRED)
# find_library(FLINT_QR_LIB flintq REQUIRED) # 若有其他依赖可继续添加

# 生成静态库
add_library(poly_impl_flint STATIC ${SRC_FILES})
target_link_libraries(poly_impl_flint PRIVATE ${FLINT_LIB} ${FLINT_QR_LIB})

# 查找tests下所有cpp文件
file(GLOB TEST_FILES "${CMAKE_SOURCE_DIR}/tests/*.cpp")

foreach(test_src ${TEST_FILES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${test_name} PRIVATE poly_impl_flint ${FLINT_LIB} ${FLINT_QR_LIB})
endforeach()
