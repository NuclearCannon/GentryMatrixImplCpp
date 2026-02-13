#pragma once
#include <utility>
#include <cstddef>


class CudaBuffer {
private:
    void* ptr_;
    size_t len_;

public:
    // 给外面用的工厂函数
    // 构造：分配 device 内存（size 字节）
    explicit CudaBuffer(size_t size_bytes);
    
    // 禁用拷贝（或实现深拷贝）
    CudaBuffer(const CudaBuffer&) = delete;
    CudaBuffer& operator=(const CudaBuffer&) = delete;

    // 移动构造/赋值（推荐支持）
    CudaBuffer(CudaBuffer&& other) noexcept;
    CudaBuffer& operator=(CudaBuffer&& other) noexcept;

    // 析构：自动释放 device 内存
    ~CudaBuffer();

    // 获取 host 可读写的临时副本（可选，用于调试/输出）
    // 注意：这不是高效路径，仅用于必要场景
    void copy_to_host(void* host_ptr) const;
    void copy_from_host(const void* host_ptr);

    // 返回字节数
    size_t size() const { return len_; }

    // 仅限.cu文件使用！
    void* get_ptr() {return ptr_;}

    // 仅限.cu文件使用！
    const void* get_ptr() const {return ptr_;}

};