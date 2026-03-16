
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_check.hpp"
#include <assert.h>
#include <mutex>
#include <unordered_map>
#include <vector>

class SimpleCudaMemoryPool {
private:
    // 核心数据结构：size -> 空闲指针栈
    std::unordered_map<size_t, std::vector<void*>> freeLists;
    
    // 锁：如果是单线程应用，可以去掉锁以进一步提升性能
    std::mutex poolMutex;

    // 配置参数
    size_t alignment;          // 对齐字节数，例如 32 或 512
    size_t maxCachePerSize;    // 每个大小类别最多缓存多少个块，防止显存膨胀

public:
    /**
     * @param align 内存对齐字节数 (建议 32 或 64，避免 warp 分支效率问题)
     * @param maxCache 每个尺寸类别的最大缓存数量 (超过此数量会真正 free)
     */
    SimpleCudaMemoryPool(size_t align = 32, size_t maxCache = 10) 
        : alignment(align), maxCachePerSize(maxCache) {}

    ~SimpleCudaMemoryPool() {
        // 析构时释放所有缓存的内存，防止泄漏
        std::lock_guard<std::mutex> lock(poolMutex);
        for (auto& pair : freeLists) {
            for (void* ptr : pair.second) {
                cudaFree(ptr);
            }
        }
        freeLists.clear();
        std::cout << "[Pool] Destroyed. All memory freed." << std::endl;
    }

    /**
     * 分配内存
     */
    void* allocate(size_t size) {
        if (size == 0) return nullptr;

        // 1. 对齐大小 (关键步骤：让 100 和 101 字节的请求能复用同一个块)
        size_t alignedSize = ((size + alignment - 1) / alignment) * alignment;

        std::lock_guard<std::mutex> lock(poolMutex);

        auto it = freeLists.find(alignedSize);
        
        // 2. 尝试从池中获取
        if (it != freeLists.end() && !it->second.empty()) {
            void* ptr = it->second.back();
            it->second.pop_back();
            // 可选：打印调试信息
            // std::cout << "[Pool] Hit cache for size " << alignedSize << std::endl;
            return ptr;
        }

        // 3. 池中没有，真正分配
        void* ptr = nullptr;
        // 注意：这里传入的是 alignedSize，保证实际分配的内存足够大
        CUDA_CHECK(cudaMalloc(&ptr, alignedSize));
        // std::cout << "[Pool] CUDA Malloc new block: " << alignedSize << " bytes" << std::endl;
        return ptr;
    }


    /**
     * 【推荐接口】如果上层知道大小，用这个函数，性能更高（无需 cudaPointerGetAttributes）
     */
    void deallocate(void* ptr, size_t originalSize) {
        if (ptr == nullptr) return;
        size_t alignedSize = ((originalSize + alignment - 1) / alignment) * alignment;
        std::lock_guard<std::mutex> lock(poolMutex);
        freeLists[alignedSize].push_back(ptr);

        std::vector<void*>& list = freeLists[alignedSize];
        if (list.size() > maxCachePerSize) {
            size_t toFreeCount = list.size() - maxCachePerSize;
            for (size_t i = 0; i < toFreeCount; ++i) {
                cudaFree(list[i]);
            }
            list.erase(list.begin(), list.begin() + toFreeCount);
        }
    }
    
    // 统计信息
    void printStats() {
        std::lock_guard<std::mutex> lock(poolMutex);
        std::cout << "--- Memory Pool Stats ---" << std::endl;
        size_t totalCachedBlocks = 0;
        for (const auto& pair : freeLists) {
            if (!pair.second.empty()) {
                std::cout << "Size: " << pair.first << " | Cached Count: " << pair.second.size() << std::endl;
                totalCachedBlocks += pair.second.size();
            }
        }
        std::cout << "Total Cached Blocks: " << totalCachedBlocks << std::endl;
        std::cout << "-----------------------" << std::endl;
    }
};

static SimpleCudaMemoryPool pool(32, 20);

CudaBuffer::CudaBuffer(
    size_t size_bytes
)
{
    // CUDA_CHECK(cudaMalloc(&ptr_, size_bytes));
    ptr_ = pool.allocate(size_bytes);
    if(ptr_ == nullptr)
    {
        std::cerr << "CudaBuffer" <<std::endl;
    }
    len_ = size_bytes;
    is_ref_ = false;
}

CudaBuffer::~CudaBuffer()
{
    if(ptr_ && (!is_ref_))
    {
        pool.deallocate(ptr_, len_);
        // cudaFree(ptr_);
    }
    ptr_ = 0;
    is_ref_ = true;
}

CudaBuffer::CudaBuffer(CudaBuffer&& other) noexcept
{
    ptr_ = other.ptr_;
    other.ptr_ = 0;
    len_ = other.len_;
    other.len_ = 0;
    is_ref_ = other.is_ref_;
    other.is_ref_ = true;
}
CudaBuffer& CudaBuffer::operator=(CudaBuffer&& other) noexcept
{
    if(ptr_ && (!is_ref_))
    {
        cudaFree(ptr_);
        ptr_ = 0;
        is_ref_ = true;
    }
    ptr_ = other.ptr_;
    other.ptr_ = 0;
    len_ = other.len_;
    other.len_ = 0;
    return *this;
}

void CudaBuffer::copy_to_host(void* host_ptr) const
{
    CUDA_CHECK(cudaMemcpy(host_ptr, ptr_, len_, cudaMemcpyDeviceToHost));

}
void CudaBuffer::copy_from_host(const void* host_ptr)
{
    CUDA_CHECK(cudaMemcpy(ptr_, host_ptr, len_, cudaMemcpyHostToDevice));
}
void CudaBuffer::copy_from_other(const CudaBuffer& other)
{
    assert(len_ == other.len_);
    CUDA_CHECK(cudaMemcpy(ptr_, other.ptr_, len_, cudaMemcpyDeviceToDevice));
}

CudaBuffer CudaBuffer::slice(size_t l, size_t r) const
{
    assert(r>l);
    return CudaBuffer(ptr_ + l, r - l, true);
}

void CudaBuffer::set_zero() const
{
    CUDA_CHECK(cudaMemset(ptr_, 0, len_));
}