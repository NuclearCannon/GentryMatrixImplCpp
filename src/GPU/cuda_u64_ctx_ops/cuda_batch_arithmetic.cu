#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include "montgomery.hpp"

constexpr size_t THREAD_PER_GROUP = 64;
constexpr size_t NUMBERS_PER_THREAD = 16;
constexpr size_t NUMBERS_PER_GROUP = THREAD_PER_GROUP * NUMBERS_PER_THREAD;

__global__ void _cuda_batch_add_kernel(
    uint64_t* dst,
    const uint64_t* src1,
    const uint64_t* src2,
    size_t batch_size,
    uint64_t M
)
{
    unsigned gid = blockIdx.x;// 线程组id
    unsigned tid = threadIdx.x; // 组内线程id
    unsigned id = gid * blockDim.x + tid;   // 全局id
    size_t stt = id * NUMBERS_PER_THREAD;
    size_t end = stt + NUMBERS_PER_THREAD;
    if(end > batch_size) end = batch_size;
    for(size_t i=stt; i<end; i++)
    {
        // _mod_add是GPU/cuda_modops.cuh中定义的__device__模加函数
        dst[i] = _mod_add(src1[i], src2[i], M);
    }

}

// 返回耗时（ms）
float cuda_batch_add(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    uint64_t M
)
{
    // 每个线程组分配THREAD_PER_GROUP个线程
    dim3 blockSize(THREAD_PER_GROUP);
    dim3 gridSize((batch_size + NUMBERS_PER_GROUP -1)/NUMBERS_PER_GROUP);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    _cuda_batch_add_kernel<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src1.get_ptr(),
        (const uint64_t*)src2.get_ptr(),
        batch_size,
        M
    );
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    return ms;
}