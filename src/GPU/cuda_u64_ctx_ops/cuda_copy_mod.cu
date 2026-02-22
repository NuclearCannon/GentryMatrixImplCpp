#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include <cassert>
#include <cuda_runtime.h>

// CUDA kernel: 对每个 uint64_t 元素取模 M
__global__ void _copy_mod_kernel(
    uint64_t* dst,
    const uint64_t* src,
    uint64_t M,
    size_t batch_size
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < batch_size) {
        dst[idx] = src[idx] % M;
    }
}

float cuda_copy_mod(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    uint64_t M,
    size_t batch_size
) {
    assert(dst.size() == batch_size * sizeof(uint64_t));
    assert(src.size() == batch_size * sizeof(uint64_t));
    assert(M > 0); // 避免除零或无意义模

    const uint64_t* srcp = static_cast<const uint64_t*>(src.get_ptr());
    uint64_t* dstp = static_cast<uint64_t*>(dst.get_ptr());

    // 配置 CUDA 执行参数
    const int block_size = 256;
    const int grid_size = (batch_size + block_size - 1) / block_size;

    // 可选：记录时间（如果返回 float 是为了表示耗时）
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    _copy_mod_kernel<<<grid_size, block_size>>>(dstp, srcp, M, batch_size);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float milliseconds = 0.0f;
    cudaEventElapsedTime(&milliseconds, start, stop);

    // 清理事件
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    // 检查 kernel 是否出错
    CUDA_CHECK(cudaGetLastError());
    return milliseconds; // 返回 kernel 执行时间（毫秒）
}