#include "GPU/cuda_buffer.hpp"
#include "../cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include "montgomery.hpp"

// 核函数：执行 dst[i] = src[i] * vec[i % vec_len]
__global__ void _cuda_batch_mul_vec_kernel_ver2(
    uint64_t* __restrict__ dst,
    const uint64_t* __restrict__ src,
    const uint64_t* __restrict__ vec,
    int vec_len,
    int total_size,
    uint64_t M,
    uint64_t N1
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < total_size) {
        int vec_idx = idx % vec_len;
        dst[idx] = _mul_cuda(src[idx], vec[vec_idx], M, N1);
    }
}



float cuda_batch_mul_vec(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    const CudaBuffer& vec,
    size_t batch_size,
    size_t vec_len,
    const MontgomeryMultiplier& mm
)
{
    assert(dst.size() == batch_size*vec_len*sizeof(uint64_t));
    assert(src.size() == batch_size*vec_len*sizeof(uint64_t));
    assert(vec.size() == vec_len*sizeof(uint64_t));
    

    // cudaEvent_t start, stop;
    // cudaEventCreate(&start);
    // cudaEventCreate(&stop);
    // cudaEventRecord(start);
    size_t total_size = batch_size * vec_len;
    int blockSize = 256;
    int gridSize = (total_size + blockSize - 1) / blockSize;

    _cuda_batch_mul_vec_kernel_ver2<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src.get_ptr(),
        (const uint64_t*)vec.get_ptr(),
        vec_len,
        total_size, 
        mm.M, mm.N1
    );

    // cudaEventRecord(stop);
    // cudaEventSynchronize(stop);
    float ms = 0;
    // cudaEventElapsedTime(&ms, start, stop);
    CUDA_CHECK(cudaGetLastError());
    // CUDA_CHECK(cudaDeviceSynchronize());

    // cudaEventDestroy(start);
    // cudaEventDestroy(stop);
    return ms;
}