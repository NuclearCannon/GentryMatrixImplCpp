#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include "montgomery.hpp"

constexpr size_t THREAD_PER_GROUP = 64;



__global__ void _cuda_i_ntt_kernel(
    uint64_t* dst,
    const uint64_t* src,
    size_t pnn,
    uint64_t I_mont,
    uint64_t M,
    uint64_t N1
)
{
    unsigned gid = blockIdx.x;// 线程组id
    unsigned tid = threadIdx.x; // 组内线程id
    size_t i = tid + THREAD_PER_GROUP*gid;
    if (i<pnn)
    {
        size_t j = i + pnn;
        uint64_t real_i = src[i];
        uint64_t image_I_i = _mul_cuda(src[j], I_mont, M, N1);
        dst[i] = _mod_add(real_i, image_I_i, M);
        dst[j] = _mod_sub(real_i, image_I_i, M);
    }
}

// 返回耗时（ms）
float cuda_i_ntt(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t pnn,
    uint64_t I_mont,
    const MontgomeryMultiplier& mm
)
{

    // 每个线程组分配THREAD_PER_GROUP个线程
    dim3 blockSize(THREAD_PER_GROUP);
    // 分配这么多线程组，向上取整，保证
    // 线程组个数*THREAD_PER_GROUP >= pnn
    dim3 gridSize((pnn + THREAD_PER_GROUP -1)/THREAD_PER_GROUP);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    _cuda_i_ntt_kernel<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src.get_ptr(),
        pnn,
        I_mont,
        mm.getM(),
        mm.getN1()
    );

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    return ms;
}