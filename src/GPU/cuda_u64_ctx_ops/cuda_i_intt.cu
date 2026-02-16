#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include "montgomery.hpp"

constexpr size_t THREAD_PER_GROUP = 64;



__global__ void _cuda_i_intt_kernel(
    uint64_t* dst,
    const uint64_t* src,
    size_t pnn,
    uint64_t I_inv_mont,
    uint64_t inv2_mont,
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
        uint64_t Pi = src[i], Ni = src[j];
        uint64_t real = _mod_add(Pi, Ni, M);
        uint64_t imag = _mul_cuda(
            _mod_sub(Pi, Ni, M),
            I_inv_mont,
            M, N1
        );
        // 除以2
        dst[i] = _mul_cuda(real, inv2_mont, M, N1);
        dst[j] = _mul_cuda(imag, inv2_mont, M, N1);
    }
}

// 返回耗时（ms）
float cuda_i_intt(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t pnn,
    uint64_t I_inv_mont,
    const MontgomeryMultiplier& mm
)
{
    // 计算2的乘法逆元
    // 不用计算了，就是(M+1)/2
    uint64_t inv2 = mm.encode((mm.getM()+1)/2);
    // 每个线程组分配THREAD_PER_GROUP个线程
    dim3 blockSize(THREAD_PER_GROUP);
    // 分配这么多线程组，向上取整，保证
    // 线程组个数*THREAD_PER_GROUP >= pnn
    dim3 gridSize((pnn + THREAD_PER_GROUP -1)/THREAD_PER_GROUP);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    _cuda_i_intt_kernel<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src.get_ptr(),
        pnn,
        I_inv_mont,
        inv2,
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