#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <cassert>
#include <cstdint>
#include "GPU/cuda_check.hpp"
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "montgomery.hpp"


static void __global__ _butterfly_inc_mont_cuda(
    uint64_t* a, 
    const uint64_t* omegas_mont,
    size_t logn,
    const uint64_t M,   // 模数
    const uint64_t N1  // 蒙哥马利约简辅助量
)
{
    // ====准备工作====
    // 获取当前线程组id
    size_t gid = blockIdx.x;
    a += gid << logn;   // 这是自己负责的片区
    // 获取组内线程id
    // 一个组内应当有n/2个线程，每个线程负责一个蝴蝶操作
    size_t tid = threadIdx.x;
    size_t wtid = tid << 1; // double tid
    // ====初始化各变量====
    size_t m = 1;
    size_t t = logn;    // t=log(n/m)
    const size_t n = 1<<logn;
    while(m<n)
    {
        size_t m_half = m;
        m <<= 1;
        t--;
        // 不是循环k,j，而是从tid中计算kj
        // j部分，不超过m_half
        size_t j = tid & (m_half-1);
        // k部分，是m的倍数
        size_t k = wtid & (~(m-1));
        size_t i1 = k | j;
        size_t i2 = i1 | m_half;
        uint64_t w = omegas_mont[(j<<t) & (n-1)];
        uint64_t u = a[i1];
        uint64_t v = _mul_cuda(a[i2], w, M, N1);
        // 同步以等待对a的读取完毕
        __syncthreads();
        a[i1] = _mod_add(u,v,M);
        a[i2] = _mod_sub(u,v,M);
        // 同步以等待对a的修改完毕
        __syncthreads();
    }
}


static void __global__ _butterfly_dec_mont_cuda(
    uint64_t* a, 
    const uint64_t* omegas_mont,
    size_t logn,
    const uint64_t M,   // 模数
    const uint64_t N1  // 蒙哥马利约简辅助量
)
{
    // ====准备工作====
    // 获取当前线程组id
    size_t gid = blockIdx.x;
    a += gid << logn;   // 这是自己负责的片区
    // 获取组内线程id
    // 一个组内应当有n/2个线程，每个线程负责一个蝴蝶操作
    size_t tid = threadIdx.x;
    size_t wtid = tid << 1; // double tid
    // ====初始化各变量====
    const size_t n = 1<<logn;
    size_t m = n;
    size_t t = 0;
    
    while(m>1)
    {
        size_t m_half = m>>1;
        // 不是循环k,j，而是从tid中计算kj
        // j部分，不超过m_half
        size_t j = tid & (m_half-1);
        // k部分，是m的倍数
        size_t k = wtid & (~(m-1));
        size_t i1 = k | j;
        size_t i2 = i1 | m_half;
        uint64_t w = omegas_mont[(j<<t) & (n-1)];
        uint64_t u = a[i1];
        uint64_t v = a[i2];
        // 同步以等待对a的读取完毕
        __syncthreads();
        a[i1] = _mod_add(u,v,M);
        a[i2] = _mul_cuda(_mod_sub(u,v,M), w, M, N1);
        // 同步以等待对a的修改完毕
        __syncthreads();
        m = m_half;
        t++;
    }
}



float cuda_ntt(
    CudaBuffer& a,
    const CudaBuffer& roots,
    size_t logn,
    const MontgomeryMultiplier& mm,
    size_t batch_size,
    bool dec
)
{
    uint32_t n = 1<<logn;
    uint64_t* ap = (uint64_t*)(a.get_ptr());
    const uint64_t* rp = (const uint64_t*)(roots.get_ptr());

    size_t bytes = n * batch_size * sizeof(uint64_t);
    assert(a.size() == bytes);

    dim3 blockSize(n/2);
    dim3 gridSize(batch_size);
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    if(dec)
    {
        _butterfly_dec_mont_cuda<<<gridSize, blockSize>>>(
            ap, rp, logn, mm.getM(), mm.getN1()
        );
    }
    else
    {
        _butterfly_inc_mont_cuda<<<gridSize, blockSize>>>(
            ap, rp, logn, mm.getM(), mm.getN1()
        );
    }
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    return milliseconds * 1000;
}