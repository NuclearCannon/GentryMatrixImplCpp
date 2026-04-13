#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_check.hpp"
#include "GPU/cuda_matmul.hpp"
#include "./cuda_modops.cuh"
#include <cstdint>
#include <cassert>
#include <vector>

constexpr int TILE = 16;
constexpr int WIDTH = 32;   // 这俩可以不等


__global__ void matmul_batch_kernel_flex(
    const uint64_t** Aseq,
    const uint64_t** Bseq,
    uint64_t** Cseq,
    uint64_t M,
    int n
)
{
    // 领取任务
    const uint64_t* A = Aseq[blockIdx.z];
    const uint64_t* B = Bseq[blockIdx.z];
    uint64_t* C = Cseq[blockIdx.z];
    int i = blockIdx.x * gridDim.x + threadIdx.x;
    int j = blockIdx.y * gridDim.y + threadIdx.y;
    int in = i*n;
    int jn = j*n;
    __uint128_t sum = 0;
    for(int k=0; k<n; k++)
    {
        sum += ((__uint128_t)(A[in+k])) * B[jn+k];
    }
    sum %= M;
    if(i<n && j<n)
    {
        C[in+j] = sum;
    }


}


/// 对于一个算例，有[n/TILE, n/TILE]个block，每个block尺寸为[TILE, TILE]
/// 总共有n^2个线程处理一个算例
__global__ void matmul_batch_kernel(
    const uint64_t** Aseq,
    const uint64_t** Bseq,
    uint64_t**  Cseq,
    uint64_t M, // 模数
    int n    // 矩阵边长
)
{
    // 我们会把输入的A,B都切分成多个小块，每个小块的尺寸是[TILE, WIDTH]，做分块矩阵乘法
    __shared__ uint64_t As[TILE][WIDTH+1];  // 用于存放一个小矩阵
    __shared__ uint64_t Bs[TILE][WIDTH+1];

    // 领取任务
    const uint64_t* A = Aseq[blockIdx.z];
    const uint64_t* B = Bseq[blockIdx.z];
    uint64_t* C = Cseq[blockIdx.z];

    // 本线程负责计算C[glb_y, glb_x]
    // glb_x = blockIdx.x * TILE + threadIdx.x;
    // glb_y = blockIdx.y * TILE + threadIdx.y;

    // 本线程组负责生成C中第[blockIdx.y, blockIdx.x]小块
    // 需要访问
    // A中的[blockIdx.y, :]块列
    // B中的[blockIdx.x, :]块列
    // 在搬运中，本块需要搬运一个块内的[threadIdx.y, threadIdx.x:WIDTH:TILE]部分
    // 预计算一些偏移量
    int offsetA = (blockIdx.y * TILE + threadIdx.y) * n;    // blockIdx.y * TILE表示自己负责哪个块，threadIdx.y表示自己负责块内哪个行
    int offsetB = (blockIdx.x * TILE + threadIdx.y) * n;

    // 直接存储128位结果，最后朴素取模。不做什么蒙哥马利乘法。
    __uint128_t sum = 0;
    uint64_t areg, breg;
    for (int t = 0; t < n / WIDTH; ++t) {
        // 将A的第[blockIdx.y, t]小块和B的第[blockIdx.x, t]小块复制到As, Bs
        // 注意，
        for(int i=threadIdx.x; i<WIDTH; i+=TILE)
        {
            // 这会执行(n / WIDTH) * (WIDTH/TILE) = n/TILE次
            // n*n 个线程
            // 因此，TILE越大，内存搬运的成本就越小
            As[threadIdx.y][i] = A[offsetA + t * WIDTH + i];
            Bs[threadIdx.y][i] = B[offsetB + t * WIDTH + i];
        }
        
        __syncthreads();
        // 每个线程各自的sum组成的矩阵 += As @ Bs.T
        for (int k = 0; k < WIDTH; ++k)
        {
            areg = As[threadIdx.y][k];
            breg = Bs[threadIdx.x][k];
            sum += (__uint128_t)areg * breg; // 延迟取模
        }
        // 交换缓冲区
        __syncthreads();
    }

    sum %= M;
    int glb_x = blockIdx.x * TILE + threadIdx.x;
    int glb_y = blockIdx.y * TILE + threadIdx.y;
    C[glb_y*n + glb_x] = sum;
}



void CudaMatmulTaskSet::run() {

    int ntask = Atasks.size();
    if(ntask == 0)return;

    CudaBuffer Aseq(ntask * sizeof(uint64_t*));
    CudaBuffer Bseq(ntask * sizeof(uint64_t*));
    CudaBuffer Cseq(ntask * sizeof(uint64_t*));
    Aseq.copy_from_host(Atasks.data());
    Bseq.copy_from_host(Btasks.data());
    Cseq.copy_from_host(Ctasks.data());

    bool aligned = (n_ % TILE == 0) && (n_ % WIDTH == 0);

    if (aligned) {
        dim3 blockSize(TILE, TILE);
        dim3 gridSize(n_/TILE, n_/TILE, ntask);
        matmul_batch_kernel<<<gridSize, blockSize>>>(
            (const uint64_t**)Aseq.get_ptr(), (const uint64_t**)Bseq.get_ptr(), (uint64_t**)Cseq.get_ptr(), M_, n_
        );
    } else {
        dim3 blockSize(TILE, TILE);
        dim3 gridSize((n_ + TILE - 1) / TILE, (n_ + TILE - 1) / TILE, ntask);
        matmul_batch_kernel_flex<<<gridSize, blockSize>>>(
            (const uint64_t**)Aseq.get_ptr(), (const uint64_t**)Bseq.get_ptr(), (uint64_t**)Cseq.get_ptr(), M_, n_
        );
    }
}


void CudaMatmulTaskSet::append(const CudaBuffer& C,const CudaBuffer& A,const CudaBuffer& B) {
    // 检查输入合法性
    uint64_t* Ap = (uint64_t*)(A.get_ptr());
    uint64_t* Bp = (uint64_t*)(B.get_ptr());
    uint64_t* Cp = (uint64_t*)(C.get_ptr());
    size_t bytes = n_ * n_ * sizeof(uint64_t);
    assert(A.size() == bytes);
    assert(B.size() == bytes);
    assert(C.size() == bytes);
    assert(Cp != Ap);
    assert(Cp != Bp);
    // 加入任务组
    Atasks.push_back(Ap);
    Btasks.push_back(Bp);
    Ctasks.push_back(Cp);

}


