#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_check.hpp"
#include "GPU/cuda_matmul.hpp"
#include "./cuda_modops.cuh"
#include "montgomery.hpp"
#include <cstdint>
#include <cassert>
#include <vector>

constexpr int TILE = 16;
constexpr int WIDTH = 32;   // 这俩可以不等



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
    __shared__ uint64_t As[TILE][WIDTH+1];    // 双缓冲区，异步加载
    __shared__ uint64_t Bs[TILE][WIDTH+1];

    // 领取任务
    const uint64_t* A = Aseq[blockIdx.z];
    const uint64_t* B = Bseq[blockIdx.z];
    uint64_t* C = Cseq[blockIdx.z];


    // 本线程组负责生成C中第[by, bx]小块的第[TILE, TILE]元素
    int row = (blockIdx.y * TILE + threadIdx.y) * n;
    int col = (blockIdx.x * TILE + threadIdx.y) * n;
    int dst = blockIdx.x * TILE + threadIdx.x + row;

    // 直接存储128位结果，最后朴素取模。不做什么蒙哥马利乘法。
    __uint128_t sum = 0;
    uint64_t areg, breg;
    for (int t = 0; t < n / WIDTH; ++t) {
        // 异步加载另一个缓冲区
        // 将A的第[by, t]小块和B的第[bx, t]小块复制到As, Bs

        for(int i=threadIdx.x; i<WIDTH; i+=TILE)
        {
            // 这会执行(n / WIDTH) * (WIDTH/TILE) = n/TILE次
            // n*n 个线程
            // 因此，TILE越大，内存搬运的成本就越小
            As[threadIdx.y][i] = A[row + t * WIDTH + i];
            Bs[threadIdx.y][i] = B[col + t * WIDTH + i];
        }
        
        __syncthreads();
        // 使用当前缓冲区
        // 每个线程各自的sum组成的矩阵 = As @ Bs
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
    C[dst] = sum;
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
    // 提交运行
    dim3 blockSize(TILE, TILE);
    dim3 gridSize(n_/TILE, n_/TILE, ntask);
    matmul_batch_kernel<<<gridSize, blockSize>>>(
        (const uint64_t**)Aseq.get_ptr(), (const uint64_t**)Bseq.get_ptr(), (uint64_t**)Cseq.get_ptr(), M_, n_
    );
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


