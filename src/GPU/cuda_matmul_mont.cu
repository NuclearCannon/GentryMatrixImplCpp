#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_check.hpp"
#include "GPU/cuda_matmul.hpp"
#include "montgomery.hpp"
#include <cstdint>
#include <cassert>


__device__ uint64_t _montgomery_reduce_cu(__uint128_t t, uint64_t M, uint64_t N1)  {
    uint64_t m = (uint64_t)t * N1;
    uint64_t v = ((t + (__uint128_t)m * M) >> 64);
    return (v<M)?v:v-M;
}

__device__ uint64_t _mul_cu(uint64_t a, uint64_t b, uint64_t M, uint64_t N1) {
    return _montgomery_reduce_cu((__uint128_t)a * b, M, N1);
}


__global__ void matmul_tiled(
    const uint64_t* __restrict__ A,
    const uint64_t* __restrict__ B,
    uint64_t* __restrict__ C,
    uint64_t M, // 模数
    uint64_t N1, // 模数
    int size    // 矩阵边长
) {
    constexpr int TILE_SIZE = 16;
    __shared__ uint64_t As[TILE_SIZE][TILE_SIZE];
    __shared__ uint64_t Bs[TILE_SIZE][TILE_SIZE];

    int row = blockIdx.y * TILE_SIZE + threadIdx.y;
    int col = blockIdx.x * TILE_SIZE + threadIdx.x;

    uint64_t sum = 0;
    for (int t = 0; t < (size + TILE_SIZE - 1) / TILE_SIZE; ++t) {
        if (row < size && t * TILE_SIZE + threadIdx.x < size)
            As[threadIdx.y][threadIdx.x] = A[row * size + t * TILE_SIZE + threadIdx.x];
        else
            As[threadIdx.y][threadIdx.x] = 0;

        if (col < size && t * TILE_SIZE + threadIdx.y < size)
            Bs[threadIdx.y][threadIdx.x] = B[col * size + t * TILE_SIZE + threadIdx.y];
        else
            Bs[threadIdx.y][threadIdx.x] = 0;

        __syncthreads();

        for (int k = 0; k < TILE_SIZE; ++k)
        {
            sum += _mul_cu(As[threadIdx.y][k], Bs[k][threadIdx.x], M, N1);
            if(sum>M)sum-=M;
        }
        __syncthreads();
    }

    if (row < size && col < size)
        C[row * size + col] = sum;
}




void matmul_gpu(
    CudaBuffer& C,
    const CudaBuffer& A,
    const CudaBuffer& B,
    int size,
    const MontgomeryMultiplier& mm
)
{
    // std::cout << "matmul_gpu" << std::endl;
    const uint64_t* Ap = (const uint64_t*)(A.get_ptr());
    const uint64_t* Bp = (const uint64_t*)(B.get_ptr());
    uint64_t* Cp = (uint64_t*)(C.get_ptr());

    size_t bytes = size * size * sizeof(uint64_t);
    assert(A.size() == bytes);
    assert(B.size() == bytes);
    assert(C.size() == bytes);
    dim3 blockSize(16, 16);
    dim3 gridSize((size + blockSize.x - 1) / blockSize.x,
                  (size + blockSize.y - 1) / blockSize.y);
    
    matmul_tiled<<<gridSize, blockSize>>>(Ap, Bp, Cp, mm.getM(), mm.getN1(), size);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}