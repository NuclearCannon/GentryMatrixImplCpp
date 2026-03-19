#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_check.hpp"
#include "GPU/cuda_matmul.hpp"
#include "GPU/cuda_modops.cuh"
#include "montgomery.hpp"
#include <cstdint>
#include <cassert>



// 调用此函数时，请令blockSize=(16, 16), gridSize=(size/16, size/16)
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
            sum += _mul_cuda(As[threadIdx.y][k], Bs[k][threadIdx.x], M, N1);
            if(sum>M)sum-=M;
        }
        __syncthreads();
    }

    if (row < size && col < size)
        C[row * size + col] = sum;
}




void matmul_gpu(
    const CudaBuffer& C,
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
    
    matmul_tiled<<<gridSize, blockSize>>>(Ap, Bp, Cp, mm.M, mm.N1, size);
    CUDA_CHECK(cudaGetLastError());
    // CUDA_CHECK(cudaDeviceSynchronize());
}

// matmul_gpu的制定流版本
void matmul_gpu_stream(
    const CudaBuffer& C,
    const CudaBuffer& A,
    const CudaBuffer& B,
    int size,
    const MontgomeryMultiplier& mm,
    cudaStream_t stream
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

    
    
    matmul_tiled<<<gridSize, blockSize, 0, stream>>>(Ap, Bp, Cp, mm.M, mm.N1, size);
    CUDA_CHECK(cudaGetLastError());
    // CUDA_CHECK(cudaDeviceSynchronize());
}

// 这个函数会打印matmul_tiled要想占满GPU
// 没人调用它是正常的
void matmul_tiled_report()
{
    int minGridSize = -1;
    int blockSize = -1;
    auto error = cudaOccupancyMaxPotentialBlockSize(
        &minGridSize, &blockSize, matmul_tiled, 0, 256
    );
    if(error != cudaSuccess)
    {
        printf("matmul_tiled_report: Error发生\n");
    }
    printf("matmul_tiled_report:\n");
    printf("blockSize  \t %d\n", blockSize);
    printf("minGridSize\t %d\n", minGridSize);
    
}