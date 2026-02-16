#include <cstdint>
#include <assert.h>
#include "GPU/cuda_check.hpp"
#include "GPU/cuda_buffer.hpp"

constexpr int TILE_SIZE = 8;


__global__ void _transpose_rect_restrict_kernel(
    uint64_t* __restrict__ dst,
    const uint64_t* __restrict__ src,
    size_t r, size_t c
)
{
    unsigned ii = blockIdx.x * TILE_SIZE;
    unsigned jj = blockIdx.y * TILE_SIZE;
    unsigned i = ii + threadIdx.x;
    unsigned j = jj + threadIdx.y;
    dst[j * r + i] = src[i * c + j];
}


float cuda_transpose_rect_restrict(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t r, size_t c
)
{
    assert(dst.size() == r*c*sizeof(uint64_t));
    assert(src.size() == r*c*sizeof(uint64_t));

    uint64_t * dstp = (uint64_t*)dst.get_ptr();
    const uint64_t * srcp = (const uint64_t*)src.get_ptr();


    assert(srcp!=dstp);
    assert(r % TILE_SIZE == 0);
    assert(c % TILE_SIZE == 0);


    dim3 blockSize(TILE_SIZE, TILE_SIZE);
    dim3 gridSize(r/TILE_SIZE, c/TILE_SIZE);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    _transpose_rect_restrict_kernel<<<gridSize, blockSize>>>(dstp, srcp, r, c);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    return milliseconds;
}