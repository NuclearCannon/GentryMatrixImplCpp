#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include "montgomery.hpp"

// 用这组参数就没问题
constexpr int THREAD_PER_GROUP = 32;
constexpr int ROWS_PER_THREAD = 32;


// 暂行方案：
// 每个线程负责一行
// 共开batch_size个线程
__global__ void _cuda_batch_mul_vec_kernel(
    uint64_t* dst,
    const uint64_t* src,
    const uint64_t* vec,
    size_t batch_size,
    size_t vec_len,
    uint64_t M,
    uint64_t N1
)
{
    unsigned gid = blockIdx.x;
    unsigned tid = threadIdx.x;
    unsigned id = gid * blockDim.x + tid;
    if(id < batch_size)
    {
        dst += id * vec_len;
        src += id * vec_len;
        for(size_t i=0; i<vec_len; i++)
        {
            dst[i] = _mul_cuda(src[i], vec[i], M, N1);
        }
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
    
    dim3 blockSize(THREAD_PER_GROUP);
    dim3 gridSize((batch_size + ROWS_PER_THREAD - 1) / ROWS_PER_THREAD);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    _cuda_batch_mul_vec_kernel<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src.get_ptr(),
        (const uint64_t*)vec.get_ptr(),
        batch_size, vec_len,
        mm.getM(), mm.getN1()
    );

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return ms;
}