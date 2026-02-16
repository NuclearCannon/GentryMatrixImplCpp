#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_modops.cuh"
#include "GPU/cuda_check.hpp"
#include "montgomery.hpp"

constexpr size_t THREAD_PER_GROUP = 64;
constexpr size_t NUMBERS_PER_THREAD = 16;
constexpr size_t NUMBERS_PER_GROUP = THREAD_PER_GROUP * NUMBERS_PER_THREAD;

template<typename Op>
__global__ void _cuda_batch_binary_op_kernel(
    uint64_t* dst,
    const uint64_t* src1,
    const uint64_t* src2,
    size_t batch_size,
    uint64_t M
)
{
    unsigned gid = blockIdx.x;
    unsigned tid = threadIdx.x;
    unsigned id = gid * blockDim.x + tid;
    size_t stt = id * NUMBERS_PER_THREAD;
    size_t end = min(stt + NUMBERS_PER_THREAD, batch_size);

    Op op; // 构造仿函数（无开销）
    for (size_t i = stt; i < end; ++i) {
        dst[i] = op(src1[i], src2[i], M);
    }
}


template<typename Op>
inline float cuda_batch_binary_op_impl(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    uint64_t M
)
{
    dim3 blockSize(THREAD_PER_GROUP);
    dim3 gridSize((batch_size + NUMBERS_PER_GROUP - 1) / NUMBERS_PER_GROUP);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    _cuda_batch_binary_op_kernel<Op><<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src1.get_ptr(),
        (const uint64_t*)src2.get_ptr(),
        batch_size,
        M
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

// 公共接口
float cuda_batch_add(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    uint64_t M
) {
    return cuda_batch_binary_op_impl<ModAdd>(dst, src1, src2, batch_size, M);
}

float cuda_batch_sub(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    uint64_t M
) {
    return cuda_batch_binary_op_impl<ModSub>(dst, src1, src2, batch_size, M);
}




__global__ void _cuda_batch_mul_mont_kernel(
    uint64_t* dst,
    const uint64_t* src1,
    const uint64_t* src2,
    size_t batch_size,
    uint64_t M,
    uint64_t N1
)
{
    unsigned gid = blockIdx.x;
    unsigned tid = threadIdx.x;
    unsigned id = gid * blockDim.x + tid;
    size_t stt = id * NUMBERS_PER_THREAD;
    size_t end = min(stt + NUMBERS_PER_THREAD, batch_size);


    for (size_t i = stt; i < end; ++i) {
        dst[i] = _mul_cuda(src1[i], src2[i], M, N1);
    }
}


float cuda_batch_mul_mont(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    const MontgomeryMultiplier& mm
)
{
    dim3 blockSize(THREAD_PER_GROUP);
    dim3 gridSize((batch_size + NUMBERS_PER_GROUP - 1) / NUMBERS_PER_GROUP);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    _cuda_batch_mul_mont_kernel<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src1.get_ptr(),
        (const uint64_t*)src2.get_ptr(),
        batch_size,
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

__global__ void _cuda_batch_neg_kernel(
    uint64_t* dst,
    const uint64_t* src,
    size_t batch_size,
    uint64_t M
)
{
    unsigned gid = blockIdx.x;
    unsigned tid = threadIdx.x;
    unsigned id = gid * blockDim.x + tid;
    size_t stt = id * NUMBERS_PER_THREAD;
    size_t end = min(stt + NUMBERS_PER_THREAD, batch_size);


    for (size_t i = stt; i < end; ++i) {
        dst[i] = _mod_sub(0, src[i], M);
    }
}


float cuda_batch_neg(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t batch_size,
    uint64_t M
)
{
    dim3 blockSize(THREAD_PER_GROUP);
    dim3 gridSize((batch_size + NUMBERS_PER_GROUP - 1) / NUMBERS_PER_GROUP);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    _cuda_batch_neg_kernel<<<gridSize, blockSize>>>(
        (uint64_t*)dst.get_ptr(),
        (const uint64_t*)src.get_ptr(),
        batch_size,
        M
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