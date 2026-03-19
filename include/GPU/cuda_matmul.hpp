#pragma once
#include "GPU/cuda_buffer.hpp"
#include <cstdint>
#include "montgomery.hpp"
#include <cuda_runtime.h>

void matmul_gpu(
    const CudaBuffer& C,
    const CudaBuffer& A,
    const CudaBuffer& B,
    int size,
    const MontgomeryMultiplier& mm
);

void matmul_gpu_stream(
    const CudaBuffer& C,
    const CudaBuffer& A,
    const CudaBuffer& B,
    int size,
    const MontgomeryMultiplier& mm,
    cudaStream_t stream
);