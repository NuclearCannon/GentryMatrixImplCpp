#pragma once
#include "GPU/cuda_buffer.hpp"
#include <cstdint>
#include "montgomery.hpp"

void matmul_gpu(
    CudaBuffer& C,
    const CudaBuffer& A,
    const CudaBuffer& B,
    int size,
    const MontgomeryMultiplier& mm
);