#pragma once
#include "GPU/cuda_buffer.hpp"
#include "montgomery.hpp"

float cuda_ntt(
    CudaBuffer& a,
    const CudaBuffer& roots,
    size_t logn,
    const MontgomeryMultiplier& mm,
    size_t batch_size,
    bool dec
);