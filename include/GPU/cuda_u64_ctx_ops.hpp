#pragma once
#include "GPU/cuda_buffer.hpp"
#include "montgomery.hpp"

float cuda_i_ntt(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t pnn,
    uint64_t I_mont,
    const MontgomeryMultiplier& mm
);

float cuda_i_intt(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t pnn,
    uint64_t I_inv_mont,
    const MontgomeryMultiplier& mm
);