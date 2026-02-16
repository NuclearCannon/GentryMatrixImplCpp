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

float cuda_batch_add(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    uint64_t M
);

float cuda_batch_sub(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    uint64_t M
);

float cuda_batch_mul_mont(
    const CudaBuffer& dst,
    const CudaBuffer& src1,
    const CudaBuffer& src2,
    size_t batch_size,
    const MontgomeryMultiplier& mm
);

float cuda_batch_neg(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t batch_size,
    uint64_t M
);

float cuda_batch_mul_scalar(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    uint64_t scalar_encoded,
    size_t batch_size,
    const MontgomeryMultiplier& mm
);

float cuda_batch_mul_vec(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    const CudaBuffer& vec,
    size_t batch_size,
    size_t vec_len,
    const MontgomeryMultiplier& mm
);

float cuda_transpose_rect_restrict(
    const CudaBuffer& dst,
    const CudaBuffer& src,
    size_t r, size_t c
);