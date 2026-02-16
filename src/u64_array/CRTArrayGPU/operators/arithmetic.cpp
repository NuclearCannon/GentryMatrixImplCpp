#include "u64_array.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"
#include <assert.h>

void CRTArrayGPU::add(const CRTArrayGPU& src1, const CRTArrayGPU& src2)
{
    size_t len = cuda_data_.size();
    assert(src1.cuda_data_.size() == len);
    assert(src2.cuda_data_.size() == len);
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        cuda_batch_add(
            *cuda_data_[i],
            *src1.cuda_data_[i],
            *src2.cuda_data_[i],
            batch_size,
            cc_->get_mods()[i]
        );
    }
}

void CRTArrayGPU::sub(const CRTArrayGPU& src1, const CRTArrayGPU& src2)
{
    size_t len = cuda_data_.size();
    assert(src1.cuda_data_.size() == len);
    assert(src2.cuda_data_.size() == len);
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        cuda_batch_sub(
            *cuda_data_[i],
            *src1.cuda_data_[i],
            *src2.cuda_data_[i],
            batch_size,
            cc_->get_mods()[i]
        );
    }
}

void CRTArrayGPU::neg_inplace()
{
    size_t len = cuda_data_.size();
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        cuda_batch_neg(
            *cuda_data_[i],
            *cuda_data_[i],
            batch_size,
            cc_->get_mods()[i]
        );
    }
}


void CRTArrayGPU::mul_mont(const CRTArrayGPU& src1, const CRTArrayGPU& src2)
{
    size_t len = cuda_data_.size();
    assert(src1.cuda_data_.size() == len);
    assert(src2.cuda_data_.size() == len);
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        cuda_batch_mul_mont(
            *cuda_data_[i],
            *src1.cuda_data_[i],
            *src2.cuda_data_[i],
            batch_size,
            cc_->get_ctx()[i]->get_multiplier()
        );
    }
}

void CRTArrayGPU::mul_scalar_inplace(uint64_t scalar)
{
    size_t len = cuda_data_.size();
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        const auto& mm = cc_->get_ctx()[i]->get_multiplier();
        uint64_t scalar_encoded = mm.encode(scalar);

        cuda_batch_mul_scalar(
            *cuda_data_[i],
            *cuda_data_[i],
            scalar_encoded,
            batch_size,
            mm
        );
    }
}

void CRTArrayGPU::mont_encode_inplace()
{
    size_t len = cuda_data_.size();
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        const auto& mm = cc_->get_ctx()[i]->get_multiplier();
        // encode 相当于乘以一个R2
        cuda_batch_mul_scalar(
            *cuda_data_[i],
            *cuda_data_[i],
            mm.getR2(),
            batch_size,
            mm
        );
    }
}

void CRTArrayGPU::mont_decode_inplace()
{
    size_t len = cuda_data_.size();
    size_t batch_size = cc_->get_size();
    for(size_t i=0; i<len; i++)
    {
        const auto& mm = cc_->get_ctx()[i]->get_multiplier();
        // decode 相当于乘以一个1
        cuda_batch_mul_scalar(
            *cuda_data_[i],
            *cuda_data_[i],
            1,
            batch_size,
            mm
        );
    }
}