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