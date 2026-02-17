#include "u64_array.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"
#include <assert.h>

void CRTArrayGPU::iw_ntt_inplace()
{
    int len = cuda_data_.size();
    for(int i=0; i<len; i++)
    {
        cc_->get_ctx()[i]->iw_ntt_cuda(
            *cuda_data_[i],
            *cuda_data_[i]
        );
    }
}

void CRTArrayGPU::iw_intt_inplace()
{
    int len = cuda_data_.size();
    for(int i=0; i<len; i++)
    {
        cc_->get_ctx()[i]->iw_intt_cuda(
            *cuda_data_[i],
            *cuda_data_[i]
        );
    }
}

void CRTArrayGPU::xy_ntt_inplace()
{
    int len = cuda_data_.size();
    for(int i=0; i<len; i++)
    {
        cc_->get_ctx()[i]->xy_ntt_cuda(
            *cuda_data_[i],
            *cuda_data_[i]
        );
    }
}

void CRTArrayGPU::xy_intt_inplace()
{
    int len = cuda_data_.size();
    for(int i=0; i<len; i++)
    {
        cc_->get_ctx()[i]->xy_intt_cuda(
            *cuda_data_[i],
            *cuda_data_[i]
        );
    }
}