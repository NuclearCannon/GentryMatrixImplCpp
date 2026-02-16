#include "u64_array.hpp"

// 未初始化构造
CRTArrayGPU::CRTArrayGPU(std::shared_ptr<const U64CtxChain> cc):
    cc_(cc),
    cuda_data_(cc->get_chain_length())
{
    // 分配空间
    size_t size = cc->get_size();
    for(int i=0; i<cc->get_chain_length(); i++)
    {
        cuda_data_[i] = std::make_unique<CudaBuffer>(size);
    }
}

// 移动
CRTArrayGPU::CRTArrayGPU(CRTArrayGPU&& other) = default;

// 析构
CRTArrayGPU::~CRTArrayGPU() = default;