#include "u64_array.hpp"

// 从内存的vv64中载入数据
void CRTArrayGPU::set_from_vv64(const vv64& data)
{
    size_t len = data.size();
    assert(cuda_data_.size() == len);
    for(int i=0; i<len; i++)
    {
        assert(cuda_data_[i]->size() == data[i].size()*sizeof(uint64_t));
        cuda_data_[i]->copy_from_host(data[i].data());
    }
}
// 将数据写入内存（vv64格式）
vv64 CRTArrayGPU::export_to_vv64() const
{
    vv64 result;
    for(const auto& ptr: cuda_data_)
    {
        vec64 r(ptr->size() / sizeof(uint64_t));
        ptr->copy_to_host(r.data());
        result.push_back(std::move(r));
    }
    return result;
}
// 设为全0
void CRTArrayGPU::set_to_zero()
{
    for(const auto& ptr: cuda_data_) ptr->set_zero();
}
// 从另一个对象处拷贝（取代拷贝赋值）
void CRTArrayGPU::copy_from(const CRTArrayGPU& other)
{
    int len = cuda_data_.size();

    assert(len == other.cuda_data_.size());
    for(int i=0; i<len; i++)
    {
        cuda_data_[i]->copy_from_other(*(cuda_data_[i]));
    }
}