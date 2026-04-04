#include "GentryPoly.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"

void GPComponentCuda::set_from_other(const GPComponentCuda& other)
{
    assert(like(other));
    data_.copy_from_other(other.data_);
}
void GPComponentCuda::set_from_cuda_buffer(const CudaBuffer& buf)
{
    assert(buf.size() == get_size()*sizeof(uint64_t));
    cuda_copy_mod(data_, buf, q_, get_size());
}
void GPComponentCuda::set_from_cpu(const GPComponent& other)
{
    assert(like(other));
    data_.copy_from_host(other.data_.data());
}
void GPComponentCuda::to_buffer(uint64_t* data) const
{
    data_.copy_to_host(data);
}
void GPComponentCuda::to_cpu(GPComponent& other) const
{
    assert(like(other));
    to_buffer(other.data_.data());
}