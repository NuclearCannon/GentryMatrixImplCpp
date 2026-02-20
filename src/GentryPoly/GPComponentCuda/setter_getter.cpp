#include "GentryPoly.hpp"


void GPComponentCuda::set_from_other(const GPComponentCuda& other)
{
    assert(like(other));
    data_.copy_from_other(other.data_);
}
void GPComponentCuda::set_from_cpu(const GPComponent& other)
{
    assert(like(other));
    data_.copy_from_host(other.data_.data());
}
void GPComponentCuda::to_buffer(uint64_t* data)
{
    data_.copy_to_host(data);
}
void GPComponentCuda::to_cpu(GPComponent& other)
{
    assert(like(other));
    to_buffer(other.data_.data());
}