#include "GentryPoly.hpp"
#include "modops.hpp"


void GPComponent::neg(GPComponent& dst, const GPComponent& src)
{
    assert(dst.like(src));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    for(size_t i=0; i<size; i++)dst.data_[i] = mod_sub(0, src.data_[i], q);
}
void GPComponent::add(GPComponent& dst, const GPComponent& src1, const GPComponent& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    for(size_t i=0; i<size; i++)dst.data_[i] = mod_add(src1.data_[i], src2.data_[i], q);
}
void GPComponent::sub(GPComponent& dst, const GPComponent& src1, const GPComponent& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    for(size_t i=0; i<size; i++)dst.data_[i] = mod_sub(src1.data_[i], src2.data_[i], q);
}
void GPComponent::mul(GPComponent& dst, const GPComponent& src1, const GPComponent& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    for(size_t i=0; i<size; i++)dst.data_[i] = mod_mul(src1.data_[i], src2.data_[i], q);
}
void GPComponent::mul_scalar(GPComponent& dst, const GPComponent& src1, uint64_t src_scalar)
{
    assert(dst.like(src1));
    size_t size = dst.get_size();
    uint64_t encoded = dst.mm_.encode(src_scalar);
    for(size_t i=0; i<size; i++)dst.data_[i] = dst.mm_.mul(src1.data_[i], encoded);
}
void GPComponent::mont_encode(GPComponent& dst, const GPComponent& src)
{
    assert(dst.like(src));
    size_t size = dst.get_size();
    for(size_t i=0; i<size; i++)dst.data_[i] = dst.mm_.encode(src.data_[i]);
}
void GPComponent::mont_decode(GPComponent& dst, const GPComponent& src)
{
    assert(dst.like(src));
    size_t size = dst.get_size();
    for(size_t i=0; i<size; i++)dst.data_[i] = dst.mm_.decode(src.data_[i]);
}
void GPComponent::mul_mont(GPComponent& dst, const GPComponent& src1, const GPComponent& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    for(size_t i=0; i<size; i++)dst.data_[i] = dst.mm_.mul(src1.data_[i], src2.data_[i]);
}