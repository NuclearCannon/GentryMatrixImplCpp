#include "GentryPoly.hpp"
#include "modops.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"

void GPComponentCuda::neg(GPComponentCuda& dst, const GPComponentCuda& src)
{
    assert(dst.like(src));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    cuda_batch_neg(dst.data_, src.data_, size, q);

}
void GPComponentCuda::add(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    cuda_batch_add(dst.data_, src1.data_, src2.data_, size, q);

}
void GPComponentCuda::sub(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    cuda_batch_sub(dst.data_, src1.data_, src2.data_, size, q);
}
void GPComponentCuda::mul(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    // 用一次蒙哥马利乘和一次编码模拟逐位乘的效果
    cuda_batch_mul_mont(dst.data_, src1.data_, src2.data_, size, q);
    GPComponentCuda::mont_encode(dst, dst);
}
void GPComponentCuda::mul_scalar(GPComponentCuda& dst, const GPComponentCuda& src1, uint64_t src_scalar)
{
    assert(dst.like(src1));
    size_t size = dst.get_size();
    uint64_t encoded = dst.mm_.encode(src_scalar);
    cuda_batch_mul_scalar(dst.data_, src1.data_, encoded, size, dst.mm_);

}
void GPComponentCuda::mont_encode(GPComponentCuda& dst, const GPComponentCuda& src)
{
    assert(dst.like(src));
    mul_scalar(dst, src, dst.mm_.R2);

}
void GPComponentCuda::mont_decode(GPComponentCuda& dst, const GPComponentCuda& src)
{
    assert(dst.like(src));
    size_t size = dst.get_size();
    mul_scalar(dst, src, 1);
}
void GPComponentCuda::mul_mont(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2)
{
    assert(dst.like(src1));
    assert(dst.like(src2));
    size_t size = dst.get_size();
    uint64_t q = dst.q_;
    cuda_batch_mul_mont(dst.data_, src1.data_, src2.data_, size, q);
}
