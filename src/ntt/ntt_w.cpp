#include "u64_array.hpp"
#include "ntt.hpp"
#include <cstring>
#include "GPU/cuda_u64_ctx_ops.hpp"
#include "modops.hpp"

TwistedNtterW64::TwistedNtterW64(int p , uint64_t q, uint64_t qroot):
    p_(p), q_(q), 
    subntter(p-1, q, qroot),
    mm(q),
    buf1(p-1), 
    b_(p-1), binv_(p-1),
    b_cuda_((p-1)*sizeof(uint64_t)), 
    b_inv_cuda_((p-1)*sizeof(uint64_t))
{
    vec64 gpp = get_powers(3, p-1, p);
    // 生成etas_
    uint64_t eta = mod_pow(qroot, (q-1)/p, q);
    for(int i=0; i<p-1; i++)b_[i] = mod_pow(eta, gpp[i], q);
    subntter.ntt(b_.data());
    // 生成逐位乘法逆元
    for(int i=0; i<p-1; i++)binv_[i] = mod_inv(b_[i], q);
    // 对b和binv都除以p-1，以中和subntter的intt副作用
    uint64_t inv = mod_inv(p-1, q);
    for(int i=0; i<p-1; i++)
    {
        b_[i] = mod_mul(b_[i], inv, q);
        binv_[i] = mod_mul(binv_[i], inv, q);
    }
    // 蒙哥马利编码
    for(auto& i: b_)i = mm.encode(i);
    for(auto& i: binv_)i = mm.encode(i);
    b_cuda_.copy_from_host(b_.data());
    b_inv_cuda_.copy_from_host(binv_.data());
}

TwistedNtterW64::~TwistedNtterW64()
{
    // do nothing
}


void TwistedNtterW64::ntt_mont(vec64& dst, const vec64& src) const
{
    memcpy(buf1.data(), src.data(), (p_-1)*sizeof(uint64_t));
    subntter.ntt(buf1.data());
    for(int i=0; i<p_-1; i++)dst[i] = mm.mul(buf1[i], b_[i]);
    subntter.intt(dst.data());    
}
void TwistedNtterW64::intt_mont(vec64& dst, const vec64& src) const
{
    memcpy(buf1.data(), src.data(), (p_-1)*sizeof(uint64_t));
    subntter.ntt(buf1.data());
    for(int i=0; i<p_-1; i++)dst[i] = mm.mul(buf1[i], binv_[i]);
    subntter.intt(dst.data());
}

void TwistedNtterW64::ntt_batch(uint64_t* dst, size_t batch_size) const
{
    subntter.ntt_batch(dst, batch_size);
    const size_t N = p_-1;
    for(int i=0; i<batch_size; i++)
    {
        for(int j=0; j<N; j++)
        {
            dst[i*N+j] = mm.mul(dst[i*N+j], b_[j]);
        }
    }
    subntter.intt_batch(dst, batch_size);
}
void TwistedNtterW64::intt_batch(uint64_t* dst, size_t batch_size) const
{
    subntter.ntt_batch(dst, batch_size);
    const size_t N = p_-1;
    for(int i=0; i<batch_size; i++)
    {
        for(int j=0; j<N; j++)
        {
            dst[i*N+j] = mm.mul(dst[i*N+j], binv_[j]);
        }
    }
    subntter.intt_batch(dst, batch_size);
}

void TwistedNtterW64::ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    assert(a.size() == (p_-1)*batch_size*sizeof(uint64_t));
    subntter.ntt_batch_cuda(a, batch_size);
    cuda_batch_mul_vec(a, a, b_cuda_, batch_size, p_-1, mm);
    subntter.intt_batch_cuda(a, batch_size);
}
void TwistedNtterW64::intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    assert(a.size() == (p_-1)*batch_size*sizeof(uint64_t));
    subntter.ntt_batch_cuda(a, batch_size);
    cuda_batch_mul_vec(a, a, b_inv_cuda_, batch_size, p_-1, mm);
    subntter.intt_batch_cuda(a, batch_size);
}