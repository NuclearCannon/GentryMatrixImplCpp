#include "u64_array.hpp"
#include "ntt.hpp"
#include <cstring>
#include "GPU/cuda_u64_ctx_ops.hpp"
#include "modops.hpp"

TwistedNtterW::TwistedNtterW(int p , uint64_t q, uint64_t qroot):
    p_(p), q_(q), 
    std_ntter_(p-1, q, qroot),
    mm_(q),
    b_(p-1), binv_(p-1),
    b_cuda_((p-1)*sizeof(uint64_t)), 
    b_inv_cuda_((p-1)*sizeof(uint64_t))
{
    vec64 gpp = get_powers(3, p-1, p);
    // 生成etas_
    uint64_t eta = mod_pow(qroot, (q-1)/p, q);
    for(int i=0; i<p-1; i++)b_[i] = mod_pow(eta, gpp[i], q);
    std_ntter_.ntt(b_.data());
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
    for(auto& i: b_)i = mm_.encode(i);
    for(auto& i: binv_)i = mm_.encode(i);
    b_cuda_.copy_from_host(b_.data());
    b_inv_cuda_.copy_from_host(binv_.data());
}

TwistedNtterW::~TwistedNtterW()
{
    // do nothing
}


void TwistedNtterW::ntt_mont(vec64& dst, const vec64& src) const
{
    if(&dst != &src)memcpy(dst.data(), src.data(), (p_-1)*sizeof(uint64_t));
    std_ntter_.ntt(dst.data());
    for(int i=0; i<p_-1; i++)dst[i] = mm_.mul(dst[i], b_[i]);
    std_ntter_.intt(dst.data());    
}
void TwistedNtterW::intt_mont(vec64& dst, const vec64& src) const
{
    if(&dst != &src)memcpy(dst.data(), src.data(), (p_-1)*sizeof(uint64_t));
    std_ntter_.ntt(dst.data());
    for(int i=0; i<p_-1; i++)dst[i] = mm_.mul(dst[i], binv_[i]);
    std_ntter_.intt(dst.data());
}

void TwistedNtterW::ntt_batch(uint64_t* dst, size_t batch_size) const
{
    std_ntter_.ntt_batch(dst, batch_size);
    const size_t N = p_-1;
    for(int i=0; i<batch_size; i++)
    {
        for(int j=0; j<N; j++)
        {
            dst[i*N+j] = mm_.mul(dst[i*N+j], b_[j]);
        }
    }
    std_ntter_.intt_batch(dst, batch_size);
}
void TwistedNtterW::intt_batch(uint64_t* dst, size_t batch_size) const
{
    std_ntter_.ntt_batch(dst, batch_size);
    const size_t N = p_-1;
    for(int i=0; i<batch_size; i++)
    {
        for(int j=0; j<N; j++)
        {
            dst[i*N+j] = mm_.mul(dst[i*N+j], binv_[j]);
        }
    }
    std_ntter_.intt_batch(dst, batch_size);
}

void TwistedNtterW::ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    assert(a.size() == (p_-1)*batch_size*sizeof(uint64_t));
    std_ntter_.ntt_batch_cuda(a, batch_size);
    cuda_batch_mul_vec(a, a, b_cuda_, batch_size, p_-1, mm_);
    std_ntter_.intt_batch_cuda(a, batch_size);
}
void TwistedNtterW::intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    assert(a.size() == (p_-1)*batch_size*sizeof(uint64_t));
    std_ntter_.ntt_batch_cuda(a, batch_size);
    cuda_batch_mul_vec(a, a, b_inv_cuda_, batch_size, p_-1, mm_);
    std_ntter_.intt_batch_cuda(a, batch_size);
}