#include "ntt.hpp"
#include "math_utils.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"
#include <cstring>
#include "modops.hpp"

TwistedNtterXY::TwistedNtterXY(int n, uint64_t q, uint64_t qroot):
    n_(n), q_(q),
    zeta_pos_pows_mont_(n),
    zeta_neg_pows_mont_(n),
    std_ntter_(n, q, qroot),
    mm_(q),
    zeta_pos_pows_mont_cuda_(n*sizeof(uint64_t)),
    zeta_neg_pows_mont_cuda_(n*sizeof(uint64_t))
{
    // 检查n：必须是power of 2
    assert(is_power_of_two(n));
    assert(q%(4*n)==1);
    uint64_t zeta = mod_pow(qroot, (q-1)/(4*n), q);
    // 检查zeta合法性
    uint64_t t = mod_pow(zeta, n*2, q);  // 理论上应该是-1
    assert(t+1==q);
    // zeta合法
    // 生成 zeta_pos_pows_
    vec64 zeta_pos_pows = get_powers(zeta, n, q);
    // 生成 zeta_neg_pows_ninv_
    uint64_t zeta_inv = mod_inv(zeta, q);
    vec64 zeta_neg_pows = get_powers(zeta_inv, n, q);

    uint64_t ninv = mod_inv(n, q);
    for(auto& i:zeta_neg_pows)
    {
        i = mod_mul(i,ninv,q);
    }
    for(int i=0; i<n; i++)
    {
        zeta_pos_pows_mont_[i] = mm_.encode(zeta_pos_pows[i]);
        zeta_neg_pows_mont_[i] = mm_.encode(zeta_neg_pows[i]);
    }

    zeta_pos_pows_mont_cuda_.copy_from_host(zeta_pos_pows_mont_.data());
    zeta_neg_pows_mont_cuda_.copy_from_host(zeta_neg_pows_mont_.data());



}

TwistedNtterXY::~TwistedNtterXY()
{
    // do nothing
}

void TwistedNtterXY::ntt_mont(vec64& dst, const vec64& src) const
{
    for(int i=0; i<n_; i++)dst[i] = mm_.mul(src[i], zeta_pos_pows_mont_[i]);
    std_ntter_.ntt(dst.data());
}
void TwistedNtterXY::intt_mont(vec64& dst, const vec64& src) const
{
    if(&src != &dst)memcpy(dst.data(), src.data(), n_*sizeof(uint64_t));
    std_ntter_.intt(dst.data());
    for(int i=0; i<n_; i++)dst[i] = mm_.mul(dst[i], zeta_neg_pows_mont_[i]);
}

void TwistedNtterXY::ntt_batch(uint64_t* dst, size_t batch_size) const
{
    for(size_t i=0; i<batch_size; i++)
    {
        for(int j=0; j<n_; j++)
        {
            dst[i*n_ + j] = mm_.mul(dst[i*n_ + j], zeta_pos_pows_mont_[j]);
        }
    }
    std_ntter_.ntt_batch(dst, batch_size);
}
void TwistedNtterXY::intt_batch(uint64_t* dst, size_t batch_size) const
{
    std_ntter_.intt_batch(dst, batch_size);
    for(size_t i=0; i<batch_size; i++)
    {
        for(int j=0; j<n_; j++)
        {
            dst[i*n_ + j] = mm_.mul(dst[i*n_ + j], zeta_neg_pows_mont_[j]);
        }
    }

}

void TwistedNtterXY::ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    cuda_batch_mul_vec(
        a, a,
        zeta_pos_pows_mont_cuda_,
        batch_size,
        n_,
        mm_
    );

    std_ntter_.ntt_batch_cuda(a, batch_size);
}
void TwistedNtterXY::intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    std_ntter_.intt_batch_cuda(a, batch_size);
    cuda_batch_mul_vec(
        a, a,
        zeta_neg_pows_mont_cuda_,
        batch_size,
        n_,
        mm_
    );

}