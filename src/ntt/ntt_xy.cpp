#include "u64_array.hpp"
#include "ntt.hpp"
#include "math_utils.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"
#include <cstring>
#include "modops.hpp"

TwistedNtterXY64::TwistedNtterXY64(int n, uint64_t q, uint64_t qroot):
    n_(n), q_(q),
    zeta_pos_pows_(n),
    zeta_neg_pows_(n),
    zeta_pos_pows_mont_(n),
    zeta_neg_pows_mont_(n),
    std_ntter(n, q, qroot),
    mm(q),
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
    get_powers(zeta_pos_pows_, zeta, n, q);
    // 生成 zeta_neg_pows_ninv_
    uint64_t zeta_inv = mod_inv(zeta, q);
    get_powers(zeta_neg_pows_, zeta_inv, n, q);

    uint64_t ninv = mod_inv(n, q);
    for(auto& i:zeta_neg_pows_)
    {
        i = mod_mul(i,ninv,q);
    }
    for(int i=0; i<n; i++)
    {
        zeta_pos_pows_mont_[i] = mm.encode(zeta_pos_pows_[i]);
        zeta_neg_pows_mont_[i] = mm.encode(zeta_neg_pows_[i]);
    }

    zeta_pos_pows_mont_cuda_.copy_from_host(zeta_pos_pows_mont_.data());
    zeta_neg_pows_mont_cuda_.copy_from_host(zeta_neg_pows_mont_.data());



}

TwistedNtterXY64::~TwistedNtterXY64()
{
    // do nothing
}

void TwistedNtterXY64::ntt_mont(vec64& dst, const vec64& src) const
{
    for(int i=0; i<n_; i++)dst[i] = mm.mul(src[i], zeta_pos_pows_mont_[i]);
    std_ntter.ntt(dst.data());
}
void TwistedNtterXY64::intt_mont(vec64& dst, const vec64& src) const
{
    if(&src != &dst)memcpy(dst.data(), src.data(), n_*sizeof(uint64_t));
    std_ntter.intt(dst.data());
    for(int i=0; i<n_; i++)dst[i] = mm.mul(dst[i], zeta_neg_pows_mont_[i]);
}

void TwistedNtterXY64::ntt_batch(uint64_t* dst, size_t batch_size) const
{
    for(size_t i=0; i<batch_size; i++)
    {
        for(int j=0; j<n_; j++)
        {
            dst[i*n_ + j] = mm.mul(dst[i*n_ + j], zeta_pos_pows_mont_[j]);
        }
    }
    std_ntter.ntt_batch(dst, batch_size);
}
void TwistedNtterXY64::intt_batch(uint64_t* dst, size_t batch_size) const
{
    std_ntter.intt_batch(dst, batch_size);
    for(size_t i=0; i<batch_size; i++)
    {
        for(int j=0; j<n_; j++)
        {
            dst[i*n_ + j] = mm.mul(dst[i*n_ + j], zeta_neg_pows_mont_[j]);
        }
    }

}

void TwistedNtterXY64::ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    cuda_batch_mul_vec(
        a, a,
        zeta_pos_pows_mont_cuda_,
        batch_size,
        n_,
        mm
    );

    std_ntter.ntt_batch_cuda(a, batch_size);
}
void TwistedNtterXY64::intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const
{
    std_ntter.intt_batch_cuda(a, batch_size);
    cuda_batch_mul_vec(
        a, a,
        zeta_neg_pows_mont_cuda_,
        batch_size,
        n_,
        mm
    );

}