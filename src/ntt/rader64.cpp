#include "ntt.hpp"



RaderNTTer64::RaderNTTer64(u64 p, u64 q, u64 qroot):
    p_(p), g_(findPrimitiveRoot(p)), eta_(mod_pow(qroot, (q-1)/p, q)), q_(q), pinv_(mod_inv(p, q)), 
    gpp(p), gnp(p), b1ntt(p-1), b2ntt(p-1), subntter(p-1, q, qroot)
{

    // 计算g的各个次幂
    get_powers(gpp, g_, p, p);
    assert(gpp[p-1] == 1);
    for(int i=0;i<p;i++)gnp[i] = gpp[p-1-i];
    // 计算eta^{g^i}
    vec64 b1(p-1), b2(p-1);
    u64 ieta = mod_inv(eta_, q);
    for(int i=0;i<p-1;i++)
    {
        b1[i] = mod_pow(eta_,  gpp[i], q);
        b2[i] = mod_pow(ieta, gpp[i], q);
    }
    // 预先NTT
    subntter.ntt(b1ntt.data(), b1.data());
    subntter.ntt(b2ntt.data(), b2.data());
}

RaderNTTer64::~RaderNTTer64() = default;

void RaderNTTer64::_rader_inner(u64* dst, const u64* src, const vec64& bntt) const
{
    vec64 a(p_-1), c(p_-1), cntt(p_-1);
    // a[i] = x[gnp[i]]
    for(int i=0; i<p_-1; i++)
    {
        a[i] = src[gnp[i]];
    }
    u64 x0 = src[0];
    vec64 antt(p_-1);
    subntter.ntt(antt.data(), a.data());
    vec_mul(cntt, antt, bntt, q_);
    subntter.intt(c.data(), cntt.data());
    for(int u=0; u<p_-1; u++)
    {
        dst[gpp[u]] = mod_add(x0, c[u], q_);
    }
    u64 sumx = 0;
    for(int i=0; i<p_; i++)sumx = mod_add(sumx, src[i], q_);
    dst[0] = sumx;

}

void RaderNTTer64::rader(u64* dst, const u64* src) const
{
    _rader_inner(dst, src, b1ntt);
}
void RaderNTTer64::irader(u64* dst, const u64* src) const
{
    _rader_inner(dst, src, b2ntt);
    // 除以p
    for(int i=0; i<p_; i++)dst[i] = mod_mul(dst[i], pinv_, q_);
}