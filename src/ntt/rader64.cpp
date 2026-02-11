#include "ntt.hpp"



RaderNTTer64::RaderNTTer64(u64 p, u64 q, u64 qroot):
    p_(p), eta_(mod_pow(qroot, (q-1)/p, q)), q_(q), pinv_(mod_inv(p, q)), 
    gpp(p), gnp(p), b1ntt(p-1), b2ntt(p-1), subntter(p-1, q, qroot), mm(q),
    buf_a_(p-1), buf_c_(p-1), buf_a_ntt_(p-1), buf_c_ntt_(p-1)
{
    assert(is_power_of_two(p-1));
    // 计算g的各个次幂
    get_powers(gpp, 3, p, p);
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
    subntter.ntt_mont(b1ntt.data(), b1.data());
    subntter.ntt_mont(b2ntt.data(), b2.data());
    
    b1ntt_mont = mm.batch_encode(b1ntt);
    b2ntt_mont = mm.batch_encode(b2ntt);
}

RaderNTTer64::~RaderNTTer64() = default;



void RaderNTTer64::_rader_inner_mont(u64* dst, const u64* src, const vec64& bntt) const
{
    vec64& a = buf_a_, &c = buf_c_, &antt = buf_a_ntt_, &cntt = buf_c_ntt_;
    // a[i] = x[gnp[i]]
    for(int i=0; i<p_-1; i++)
    {
        a[i] = src[gnp[i]];
    }
    u64 x0 = src[0];
    subntter.ntt_mont(antt.data(), a.data());
    mm.vec_mul_mont(cntt, antt, bntt);
    subntter.intt_mont(c.data(), cntt.data());
    for(int u=0; u<p_-1; u++)
    {
        dst[gpp[u]] = mod_add(x0, c[u], q_);
    }
    u64 sumx = 0;
    for(int i=0; i<p_; i++)sumx = mod_add(sumx, src[i], q_);
    dst[0] = sumx;

}



void RaderNTTer64::rader_mont(u64* dst, const u64* src) const
{
    _rader_inner_mont(dst, src, b1ntt_mont);
}
void RaderNTTer64::irader_mont(u64* dst, const u64* src) const
{
    _rader_inner_mont(dst, src, b2ntt_mont);
    // 除以p
    mm.vec_scalar_mul_mont_ptr(dst, dst, p_, pinv_);
}