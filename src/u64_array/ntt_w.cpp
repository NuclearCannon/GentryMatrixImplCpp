#include "u64_array.hpp"
#include "ntt.hpp"
#include <cstring>


TwistedNtterW64::TwistedNtterW64(int p , u64 q, u64 qroot):
    p_(p), q_(q), 
    subntter(p-1, q, qroot),
    mm(q),
    buf1(p-1), 
    buf2(p-1),
    b_(p-1), binv_(p-1)
{
    vec64 gpp = get_powers(3, p-1, p);
    // 生成etas_
    u64 eta = mod_pow(qroot, (q-1)/p, q);
    vec64 etas(p-1);
    for(int i=0; i<p-1; i++)etas[i] = mod_pow(eta, gpp[i], q);
    subntter.ntt_mont(b_.data(), etas.data());
    // 生成逐位乘法逆元
    for(int i=0; i<p-1; i++)binv_[i] = mod_inv(b_[i], q);

    // 蒙哥马利编码
    mm.batch_encode_inplace(b_);
    mm.batch_encode_inplace(binv_);

}

TwistedNtterW64::~TwistedNtterW64()
{
    // do nothing
}


void TwistedNtterW64::ntt_mont(vec64& dst, const vec64& src) const
{
    vec64& nttg = buf1, &nttG = buf2;
    subntter.ntt_mont(nttg.data(), src.data());
    mm.vec_mul_mont(nttG, nttg, b_);
    subntter.intt_mont(dst.data(), nttG.data());

    
}
void TwistedNtterW64::intt_mont(vec64& dst, const vec64& src) const
{
    vec64& nttg = buf1, &nttG = buf2;
    subntter.ntt_mont(nttG.data(), src.data());
    mm.vec_mul_mont(nttg, nttG, binv_);
    subntter.intt_mont(dst.data(), nttg.data());
}
