#include "u64_array.hpp"
#include "ntt.hpp"
#include <cstring>


TwistedNtterW64::TwistedNtterW64(int p , uint64_t q, uint64_t qroot):
    p_(p), q_(q), 
    subntter(p-1, q, qroot),
    mm(q),
    buf1(p-1), 
    b_(p-1), binv_(p-1)
{
    vec64 gpp = get_powers(3, p-1, p);
    // 生成etas_
    uint64_t eta = mod_pow(qroot, (q-1)/p, q);
    for(int i=0; i<p-1; i++)b_[i] = mod_pow(eta, gpp[i], q);
    subntter.ntt(b_.data());
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
    memcpy(buf1.data(), src.data(), (p_-1)*sizeof(uint64_t));
    subntter.ntt(buf1.data());
    mm.vec_mul_mont(dst, buf1, b_);
    subntter.intt(dst.data());    
}
void TwistedNtterW64::intt_mont(vec64& dst, const vec64& src) const
{
    memcpy(buf1.data(), src.data(), (p_-1)*sizeof(uint64_t));
    subntter.ntt(buf1.data());
    mm.vec_mul_mont(dst, buf1, binv_);
    subntter.intt(dst.data());
}
