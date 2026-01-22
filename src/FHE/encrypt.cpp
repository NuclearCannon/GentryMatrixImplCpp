#include "FHE/encrypt.hpp"
#include <cassert>



std::pair<ZiqArray, ZiqArray> encrypt(const ZiqArray& message, const ZiqArray& sk)
{
    
    assert (message.ctx() == sk.ctx());
    const ZiqArrayContext* ctx = message.ctx();
    // 随机生成一个a
    ZiqArray a = ctx->uniform();
    // 随机生成一个e
    ZiqArray e = ctx->dg();
    // b = message + e - a*s
    ZiqArray a_ntt = a.iw_ntt().xy_ntt();
    ZiqArray s_ntt = sk.iw_ntt().xy_ntt();
    ZiqArray as_ntt = a_ntt.mul(s_ntt);
    ZiqArray as = as_ntt.iw_intt().xy_intt();
    ZiqArray b = message.add(e).sub(as);
    // ZiqArray b = message.sub(as);
    return std::make_pair(std::move(a), std::move(b));
}


ZiqArray decrypt(const ZiqArray& ct_a, const ZiqArray& ct_b, const ZiqArray& sk)
{
    // return a*s+b
    ZiqArray a_ntt = ct_a.iw_ntt().xy_ntt();
    ZiqArray s_ntt = sk.iw_ntt().xy_ntt();
    ZiqArray as_ntt = a_ntt.mul(s_ntt);
    ZiqArray as = as_ntt.iw_intt().xy_intt();
    return as.add(ct_b);
}