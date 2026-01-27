#include "u64_array.hpp"
#include "ntt.hpp"


TwistedNtterW64::TwistedNtterW64(int p, u64 q, u64 eta):
    p_(p), q_(q), 
    eta_powers_(p), 
    p_inv_(mod_inv(p, q))
{
    // 检查eta合法性
    assert(mod_pow(eta, p, q) == 1);
    // 计算eta的各个次幂
    get_powers(eta_powers_, eta, p, q);
}

TwistedNtterW64::~TwistedNtterW64()
{
    // do nothing
}

void TwistedNtterW64::ntt(vec64& dst, const vec64& src) const
{
    if (&src == &dst)
    {
        vec64 buf(dst.size());
        _ntt_no_alias(buf, src);
        dst.swap(buf);
    }
    else
    {
        _ntt_no_alias(dst, src);
    }
}
void TwistedNtterW64::intt(vec64& dst, const vec64& src) const
{
    if (&src == &dst)
    {
        vec64 buf(dst.size());
        _intt_no_alias(buf, src);
        dst.swap(buf);
    }
    else
    {
        _intt_no_alias(dst, src);
    }
    
}


void TwistedNtterW64::_ntt_no_alias(vec64& dst, const vec64& src) const
{
    // 朴素实现
    assert(src.size() == p_-1);
    assert(dst.size() == p_-1);
    assert(&src != &dst);
    for(int k=1;k<p_;k++)
    {
        u64 Xk = 0;
        for(int j=0;j<p_-1;j++)
        {
            Xk += mod_mul(src[j], eta_powers_[(k*j)%p_], q_);
            Xk %= q_;
        }
        dst[k-1] = Xk;
    }
}
void TwistedNtterW64::_intt_no_alias(vec64& dst, const vec64& src) const
{
    assert(src.size() == p_-1);
    assert(dst.size() == p_-1);
    assert(&src != &dst);
    u64 sum_t = 0;
    for(int j=0;j<p_-1;j++)
    {
        u64 tj = 0;
        for(int k=1;k<p_;k++)
        {
            tj += mod_mul(src[k-1], eta_powers_[(k*(p_-j))%p_], q_);
            tj %= q_;
        }
        dst[j] = tj;
        sum_t = (sum_t + tj) % q_;
    }
    for(int j=0;j<p_-1;j++)
    {
        dst[j] = mod_mul(
            (dst[j]+sum_t) % q_, 
            p_inv_,
            q_
        );
    }
    
}