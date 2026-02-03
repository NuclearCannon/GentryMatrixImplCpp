#include "u64_array.hpp"
#include "ntt.hpp"


TwistedNtterW64::TwistedNtterW64(int p , u64 q, u64 eta):
    p_(p), q_(q), 
    rader(p, q, eta)
{
    // 检查eta合法性
    assert(mod_pow(eta, p, q) == 1);
}

TwistedNtterW64::~TwistedNtterW64()
{
    // do nothing
}

void TwistedNtterW64::ntt(vec64& dst, const vec64& src) const
{
    vec64 a2(src), A2(p_);
    a2.push_back(0);
    rader.rader(A2.data(), a2.data());
    for(int i=0;i<p_-1;i++)
    {
        dst[i] = A2[i+1];
    }
}
void TwistedNtterW64::intt(vec64& dst, const vec64& src) const
{
    vec64 A2(p_), a2(p_);
    for(int i=0;i<p_-1;i++)
    {
        A2[i+1] = src[i];
    }
    A2[0] = 0;
    rader.irader(a2.data(), A2.data());
    u64 delta = a2[p_-1];
    for(int i=0;i<p_-1;i++)
    {
        dst[i] = mod_sub(a2[i], delta, q_);
    }
}

