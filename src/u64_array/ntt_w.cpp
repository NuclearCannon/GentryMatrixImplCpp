#include "u64_array.hpp"
#include "ntt.hpp"
#include <cstring>


TwistedNtterW64::TwistedNtterW64(int p , u64 q, u64 qroot):
    p_(p), q_(q), 
    rader(p, q, qroot),
    mm(q),
    buf1(std::make_unique<vec64>(p)), 
    buf2(std::make_unique<vec64>(p))
{

}

TwistedNtterW64::~TwistedNtterW64()
{
    // do nothing
}

void TwistedNtterW64::ntt(vec64& dst, const vec64& src) const
{
    vec64& a = *buf1, &A = *buf2;
    memcpy(a.data(), src.data(), (p_-1)*sizeof(u64));
    a[p_-1] = 0;
    mm.batch_encode_inplace(a);
    rader.rader_mont(A.data(), a.data());
    mm.batch_decode_inplace(A);
    memcpy(dst.data(), A.data()+1, (p_-1)*sizeof(u64));
}
void TwistedNtterW64::intt(vec64& dst, const vec64& src) const
{
    vec64& a = *buf1, &A = *buf2;
    memcpy(A.data()+1, src.data(), (p_-1)*sizeof(u64));
    A[0] = 0;
    mm.batch_encode_inplace(A);
    rader.irader_mont(a.data(), A.data());
    mm.batch_decode_inplace(a);
    u64 delta = mod_sub(0, a[p_-1], q_);
    for(int i=0;i<p_-1;i++)
    {
        dst[i] = mod_add(a[i], delta, q_);
    }
}

void TwistedNtterW64::ntt_mont(vec64& dst, const vec64& src) const
{
    vec64& a = *buf1, &A = *buf2;
    memcpy(a.data(), src.data(), (p_-1)*sizeof(u64));
    a[p_-1] = 0;
    rader.rader_mont(A.data(), a.data());
    memcpy(dst.data(), A.data()+1, (p_-1)*sizeof(u64));
}
void TwistedNtterW64::intt_mont(vec64& dst, const vec64& src) const
{
    vec64& a = *buf1, &A = *buf2;
    memcpy(A.data()+1, src.data(), (p_-1)*sizeof(u64));
    A[0] = 0;
    rader.irader_mont(a.data(), A.data());
    u64 delta = mod_sub(0, a[p_-1], q_);
    for(int i=0;i<p_-1;i++)
    {
        dst[i] = mod_add(a[i], delta, q_);
    }
}
