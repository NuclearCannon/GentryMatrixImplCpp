#include "ntt.hpp"
#include "math_utils.hpp"
#include <cstring>


StandardNTTer::StandardNTTer(size_t n, uint64_t q, uint64_t qroot):
    n_(n), q_(q), mm(q)
{
    logn_ = Log2(n);
    assert(q%n==1);
    uint64_t nroot = mod_pow(qroot, (q-1)/n, q);
    roots_ = get_powers(nroot, n, q);
    iroots_ = get_powers(mod_inv(nroot, q), n, q);
    ninv_ = mod_inv(n, q);

    roots_mont_ = mm.batch_encode(roots_);
    iroots_mont_ = mm.batch_encode(iroots_);
    ninv_mont_ = mm.encode(ninv_);
}

StandardNTTer::~StandardNTTer() = default;

static void _butterfly_inc_mont(
    uint64_t* a, 
    const uint64_t* omegas_mont,
    size_t logn,
    const MontgomeryMultiplier& mm
)
{
    size_t m = 1;
    size_t t = logn;
    const size_t n = 1<<logn;
    const uint64_t q = mm.getM();
    while(m<n)
    {
        size_t m_half = m;
        m <<= 1;
        t--;
        for(size_t k=0; k<n; k+=m)
        {
            for(size_t j=0; j<m_half; j++)
            {
                uint64_t w = omegas_mont[(j<<t) & (n-1)];
                size_t i1 = k+j;
                size_t i2 = i1 + m_half;
                uint64_t u = a[i1];
                uint64_t v = mm.mul(a[i2],w);
                a[i1] = mod_add(u,v,q);
                a[i2] = mod_sub(u,v,q);
            }
        }
    }
}


static void _butterfly_dec_mont(
    uint64_t* a, 
    const uint64_t* omegas_mont,
    size_t logn,
    const MontgomeryMultiplier& mm
)
{
    const size_t n = 1<<logn;
    size_t m = n;
    size_t t = 0;
    const uint64_t q = mm.getM();
    // 恒有m*t == n/2
    while (m>1) {
        size_t m_half = m>>1;
        for(size_t k=0; k<n; k+=m)
        {
            for(size_t j=k; j<k+m_half; j++)
            {
                uint64_t w = omegas_mont[(j<<t) & (n-1)];
                size_t i1 = j;
                size_t i2 = j + m_half;
                uint64_t u = a[i1];
                uint64_t v = a[i2];
                a[i1] = mod_add(u,v,q);
                a[i2] = mm.mul(mod_sub(u,v,q),w);
            }
        }
        m = m_half;
        t++;
    }
}


void StandardNTTer::ntt(uint64_t* dst) const
{
    _butterfly_dec_mont(dst, roots_mont_.data(), logn_, mm);
}
void StandardNTTer::intt(uint64_t* dst) const
{
    _butterfly_inc_mont(dst, iroots_mont_.data(), logn_, mm);
    for(int i=0;i<n_;i++)dst[i] = mm.mul(dst[i], ninv_mont_);
}