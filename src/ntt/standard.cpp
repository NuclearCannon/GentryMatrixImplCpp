#include "ntt.hpp"
#include "math_utils.hpp"


StandardNTTer::StandardNTTer(size_t n, uint64_t q, uint64_t qroot):
    n_(n), q_(q), mm(q)
{
    logn_ = Log2(n);
    assert(q%n==1);
    uint64_t nroot = mod_pow(qroot, (q-1)/n, q);
    roots_ = get_powers(nroot, n/2, q);
    iroots_ = get_powers(mod_inv(nroot, q), n/2, q);
    ninv_ = mod_inv(n, q);

    roots_mont_ = mm.batch_encode(roots_);
    iroots_mont_ = mm.batch_encode(iroots_);
    ninv_mont_ = mm.encode(ninv_);
}

StandardNTTer::~StandardNTTer() = default;


void StandardNTTer::_ntt_standard_inner_mont(uint64_t* dst, const uint64_t* src, const uint64_t* roots) const
{
    assert(src != nullptr && dst != nullptr);
    assert(src != dst);   // 暂不支持……
    const auto& rev = get_bit_reverse_table_by_logn(logn_);
    // dst[i] = a[rev[i]] mod ctx
    for (size_t i = 0; i < n_; ++i) {
        dst[i] = src[rev[i]];
    }
    size_t m = 1;
    uint64_t t = n_>>1;
    // 恒有m*t == n/2
    while (t) {
        for (size_t i = 0; i < n_; i += 2 * m) {
            size_t end = i+m;
            for (size_t k=0, j=i; j < end; k+=t, ++j) {
                uint64_t w = roots[k];
                uint64_t u = dst[j];
                uint64_t v = mm.mul(dst[j+m], w);
                dst[j] = mod_add(u, v, q_);
                dst[j+m] = mod_sub(u, v, q_);         
            }
        }
        m <<= 1;    // m *= 2
        t >>= 1;
    }
}

void StandardNTTer::ntt_mont(uint64_t* dst, const uint64_t* src) const
{
    _ntt_standard_inner_mont(dst, src, roots_mont_.data());
}
void StandardNTTer::intt_mont(uint64_t* dst, const uint64_t* src) const
{
    _ntt_standard_inner_mont(dst, src, iroots_mont_.data());
    for(int i=0;i<n_;i++)dst[i] = mm.mul(dst[i], ninv_mont_);
}