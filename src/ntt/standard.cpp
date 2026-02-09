#include "ntt.hpp"



StandardNTTer::StandardNTTer(size_t n, u64 q, u64 qroot):
    n_(n), q_(q), qroot_(qroot)
{
    logn_ = log2(n);
    assert(q%n==1);
    nroot_ = mod_pow(qroot, (q-1)/n, q);
    roots_ = get_powers(nroot_, n/2, q);
    iroots_ = get_powers(mod_inv(nroot_, q), n/2, q);
    ninv_ = mod_inv(n, q);
}

StandardNTTer::~StandardNTTer() = default;

void StandardNTTer::_ntt_standard_inner(u64* dst, const u64* src, const u64* roots) const
{
    assert(src != nullptr && dst != nullptr);
    assert(src != dst);   // 暂不支持……
    const auto& rev = get_bit_reverse_table_by_logn(logn_);
    // dst[i] = a[rev[i]] mod ctx
    for (size_t i = 0; i < n_; ++i) {
        dst[i] = src[rev[i]];
    }
    size_t m = 1;
    int t = logn_-1;
    while (t>=0) {
        for (size_t i = 0; i < n_; i += 2 * m) {
            for (size_t k=0; k < m; ++k) {
                size_t j = i+k;
                // printf("k=%zu, t=%d, k<<t=%zu\n", k, t, k<<t);
                u64 w = roots[k<<t];
                u64 u = dst[j];
                u64 v = mod_mul(dst[j+m], w, q_);
                dst[j] = mod_add(u, v, q_);
                dst[j+m] = mod_sub(u, v, q_);         
            }
        }
        m <<= 1;    // m *= 2
        t--;
    }
}

void StandardNTTer::ntt(u64* dst, const u64* src) const
{
    _ntt_standard_inner(dst, src, roots_.data());
}
void StandardNTTer::intt(u64* dst, const u64* src) const
{
    _ntt_standard_inner(dst, src, iroots_.data());
    for(int i=0;i<n_;i++)dst[i] = mod_mul(dst[i], ninv_, q_);
}