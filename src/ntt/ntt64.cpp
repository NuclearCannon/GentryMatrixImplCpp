#include "ntt.hpp"
#include "uint64.hpp"


void ntt_standard_64(
    const u64* a, 
    u64* dst, 
    u64 root,
    size_t n, 
    const u64 mod
)
{
    std::vector<u64> roots(n);
    roots[0] = 1;
    for(int i=1;i<n;i++)roots[i] = mod_mul(roots[i-1], root, mod);
    ntt_standard_64_with_roots(a,dst,roots.data(), n, mod);
}


void ntt_standard_64_with_roots(
    const u64* a, 
    u64* dst, 
    const u64* roots,   // 需要提供root的至少[0,n/2)次方 
    size_t n, 
    const u64 mod
)
{
    assert(a != nullptr && dst != nullptr);
    assert((n & (n - 1)) == 0 && n >= 1);
    assert(a != dst);   // 暂不支持……
    const auto& rev = get_bit_reverse_table(n);
    // dst[i] = a[rev[i]] mod ctx
    for (size_t i = 0; i < n; ++i) {
        dst[i] = a[rev[i]] % mod;
    }
    size_t m = 1;
    int t = log2(n>>1);
    while (t>=0) {
        for (size_t i = 0; i < n; i += 2 * m) {
            for (size_t k=0; k < m; ++k) {
                size_t j = i+k;
                // printf("k=%zu, t=%d, k<<t=%zu\n", k, t, k<<t);
                u64 w = roots[k<<t];
                u64 u = dst[j];
                u64 v = mod_mul(dst[j+m], w, mod);
                dst[j] = mod_add(u, v, mod);
                dst[j+m] = mod_sub(u, v, mod);         
            }
        }
        m <<= 1;    // m *= 2
        t--;
    }
}