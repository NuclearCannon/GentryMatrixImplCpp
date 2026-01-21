#include "ntt.hpp"
#include "flints.hpp"
#include <vector>
#include <cassert>



// 新增: 使用 fmpz_mod_ctx_t 取模的 NTT 实现
void ntt_standard_flint(const fmpz* a, fmpz* dst, const fmpz_t root, size_t n, const fmpz_mod_ctx_t ctx)
{
    assert(a != nullptr && dst != nullptr);
    assert((n & (n - 1)) == 0 && n >= 1);
    assert(a != dst);   // 暂不支持……

    const auto& rev = get_bit_reverse_table(n);

    // dst[i] = a[rev[i]] mod ctx
    for (size_t i = 0; i < n; ++i) {
        fmpz_mod_set_fmpz(dst+i, a+rev[i], ctx);
    }

    size_t m = 1;
    fmpz_t exponent, w_m, w, u, v;
    fmpz_init(exponent); fmpz_init(w_m); fmpz_init(w);
    fmpz_init(u); fmpz_init(v);

    while (m < n) {
        fmpz_set_ui(exponent, n / (2 * m));
        fmpz_mod_pow_fmpz(w_m, root, exponent, ctx);
        for (size_t i = 0; i < n; i += 2 * m) {
            fmpz_set_ui(w, 1);
            for (size_t j = i; j < i + m; ++j) {
                /*
                u = dst[j]
                v = dst[j+m]*w
                dst[j] = (u + v) mod q
                dst[j + m] = (u - v) mod q
                w = w * w_m mod q
                但是我们没必要真的把u读出来
                */
                fmpz* dst_j = &dst[j];
                fmpz* dst_j_m = &dst[j + m];
                // let v = dst[j+m]*w
                fmpz_mod_mul(v, dst_j_m, w, ctx);
                // let dst[j + m] = dst[j] - v
                fmpz_mod_sub(dst_j_m, dst_j, v, ctx);
                // let dst[j] += v
                fmpz_mod_add(dst_j, dst_j, v, ctx);                
                // w *= w_m
                fmpz_mod_mul(w, w, w_m, ctx);
            }
        }
        m *= 2;
    }

    fmpz_clear(exponent); fmpz_clear(w_m); fmpz_clear(w);
    fmpz_clear(u); fmpz_clear(v);
}


// 新增: fmpz_vector 版本
void ntt_standard_flint(const fmpz_vector& a, fmpz_vector& dst, const fmpz_t root, size_t n, const fmpz_mod_ctx_t ctx)
{
    ntt_standard_flint(
        a.raw(), dst.raw(), root, n, ctx
    );
}