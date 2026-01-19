#include "ntt.hpp"
#include "flints.hpp"
#include <vector>
#include <cassert>


void ntt_standard_flint(const fmpz* a, fmpz* dst, fmpz_t root, size_t n, const fmpz_t q)
{
    // 输入校验
    assert(a != nullptr && dst != nullptr);
    assert((n & (n - 1)) == 0 && n >= 1);
    assert(a != dst);   // 暂不支持inplace NTT

    // Step 1: Bit-reverse copy
    const auto& rev = get_bit_reverse_table(n);


    // let dst[i] = a[rev[i]] % q
    for (size_t i = 0; i < n; ++i) {
        fmpz_mod(dst+i, a+rev[i], q);
    }

    // Step 2: Cooley-Tukey iterative NTT
    size_t m = 1;
    fmpz_t exponent, w_m, w, u, v, tmp1, tmp2;
    fmpz_init(exponent); fmpz_init(w_m); fmpz_init(w);
    fmpz_init(u); fmpz_init(v); fmpz_init(tmp1); fmpz_init(tmp2);

    while (m < n) {
        // exponent = n / (2 * m)
        fmpz_set_ui(exponent, n / (2 * m));
        // w_m = root^exponent mod q
        fmpz_powm(w_m, root, exponent, q);

        for (size_t i = 0; i < n; i += 2 * m) {
            fmpz_set_ui(w, 1);
            for (size_t j = i; j < i + m; ++j) {
                fmpz_set(u, &dst[j]);
                // v = dst[j + m] * w mod q
                fmpz_mul(v, &dst[j + m], w);
                fmpz_mod(v, v, q);

                // dst[j] = (u + v) mod q
                fmpz_add(tmp1, u, v);
                fmpz_mod(&dst[j], tmp1, q);

                // dst[j + m] = (u - v) mod q
                fmpz_sub(tmp2, u, v);
                fmpz_mod(&dst[j + m], tmp2, q);

                // w = w * w_m mod q
                fmpz_mul(w, w, w_m);
                fmpz_mod(w, w, q);
            }
        }
        m *= 2;
    }

    fmpz_clear(exponent); fmpz_clear(w_m); fmpz_clear(w);
    fmpz_clear(u); fmpz_clear(v); fmpz_clear(tmp1); fmpz_clear(tmp2);
}


void ntt_standard_flint(const fmpz_vector& a, fmpz_vector& dst, fmpz_t root, size_t n, const fmpz_t q)
{
    ntt_standard_flint(
        a.raw(), dst.raw(), root, n, q
    );
}