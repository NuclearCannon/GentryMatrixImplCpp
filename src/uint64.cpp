#include "uint64.hpp"

u64 mod_mul(u64 a, u64 b, u64 mod) {
    __uint128_t product = (__uint128_t)a * b;
    return (u64)(product % mod);
}


u64 mod_pow(u64 base, u64 e, u64 mod)
{
    u64 result = 1;
    base = base % mod;
    while (e > 0) {
        if (e & 1) {
            result = mod_mul(result, base, mod);
        }
        base = mod_mul(base, base, mod);
        e >>= 1;
    }
    return result;
}

u64 mod_inv(u64 x, u64 mod)
{
    return mod_pow(x, mod-2, mod);
}
