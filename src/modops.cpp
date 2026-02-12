#include <cstdint>
#include "modops.hpp"
#include "vec64.hpp"


uint64_t mod_pow(uint64_t base, uint64_t e, uint64_t mod)
{
    uint64_t result = 1;
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

uint64_t mod_inv(uint64_t x, uint64_t mod)
{
    return mod_pow(x, mod-2, mod);
}

