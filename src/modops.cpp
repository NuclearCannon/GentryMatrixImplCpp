
#include <cassert>
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


void get_powers(vec64& dst, uint64_t x, size_t len, uint64_t mod)
{
    assert(dst.size() == len);
    dst[0] = 1;
    x %= mod;
    for(int i=1;i<len;i++)
    {
        dst[i] = mod_mul(dst[i-1], x, mod);
    }
}

vec64 get_powers(uint64_t x, size_t len, uint64_t mod)
{
    vec64 result(len);
    get_powers(result, x, len, mod);
    return result;
}



std::size_t vec64hash(const vec64& v) {
    std::size_t seed = 0;
    for (uint64_t x : v) {
        // 使用 boost 风格的 hash_combine 简化版
        seed += x + (seed << 6) + (seed >> 2);
    }
    return seed;
}


std::size_t vv64hash(const vv64& v) {
    std::size_t seed = 0;
    for (const auto& inner : v) {
        seed += vec64hash(inner) + (seed << 6) + (seed >> 2);
    }
    return seed;
}