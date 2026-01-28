#include "uint64.hpp"
#include <cassert>

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

void vec_mul(vec64& dst, const vec64& src1, const vec64& src2, u64 mod)
{
    std::size_t len = dst.size();
    assert(src1.size() == len);
    assert(src2.size() == len);
    assert(mod);
    for(std::size_t i=0;i<len;i++)
    {
        dst[i] = mod_mul(src1[i], src2[i], mod);
    }
}

// 逐位乘
void vec_scalar_mul(vec64& dst, const vec64& src1, u64 src2, u64 mod)
{
    std::size_t len = dst.size();
    assert(src1.size() == len);
    assert(mod);
    for(std::size_t i=0;i<len;i++)
    {
        dst[i] = mod_mul(src1[i], src2, mod);
    }
}

void get_powers(vec64& dst, u64 x, size_t len, u64 mod)
{
    assert(dst.size() == len);
    dst[0] = 1;
    for(int i=1;i<len;i++)
    {
        dst[i] = mod_mul(dst[i-1], x, mod);
    }
}


std::vector<u64> copy_from(const std::vector<u64>& src, size_t begin, size_t length)
{
    // 边界检查：确保 begin 不越界
    assert(begin + length <= src.size());
    // 使用迭代器构造新 vector
    return std::vector<u64>(src.begin() + begin, src.begin() + begin + length);
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