#include <cassert>
#include <cstddef>
#include "math_utils.hpp"
#include <bits/stl_pair.h>


bool is_power_of_two(size_t x) noexcept {
    return ((x>0) && ((x & (x - 1)) == 0));
}

size_t Log2(size_t x)
{
    assert(is_power_of_two(x));
    size_t i=-1;
    while(x)
    {
        x>>=1;
        i++;
    }
    return i;
}

uint64_t modinv64(uint64_t a) {
    // 要求 a 是奇数，否则无逆元
    if ((a & 1) == 0) {
        return 0; 
    }

    uint64_t x = 1; // a ≡ 1 (mod 2), so inverse is 1 mod 2
    // 牛顿迭代：每次加倍精度
    x = x * (2 - a * x); // 2 bits
    x = x * (2 - a * x); // 4 bits
    x = x * (2 - a * x); // 8 bits
    x = x * (2 - a * x); // 16 bits
    x = x * (2 - a * x); // 32 bits
    x = x * (2 - a * x); // 64 bits

    return x;
}
void transpose_auto(uint64_t* dst, const uint64_t* src, int n)
{
    if(dst == src)transpose_inplace(dst, n);
    else transpose_restrict(dst, src, n);
}

void transpose_restrict(uint64_t* __restrict__ dst, const uint64_t* __restrict__ src, int n)
{
    assert(dst != src);
    for(int x=0;x<n;x++)
    {
        for(int y=0;y<n;y++)
        {
            dst[x*n+y] = src[y*n+x];
        }
    }
}

void transpose_inplace(uint64_t* dst, int n)
{
    for(int x=0;x<n;x++)
    {
        for(int y=x+1;y<n;y++)
        {
            std::swap(dst[x*n+y], dst[y*n+x]);
        }
    }
}