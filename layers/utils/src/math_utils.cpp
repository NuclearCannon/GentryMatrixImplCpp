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

void transpose_auto(uint64_t* dst, const uint64_t* src, int n)
{
    if(dst == src)transpose_inplace(dst, n);
    else transpose_restrict(dst, src, n);
}
