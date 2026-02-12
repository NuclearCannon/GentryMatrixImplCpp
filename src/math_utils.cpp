#include <cassert>
#include <cstddef>


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
