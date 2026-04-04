#include "math_utils.hpp"

void transpose_rect_restrict(
    uint64_t* __restrict__ dst,
    const uint64_t* __restrict__ src,
    size_t r, size_t c
)
{
    assert(src!=dst);
    
    constexpr int TILE_SIZE = 8;
    assert(r % TILE_SIZE == 0);
    assert(c % TILE_SIZE == 0);
    // 分块转置
    for (size_t ii = 0; ii < r; ii += TILE_SIZE) {
        for (size_t jj = 0; jj < c; jj += TILE_SIZE) {
            size_t imax = ii + TILE_SIZE;
            size_t jmax = jj + TILE_SIZE;

            for (size_t i = ii; i < imax; ++i) {
                for (size_t j = jj; j < jmax; ++j) {
                    dst[j * r + i] = src[i * c + j];
                }
            }
        }
    }
}