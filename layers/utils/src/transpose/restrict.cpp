#include "math_utils.hpp"

void transpose_restrict(uint64_t* __restrict__ dst,
                        const uint64_t* __restrict__ src,
                        int n)
{
    assert(dst != src);
    assert((n & (n - 1)) == 0); // 确保 n 是 power of 2
    constexpr int TILE = 8; // 必须是2的幂，且 <= n
    uint64_t tile[TILE][TILE];
    for (int ii = 0; ii < n; ii += TILE) {
        for (int jj = 0; jj < n; jj += TILE) {
            // 从 src 读取 [ii, ii+TILE) × [jj, jj+TILE) 块（注意：src[y][x]）
            for (int i = 0; i < TILE && ii + i < n; ++i) {
                for (int j = 0; j < TILE && jj + j < n; ++j) {
                    tile[i][j] = src[(jj + j) * n + (ii + i)]; // src[y][x]
                }
            }
            for (int i = 0; i < TILE && ii + i < n; ++i) {
                for (int j = 0; j < TILE && jj + j < n; ++j) {
                    dst[(ii + i) * n + (jj + j)] = tile[i][j];
                }
            }
        }
    }
}



