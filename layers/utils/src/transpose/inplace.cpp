#include "math_utils.hpp"
#include <algorithm>
#include <cstring>


// #define NAIVE_TRANSPOSE_INPLACE

#ifdef NAIVE_TRANSPOSE_INPLACE

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


#else

constexpr int BLOCK = 16; // 可尝试 32/64/128


// 原地转置一个 B×B 的小方阵（B 较小，用朴素或分治）
static inline void transpose_block_inplace(uint64_t* mat, int B, int stride) {
    // stride 通常是 n（主矩阵列数）
    for (int i = 0; i < B; ++i) {
        for (int j = i + 1; j < B; ++j) {
            std::swap(mat[i * stride + j], mat[j * stride + i]);
        }
    }
}

// if B == BLOCK
static inline void transpose_block_inplace_block(uint64_t* mat, int stride) {
    // stride 通常是 n（主矩阵列数）
    for (int i = 0; i < BLOCK; ++i) {
        for (int j = i + 1; j < BLOCK; ++j) {
            std::swap(mat[i * stride + j], mat[j * stride + i]);
        }
    }
}


void transpose_inplace_aligned(uint64_t* A, int n) {
    // 分块处理
    assert((n & (n - 1)) == 0);
    assert(n % BLOCK == 0);
    uint64_t temp[BLOCK * BLOCK];

    for (int ii = 0; ii < n; ii += BLOCK) {
        // 先处理对角线块：原地转置 A[ii:imax, ii:imax]
        transpose_block_inplace_block(&A[ii * n + ii], n);

        // 再处理非对角线块：交换 (ii, jj) 与 (jj, ii)
        for (int jj = ii + BLOCK; jj < n; jj += BLOCK) {
            // 复制 A[ii:imax, jj:jmax] 到 temp
            for (int i = 0; i < BLOCK; ++i) {
                memcpy(&temp[i * BLOCK],
                       &A[(ii + i) * n + jj],
                       BLOCK * sizeof(uint64_t));
            }
            // 将 A[jj:jmax, ii:imax] 复制到 A[ii:imax, jj:jmax]
            for (int i = 0; i < BLOCK; ++i) {
                for (int j = 0; j < BLOCK; ++j) {
                    A[(ii + i) * n + jj + j] = A[(jj + j) * n + ii + i];
                }
            }
            // 将 temp 复制到 A[jj:jmax, ii:imax]（要转置）
            for(int i=0; i<BLOCK; i++)
            {
                for(int j=0; j<BLOCK; j++)
                {
                    A[(jj + j)*n + ii + i] = temp[i*BLOCK + j];
                }
            }
            
        }
    }
}

void transpose_inplace(uint64_t* A, int n) {
    // 分块处理
    assert((n & (n - 1)) == 0);
    if (n % BLOCK == 0)return transpose_inplace_aligned(A, n);
    uint64_t temp[BLOCK * BLOCK];

    for (int ii = 0; ii < n; ii += BLOCK) {
        int imax = std::min(ii + BLOCK, n);
        int Bi = imax - ii;

        // 先处理对角线块：原地转置 A[ii:imax, ii:imax]
        transpose_block_inplace(&A[ii * n + ii], Bi, n);

        // 再处理非对角线块：交换 (ii, jj) 与 (jj, ii)
        for (int jj = ii + BLOCK; jj < n; jj += BLOCK) {
            int jmax = std::min(jj + BLOCK, n);
            int Bj = jmax - jj;

            // 复制 A[ii:imax, jj:jmax] 到 temp
            for (int i = 0; i < Bi; ++i) {
                memcpy(&temp[i * Bj],
                       &A[(ii + i) * n + jj],
                       Bj * sizeof(uint64_t));
            }

            // 将 A[jj:jmax, ii:imax] 复制到 A[ii:imax, jj:jmax]
            for (int i = 0; i < Bi; ++i) {
                for (int j = 0; j < Bj; ++j) {
                    A[(ii + i) * n + jj + j] = A[(jj + j) * n + ii + i];
                }
            }

            // 将 temp 复制到 A[jj:jmax, ii:imax]（要转置）
            for(int i=0; i<Bi; i++)
            {
                for(int j=0; j<Bj; j++)
                {
                    A[(jj + j)*n + ii + i] = temp[i*Bj + j];
                }
            }
            
        }
    }
}

#endif