#include "complecx_matrix.hpp"


// let C = A @ B.conj().T
// 矩阵行优先编码，边长为n
// using complex = std::complex<double>;
void matmul_ABT_core(complex* C, const complex* A, const complex* B, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            complex sum = 0;
            for (size_t k = 0; k < n; ++k) {
                sum += A[i * n + k] * std::conj(B[j * n + k]);
            }
            C[i * n + j] = sum;
        }
    }
}


ComplexMatrixGroup ComplexMatrixGroup::matmul_ABT(const ComplexMatrixGroup& other) const
{
    ComplexMatrixGroup result(n_, p_);
    size_t nn = n_*n_;

    for(int w = 0; w<p_-1; w++)
    {
        matmul_ABT_core(
            result.data_.data()+w*nn,
            data_.data()+w*nn,
            other.data_.data()+w*nn,
            n_
        );
    }
    return result;
}