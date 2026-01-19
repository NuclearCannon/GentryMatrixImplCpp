#include "ziq_array.hpp"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_mat.h"



MatmulContext::MatmulContext(int n, fmpz_t q):
    n_(n)
{
    fmpz_init_set(q_, q);
    fmpz_mod_mat_init(A_mat, n_, n_, q_);
    fmpz_mod_mat_init(B_mat, n_, n_, q_);
    fmpz_mod_mat_init(C_mat, n_, n_, q_);

}
MatmulContext::~MatmulContext()
{
    fmpz_mod_mat_clear(A_mat);
    fmpz_mod_mat_clear(B_mat);
    fmpz_mod_mat_clear(C_mat);
    fmpz_clear(q_);
}

void MatmulContext::matmul_transpose(fmpz* C, const fmpz* A, const fmpz* B)
{
    // 填充
	for(int i=0;i<n_;i++)for(int j=0;j<n_;j++)
	{
		
		fmpz_set(fmpz_mod_mat_entry(A_mat, i, j), (A + i*n_ + j));
        fmpz_set(fmpz_mod_mat_entry(B_mat, j, i), (B + i*n_ + j));   // 在这转置！
	}
	fmpz_mod_mat_mul(C_mat, A_mat, B_mat);
    // 写回
    for(int i=0;i<n_;i++)for(int j=0;j<n_;j++)
	{
		fmpz_set((C + i*n_ + j), fmpz_mod_mat_entry(C_mat, i, j));
	}
}


fmpz_vector ZiqArrayContext::circledast(const fmpz_vector& A, const fmpz_vector& B) const
{
    fmpz_vector result(size_);
    for(size_t w=0; w<p_-1; w++)
    {
        size_t base_a = w*nn_;
        size_t base_b = (p_-2-w)*nn_;
        size_t base_ap = base_a;
        size_t base_an = base_a + pnn_;
        size_t base_bp = base_b + pnn_;
        size_t base_bn = base_b;
        // 做矩阵乘法
        // rp = ap @ bp.T
        // rn = an @ bn.T
        mm_ctx->matmul_transpose(
            result[base_ap], A[base_ap], B[base_bp]
        );
        mm_ctx->matmul_transpose(
            result[base_an], A[base_an], B[base_bn]
        );
    }
    return result;
}