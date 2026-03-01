#include "matmul.hpp"
#include "vec64.hpp"
#include "montgomery.hpp"
#include "GPU/cuda_matmul.hpp"




MatmulContext::MatmulContext(int n, const fmpz_t q):
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

void MatmulContext::matmul_transpose(fmpz* C, const fmpz* A, const fmpz* B) const
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

void MatmulContext::matmul_transpose_u64(uint64_t* C, const uint64_t* A, const uint64_t* B) const
{
    // 填充
	for(int i=0;i<n_;i++)for(int j=0;j<n_;j++)
	{
		
		fmpz_set_ui(fmpz_mod_mat_entry(A_mat, i, j), *(A + i*n_ + j));
        fmpz_set_ui(fmpz_mod_mat_entry(B_mat, j, i), *(B + i*n_ + j));   // 在这转置！
	}
	fmpz_mod_mat_mul(C_mat, A_mat, B_mat);
    // 写回
    for(int i=0;i<n_;i++)for(int j=0;j<n_;j++)
	{
        *(C + i*n_ + j) = fmpz_get_ui(fmpz_mod_mat_entry(C_mat, i, j));
	}
}


fmpz_vector MatmulContext::circledast_fmpz(const fmpz_vector& A, const fmpz_vector& B, size_t n, size_t p) 
{
    assert(n == this->n_);
    size_t nn = n*n;
    size_t pnn = (p-1)*nn;

    fmpz_vector result(2*pnn);

    vec64 gpp = get_powers(3,p-1,p), gpp_backward(p);
    for(int i=0; i<p-1; i++)gpp_backward[gpp[i]] = i;

    // [w]位置是多项式在(eta^{3^w})上的取值
    // 它需要和(eta^{p-3^w})的那一份相乘

    for(size_t w=0; w<p-1; w++)
    {
        size_t w2 = gpp_backward[p-gpp[w]];
        size_t base_a = w*nn;
        size_t base_b = w2*nn;
        size_t base_ap = base_a;
        size_t base_an = base_a + pnn;
        size_t base_bp = base_b + pnn;
        size_t base_bn = base_b;
        // 做矩阵乘法
        // rp = ap @ bp.T
        // rn = an @ bn.T
        matmul_transpose(
            result[base_ap], A[base_ap], B[base_bp]
        );
        matmul_transpose(
            result[base_an], A[base_an], B[base_bn]
        );
    }
    return result;
}

void MatmulContext::circledast_u64(uint64_t* dst, const uint64_t* A, const uint64_t* B, size_t n, size_t p)
{
    assert(n == this->n_);
    size_t nn = n*n;
    size_t pnn = (p-1)*nn;

    std::vector<size_t> gpp = get_powers(3,p-1,p), gpp_backward(p);
    for(int i=0; i<p-1; i++)gpp_backward[gpp[i]] = i;

    // [w]位置是多项式在(eta^{3^w})上的取值
    // 它需要和(eta^{p-3^w})的那一份相乘

    for(size_t w=0; w<p-1; w++)
    {
        size_t w2 = gpp_backward[p-gpp[w]];
        size_t base_a = w*nn;
        size_t base_b = w2*nn;
        size_t base_ap = base_a;
        size_t base_an = base_a + pnn;
        size_t base_bp = base_b + pnn;
        size_t base_bn = base_b;
        // 做矩阵乘法
        // rp = ap @ bp.T
        // rn = an @ bn.T
        matmul_transpose_u64(
            dst+base_ap, A + base_ap, B + base_bp
        );
        matmul_transpose_u64(
            dst+base_an, A+base_an, B+base_bn
        );
    }
}

void circledast_u64_gpu(uint64_t* dst, const uint64_t* A, const uint64_t* B, size_t n, size_t p, const MontgomeryMultiplier& mm)
{
    size_t nn = n*n;
    size_t pnn = (p-1)*nn;
    std::vector<size_t> gpp = get_powers(3,p-1,p), gpp_backward(p);
    for(int i=0; i<p-1; i++)gpp_backward[gpp[i]] = i;

    // [w]位置是多项式在(eta^{3^w})上的取值
    // 它需要和(eta^{p-3^w})的那一份相乘

    CudaBuffer Abuf(nn*sizeof(uint64_t)), Bbuf(nn*sizeof(uint64_t)), Cbuf(nn*sizeof(uint64_t));
    for(size_t w=0; w<p-1; w++)
    {
        size_t w2 = gpp_backward[p-gpp[w]];
        size_t base_a = w*nn;
        size_t base_b = w2*nn;
        size_t base_ap = base_a;
        size_t base_an = base_a + pnn;
        size_t base_bp = base_b + pnn;
        size_t base_bn = base_b;
        // 做矩阵乘法
        // rp = ap @ bp.T
        Abuf.copy_from_host(A+base_ap);
        Bbuf.copy_from_host(B+base_bp);
        matmul_gpu(Cbuf, Abuf, Bbuf, n, mm);
        Cbuf.copy_to_host(dst+base_ap);
        // rn = an @ bn.T
        Abuf.copy_from_host(A+base_an);
        Bbuf.copy_from_host(B+base_bn);
        matmul_gpu(Cbuf, Abuf, Bbuf, n, mm);
        Cbuf.copy_to_host(dst+base_an);
        
    }
}

void circledast_u64_gpu2(const CudaBuffer& C, const CudaBuffer& A, const CudaBuffer& B, size_t n, size_t p, const MontgomeryMultiplier& mm)
{
    size_t nn = n*n;
    size_t pnn = (p-1)*nn;
    std::vector<size_t> gpp = get_powers(3,p-1,p), gpp_backward(p);
    for(int i=0; i<p-1; i++)gpp_backward[gpp[i]] = i;

    // [w]位置是多项式在(eta^{3^w})上的取值
    // 它需要和(eta^{p-3^w})的那一份相乘

    for(size_t w=0; w<p-1; w++)
    {
        size_t w2 = gpp_backward[p-gpp[w]];
        size_t base_a = w*nn;
        size_t base_b = w2*nn;
        size_t base_ap = base_a;
        size_t base_an = base_a + pnn;
        size_t base_bp = base_b + pnn;
        size_t base_bn = base_b;

        constexpr size_t s64 = sizeof(uint64_t);
        // 做矩阵乘法
        // rp = ap @ bp.T
        matmul_gpu(
            C.slice(base_ap*s64, (base_ap+nn)*s64), 
            A.slice(base_ap*s64, (base_ap+nn)*s64), 
            B.slice(base_bp*s64, (base_bp+nn)*s64), 
            n, mm
        );

        matmul_gpu(
            C.slice(base_an*s64, (base_an+nn)*s64), 
            A.slice(base_an*s64, (base_an+nn)*s64), 
            B.slice(base_bn*s64, (base_bn+nn)*s64), 
            n, mm
        );

    }
}