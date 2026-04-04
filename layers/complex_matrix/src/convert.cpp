#include "complecx_matrix.hpp"
#include <cmath>

// 取一个 fmpz_t 整数对应的实数表示
double flint_to_double(const fmpz_t src)
{
    return fmpz_get_d(src);
}

// 取一个实数最接近的整数（四舍五入）
void double_to_flint(fmpz_t dst, double src)
{
    // 处理 NaN 和无穷大？
    assert(std::isfinite(src));
    // FLINT 的 fmpz_set_d 内部使用 round()，即四舍五入到最近整数
    fmpz_set_d(dst, src);
}

fmpz_vector ComplexMatrixGroup::to_fmpz_vector(double delta) const
{
    ssize_t pnn = (p_-1)*n_*n_;
    fmpz_vector result(2*pnn);
    for(int i=0; i<pnn; i++)
    {
        double_to_flint(result[i],     data_[i].real() * delta);
        double_to_flint(result[i+pnn], data_[i].imag() * delta);
    }
    return result;
}

ComplexMatrixGroup ComplexMatrixGroup::from_fmpz_vector(const fmpz_vector& src, double delta, size_t n, size_t p)
{
    ssize_t pnn = (p-1)*n*n;
    assert(src.len() == 2*pnn);
    std::vector<complex> result(pnn);

    for(int i=0; i<pnn; i++)
    {
        result[i].real(flint_to_double(src[i])/delta);
        result[i].imag(flint_to_double(src[i+pnn])/delta);
    }
    return ComplexMatrixGroup(n, p, result);
}