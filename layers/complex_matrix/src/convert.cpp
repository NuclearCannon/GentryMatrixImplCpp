#include "complex_matrix.hpp"
#include <cmath>

#include "flints.hpp"
#include <mpfr.h>
#include <flint/fmpz.h>

long double fmpz_get_ld_via_str(const fmpz_t x) {
    long double res;
    char *str;

    // 将 fmpz 转换为 10 进制字符串
    // NULL 表示让 FLINT 自己分配内存
    str = fmpz_get_str(NULL, 10, x);

    // 使用标准库函数将字符串转为 long double
    res = strtold(str, NULL);

    // 释放 FLINT 分配的内存 (flint_free 是 FLINT 内部的释放函数)
    flint_free(str);

    return res;
}

// 取一个 fmpz_t 整数对应的实数表示
long double flint_to_double(const fmpz_t src)
{
    return fmpz_get_ld_via_str(src);
}

// 取一个实数最接近的整数（四舍五入）
void double_to_flint(fmpz_t dst, long double src)
{
    // 处理 NaN 和无穷大？
    assert(std::isfinite(src));
    // FLINT 的 fmpz_set_d 内部使用 round()，即四舍五入到最近整数
    fmpz_set_d(dst, src);
}

fmpz_vector ComplexMatrixGroup::to_fmpz_vector(long double delta) const
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

ComplexMatrixGroup ComplexMatrixGroup::from_fmpz_vector(const fmpz_vector& src, long double delta, size_t n, size_t p)
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