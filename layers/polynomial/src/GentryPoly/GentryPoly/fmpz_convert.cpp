#include "GentryPoly.hpp"

fmpz_vector GentryPoly::to_fmpz_vector() const
{
    auto data = this->to_coeffs();
    int n = this->n(), p = this->p();
    int len = 2*(p-1)*n*n;
    fmpz_vector result(len);
    icrt(result, data, moduli_);
    // 中心化
    // 首先，计算模数乘积
    fmpz_scalar Q = fmpz_scalar::from_si(1);
    for(auto i : moduli_)
    {
        fmpz_mul_ui(Q.raw(), Q.raw(), i);
    }
    return result.mod_centered(Q.raw());
}
GentryPoly GentryPoly::from_fmpz(size_t n, size_t p, const std::vector<uint64_t>& moduli, const fmpz_vector& data)
{
    auto crt_ed = crt(data, moduli);
    return GentryPoly::from_coeffs(n, p, moduli, crt_ed);
}