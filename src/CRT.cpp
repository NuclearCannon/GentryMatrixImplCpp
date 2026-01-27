#include "CRT.hpp"
#include "uint64.hpp"


// CRT分解
std::vector<std::vector<u64>> crt(const fmpz_vector& src, const std::vector<u64>& mods)
{
    std::vector<std::vector<u64>> result;
    int len = src.len();
    fmpz_vector buf(len);
    for(u64 mod : mods)
    {
        fmpz_scalar mod_mpz(mod);
        _fmpz_vec_scalar_mod_fmpz(buf.raw(), src.raw(), len, mod_mpz.raw());
        result.push_back(buf.to_uint64());
    }
    return result;
}

// CRT合并
void icrt(fmpz_vector& dst, const std::vector<std::vector<u64>>& src, const std::vector<u64>& mods)
{
    int mod_len = mods.size();
    assert(mod_len >= 1);
    int data_len = src[0].size();

    // 初始化结果向量 x 为第一个模下的值
    fmpz_vector x(data_len);
    for (int j = 0; j < data_len; ++j) {
        fmpz_set_ui(x[j], src[0][j]);  // x_j ≡ src[0][j] (mod mods[0])
    }

    // 当前累积模数 M = mods[0]
    fmpz_scalar M;
    fmpz_set_ui(M.raw(), mods[0]);

    // 临时变量用于 fmpz_crt
    fmpz_t r1, r2, m1, m2, res;
    fmpz_init(r1); fmpz_init(r2);
    fmpz_init(m1); fmpz_init(m2);
    fmpz_init(res);

    // 逐个合并后续模数
    for (int i = 1; i < mod_len; ++i) {
        u64 mi = mods[i];
        fmpz_set_ui(m2, mi);          // m2 = mods[i]
        fmpz_set(m1, M.raw());        // m1 = 当前累积模数 M

        for (int j = 0; j < data_len; ++j) {
            fmpz_set(r1, x[j]);               // r1 = 当前解 x_j
            fmpz_set_ui(r2, src[i][j]);       // r2 = 新余数 src[i][j]

            // 调用 FLINT 的 CRT：res ≡ r1 (mod m1), res ≡ r2 (mod m2)
            fmpz_CRT(res, r1, m1, r2, m2, 0);  // 0 => 非负最小解

            fmpz_set(x[j], res);              // 更新 x[j]
        }

        // 更新累积模数 M = M * mods[i]
        fmpz_mul_ui(M.raw(), M.raw(), mi);
    }

    // 写入输出
    assert(dst.len() == data_len);
    _fmpz_vec_set(dst.raw(), x.raw(), data_len);

    // 清理临时变量
    fmpz_clear(r1); fmpz_clear(r2);
    fmpz_clear(m1); fmpz_clear(m2);
    fmpz_clear(res);
}
