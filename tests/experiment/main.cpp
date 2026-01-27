#include "ntt.hpp"
#include "CRT.hpp"
/*

本论文做这样的试验
进行同样一个多项式乘法
方法1.ntt_standard
方法2.先CRT，然后做NTT64，然后乘，然后iNTT64，然后iCRT

*/



const int n = 256;
// 1393796574908163946345982392040522594125057 1194085069146539581683975312461909755232262 112538967229589839993328056454558914054113

fmpz_scalar Q("1393796574908163946345982392040522594125057");
fmpz_scalar Qroot("1194085069146539581683975312461909755232262");
fmpz_scalar root_inv("112538967229589839993328056454558914054113");

fmpz_mod_ctx_t Qctx;
// 选取几个约为2^60的质数
// 质数: 1152921504606847009       原根: 13
// 质数: 1152921504606847067       原根: 6
// 质数: 1152921504606847081       原根: 11
// 质数: 1152921504606847123       原根: 2
// 质数: 1152921504606847127       原根: 5

u_int64_t mod_pow(u_int64_t base, u_int64_t e, u_int64_t mod)
{
    u_int64_t result = 1;
    base = base % mod;
    while (e > 0) {
        if (e & 1) {
            result = mod_mul(result, base, mod);
        }
        base = mod_mul(base, base, mod);
        e >>= 1;
    }
    return result;
}

u_int64_t mod_inv(u_int64_t x, u_int64_t mod)
{
    return mod_pow(x, mod-1, mod);
}

/*
质数                   256阶本原单位根      单位根的逆元
1152921504606876929 302590808370805342 203576432732386606
1152921504606877697 1039092709060134491 550634830399138777
1152921504606889217 272532665442783039 806171715214065977
1152921504606894593 568358301744714721 431433276899556624
1152921504606899969 898653445164898652 912894953107837243
*/


std::vector<u_int64_t> mods = {1152921504606876929, 1152921504606877697, 1152921504606889217, 1152921504606894593, 1152921504606899969};
// 原根
std::vector<u_int64_t> roots = {302590808370805342, 1039092709060134491, 272532665442783039, 568358301744714721, 898653445164898652};
// N单位根
std::vector<u_int64_t> roots_inv = {203576432732386606, 550634830399138777, 806171715214065977, 431433276899556624, 912894953107837243};




// fmpz_vector 



int main()
{
    // 准备阶段
    fmpz_mod_ctx_init(Qctx, Q.raw());
    // 检查单位根选择正确
    for(int i=0;i<5;i++)
    {
        assert(mod_mul(roots[i], roots_inv[i], mods[i]) == 1);
    }
    {
        fmpz_scalar tmp;
        fmpz_mod_mul(tmp.raw(), Qroot.raw(), root_inv.raw(), Qctx);
        assert(fmpz_is_one(tmp.raw()));
    }

    // 生成数据
    fmpz_vector a = fmpz_vector::uniform(n, Q.raw());
    fmpz_vector b = fmpz_vector::uniform(n, Q.raw());


    // 方法1
    // 对a b NTT
    fmpz_vector a_ntt(n),b_ntt(n), c_ntt(n), c(n);
    ntt_standard_flint(a, a_ntt, Qroot.raw(), n, Qctx);
    ntt_standard_flint(b, b_ntt, Qroot.raw(), n, Qctx);
    // 对a_ntt, b_ntt相乘
    _fmpz_mod_vec_mul(c_ntt.raw(), a_ntt.raw(), b_ntt.raw(), n, Qctx);
    // 逆NTT（不含除以N）
    ntt_standard_flint(c_ntt, c, root_inv.raw(), n, Qctx);


    // 方法2
    // 首先对原数据进行CRT
    auto a_crt = crt(a, mods);
    auto b_crt = crt(a, mods);
    // 对于crt数据，每个进行NTT
    std::vector<std::vector<u_int64_t>> a_crt_ntt, b_crt_ntt, c_crt_ntt, c_crt;
    for(int i=0;i<5;i++)
    {
        std::vector<u_int64_t> tmp(n);
        ntt_standard_64(a_crt[i].data(), tmp.data(), roots[i], n, mods[i]);
        a_crt_ntt.push_back(std::move(tmp));
    }
    for(int i=0;i<5;i++)
    {
        std::vector<u_int64_t> tmp(n);
        ntt_standard_64(b_crt[i].data(), tmp.data(), roots[i], n, mods[i]);
        b_crt_ntt.push_back(std::move(tmp));
    }
    // 相乘
    for(int i=0;i<5;i++)
    {
        u_int64_t m = mods[i];

        std::vector<u_int64_t> tmp(n);
        std::vector<u_int64_t> &s1 = a_crt_ntt[i], &s2 = b_crt_ntt[i];

        for(int j=0;j<n;j++)tmp[j] = mod_mul(s1[j], s2[j], m);
        c_crt_ntt.push_back(std::move(tmp));
    }
    // 逆NTT
    for(int i=0;i<5;i++)
    {
        std::vector<u_int64_t> tmp(n);
        ntt_standard_64(c_crt_ntt[i].data(), tmp.data(), roots_inv[i], n, mods[i]);
        c_crt.push_back(std::move(tmp));
    }
    // CRT还原
    fmpz_vector c2(n);
    icrt(c2, c_crt, mods);
    // 手动取模Q
    _fmpz_mod_vec_set_fmpz_vec(c2.raw(), c2.raw(), n, Qctx);
    // 比较结果: c vs c2
    // c.print();
    // c2.print();
    int error = 0;
    for(int i=0;i<n;i++)
    {
        if(fmpz_cmp(c[i], c2[i]))
        {
            printf("不一致：%d\n", i);
            error = 1;
        }
    }
    if (error == 0)
    {
        printf("完全一致！\n");
    }


    // 收尾阶段
    fmpz_mod_ctx_clear(Qctx);
}