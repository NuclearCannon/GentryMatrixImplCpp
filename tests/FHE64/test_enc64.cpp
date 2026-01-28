#include "FHE/encrypt64.hpp"

/*

一组满足k*256*17*5+1的2^60级质数，
[1152921504606973441, 1152921504607191041, 1152921504607299841, 1152921504607691521, 1152921504608104961]
及其原根
[13, 3, 13, 11, 3]

*/

int test_enc64()
{
    // 准备参数
    int n = 16;
    int p = 5;
    // 质数链
    vec64 mods = {1152921504606973441, 1152921504607191041, 1152921504607299841, 1152921504607691521, 1152921504608104961};
    // 原根链
    vec64 roots = {13, 3, 13, 11, 3};
    std::shared_ptr<U64CtxChain> cc = std::make_shared<U64CtxChain>(n, p, mods, roots);
    // 构造明文
    auto m = CRTArray::randint(cc, 1000, 2000);
    // 构造私钥
    auto sk = CRTArray::sk(cc);
    // 加密
    auto [cta, ctb] = encrypt64(m, sk);
    // 解密
    auto res = decrypt64(cta, ctb, sk);
    auto error = res.sub(m);
    // CRTArray error2 = CRTArray::dg(cc);
    fmpz_vector error_mpz(2*(p-1)*n*n);
    error.to_fmpz_vector(error_mpz);
    
    fmpz_t Q;   // 大模数
    fmpz_init_set_ui(Q, mods[0]);
    fmpz_mul_ui(Q, Q, mods[1]);
    fmpz_mul_ui(Q, Q, mods[2]);
    fmpz_mul_ui(Q, Q, mods[3]);
    fmpz_mul_ui(Q, Q, mods[4]);
    fmpz_vector error_mpz_center = error_mpz.mod_centered(Q);
    error_mpz_center.print();
    long error_abs = error_mpz_center.max_abs();
    printf("test_enc64: error_abs=%ld\n", error_abs);

    
    fmpz_clear(Q);
    return 1;

}