#include "FHE/encrypt64.hpp"
#include "FHE/key_switch_64.hpp"
/*

一组满足k*256*17*5+1的2^60级质数，
[1152921504606973441, 1152921504607191041, 1152921504607299841, 1152921504607691521, 1152921504608104961]
及其原根
[13, 3, 13, 11, 3]

*/

int test_ks64C()
{
    // 准备参数
    int n = 4;
    int p = 5;
    // 质数链
    vec64 mods = {1152921504606973441, 1152921504607191041, 1152921504607299841, 1152921504607691521};
    // 原根链
    vec64 roots = {13, 3, 13, 11};
    u64 qo = 1152921504608104961;
    u64 qor = 3;


    std::shared_ptr<U64CtxChain> cc = std::make_shared<U64CtxChain>(n, p, mods, roots);
    // 构造明文
    auto m = CRTArray::randint(cc, 1000, 2000);
    // 构造私钥
    auto sk = CRTArray::sk(cc);
    auto sk2 = CRTArray::sk(cc);
    // 加密（KS测试不考虑加密本身的噪声）
    auto [cta, ctb] = encrypt64_no_e(m, sk);
    // 构造KSK到sk
    KeySwitchKey64CRT ksk(sk, sk2, qo, qor);
    auto [cta2, ctb2] = ksk.key_switch_big_2(cta, ctb);
    // 解密
    auto res = decrypt64(cta2, ctb2, sk2);
    auto error = res.sub(m);
    // CRTArray error2 = CRTArray::dg(cc);
    fmpz_vector error_mpz = error.to_fmpz_vector_centered();
    // error_mpz.print();
    long error_abs = error_mpz.max_abs();
    printf("test_ks64C: error_abs=%ld\n", error_abs);
    return 1;

}