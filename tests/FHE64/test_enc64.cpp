#include "FHE/encrypt64.hpp"



int test_enc64()
{
    // 准备参数
    int n = 16;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
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
    fmpz_vector error_mpz = error.to_fmpz_vector_centered();
    // error_mpz.print();
    long error_abs = error_mpz.max_abs();
    printf("test_enc64: error_abs=%ld\n", error_abs);
    return 1;

}