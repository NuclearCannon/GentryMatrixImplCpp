#include "FHE/encrypt64.hpp"
#include "FHE/key_switch_64.hpp"
#include <chrono>
#include <gperftools/profiler.h>

int test_ks64(bool test_base, bool test_crt)
{
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    u64 qo = 576460752303421441;
    u64 qor = 19;
    std::shared_ptr<U64CtxChain> cc = std::make_shared<U64CtxChain>(n, p, mods, roots);
    // 构造明文
    auto m = CRTArray::randint(cc, 1000, 2000);
    // 构造私钥
    auto sk = CRTArray::sk(cc);
    auto sk2 = CRTArray::sk(cc);
    // 加密（KS测试不考虑加密本身的噪声）
    auto [cta, ctb] = encrypt64_no_e(m, sk);
    // 构造KSK到sk
    if (test_base)
    {
        printf("生成kskb\n");
        KeySwitchKey64Base kskb(sk, sk2, 1UL<<50, 3, qo, qor);
        auto t1 = std::chrono::high_resolution_clock::now();
        ProfilerStart("bench64b.prof");
        auto [cta2, ctb2] = kskb.key_switch_big_2(cta, ctb);
        ProfilerStop();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration_ksk = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        printf("ks64B: %ld us\n", duration_ksk);
        // 解密
        auto res = decrypt64(cta2, ctb2, sk2);
        auto error = res.sub(m);
        fmpz_vector error_mpz = error.to_fmpz_vector_centered();
        long error_abs = error_mpz.max_abs();
        printf("ks64B: error_abs=%ld\n", error_abs);
    }
    if (test_crt)
    {
        printf("生成kskc\n");
        KeySwitchKey64CRT kskc(sk, sk2, qo, qor);
        auto t1 = std::chrono::high_resolution_clock::now();
        ProfilerStart("bench64c.prof");
        auto [cta2, ctb2] = kskc.key_switch_big_2(cta, ctb);
        ProfilerStop();
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration_ksk = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        printf("ks64C: %ld us\n", duration_ksk);
        // 解密
        auto res = decrypt64(cta2, ctb2, sk2);
        auto error = res.sub(m);
        fmpz_vector error_mpz = error.to_fmpz_vector_centered();
        long error_abs = error_mpz.max_abs();
        printf("ks64C: error_abs=%ld\n", error_abs);
    }
    return 1;

}