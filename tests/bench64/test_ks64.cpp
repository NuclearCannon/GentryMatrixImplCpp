#include "FHE/encrypt64.hpp"
#include "FHE/key_switch_64.hpp"
#include <chrono>

int test_ks64()
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
    int logB = 20;
    printf("生成kskb\n");
    KeySwitchKey64Base kskb(sk, sk2, 1UL<<logB, ((61*4)/logB), qo, qor);
    printf("生成kskc\n");
    KeySwitchKey64CRT kskc(sk, sk2, qo, qor);
    printf("开始kskb\n");
    auto t1 = std::chrono::high_resolution_clock::now();
    auto [cta2, ctb2] = kskb.key_switch_big_2(cta, ctb);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration_ksk = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    printf("ks64B: %ld us\n", duration_ksk);
    printf("开始kskc\n");
    t1 = std::chrono::high_resolution_clock::now();
    auto [cta3, ctb3] = kskc.key_switch_big_2(cta, ctb);
    t2 = std::chrono::high_resolution_clock::now();
    duration_ksk = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    printf("ks64C: %ld us\n", duration_ksk);
    return 1;

}