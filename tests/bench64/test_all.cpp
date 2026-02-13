#include "FHE/encrypt64.hpp"
#include "FHE/key_switch_64.hpp"
#include "FHE/circledast.hpp"
#include <chrono>
#include <gperftools/profiler.h>

int test_all()
{
    printf("test_all\n");
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::shared_ptr<U64CtxChain> cc = std::make_shared<U64CtxChain>(n, p, mods, roots);
    printf("构造明文\n");
    auto u = CRTArray::randint(cc, 1000, 2000);
    auto v = CRTArray::randint(cc, 1000, 2000);
    
    // 构造私钥
    printf("构造私钥\n");
    auto sk = CRTArray::sk(cc);
    printf("构造密文\n");

    auto [ua, ub] = encrypt64_no_e(u, sk);
    auto [va, vb] = encrypt64_no_e(v, sk);
    printf("构造KSK\n");
    auto [ksk1, ksk2] = create_ksks_for_circledast_ct(sk, qo, qor);
    printf("执行运算\n");


    auto t1 = std::chrono::steady_clock::now();
    ProfilerStart("all.prof");
    auto [wa, wb] = circledast_ct(ua, ub, va, vb, ksk1, ksk2);
    ProfilerStop();
    auto t2 = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    printf("all: %ld us\n", duration);


    


    return 1;

}