#include "FHE/encrypt_gp.hpp"
#include "FHE/key_switch_gp.hpp"
#include "FHE/circledast.hpp"
#include "flints.hpp"
#include <gperftools/profiler.h>
#include "timer.hpp"
#include "GPU/cuda_interf.hpp"
#include <iostream>

struct Params {
    uint64_t qo, qor;
    std::vector<std::pair<uint64_t, uint64_t>> qrp;
};

struct Params param_517 = {
    .qo = 576460752303421441, .qor = 19,
    .qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {576460752303421441, 19},
    }
};

struct Params param_257 = {
    .qo = 288230376153026561, .qor = 6,
    .qrp = {
        {1099513535489, 3},
        {1099516956673, 10},
        {1099519062017, 3},
        {288230376153026561, 6},
    }
};

void bench_circledast_param(
    int n, int p,
    bool cuda, int rpt
)
{
    const struct Params& param = (p==5 || p==17)?param_517:param_257;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = param.qrp;
    uint64_t qo = param.qo, qor = param.qor;
    vec64 mods, roots;
    for(int i=0; i<3;i++)
    {
        mods.push_back(qrp[i].first);
        roots.push_back(qrp[i].second);
    }
    // ================
    
    GentryPolyCtx ctx(n, p, qrp);
    GentryPoly sk = GentryPoly::sk(n, p, mods);
    auto ksk_pair = create_ksks_for_circledast_ct(sk, qo, ctx);
    
    double total_time = 0;
    // printf("[n=%d, p=%d, cuda=%d, repeat=%d]:start ", n, p, int(cuda), rpt);
    for(int rp=0; rp<rpt; rp++)
    {
        // printf("rp=%d ", rp);
        GentryPoly ua = GentryPoly::uniform(n, p, mods);
        GentryPoly ub = GentryPoly::uniform(n, p, mods);
        GentryPoly va = GentryPoly::uniform(n, p, mods);
        GentryPoly vb = GentryPoly::uniform(n, p, mods);
        // printf("随机生成完成\n");

        // 生成含有qo的sk1, sk2
        // 生成KSK
        if(!cuda)
        {
            HighResolutionTimer timer;
            timer.start();
            auto [ra, rb] = circledast_ct(ua, ub, va, vb, ksk_pair.first, ksk_pair.second, ctx);
            total_time += timer.stop();
        }
        else
        {
            GentryPoly uac = ua.to_cuda();
            GentryPoly ubc = ub.to_cuda();
            GentryPoly vac = va.to_cuda();
            GentryPoly vbc = vb.to_cuda();
            auto [ra_, rb_] = circledast_ct(uac, ubc, vac, vbc, ksk_pair.first, ksk_pair.second, ctx);
            CudaInterf::cuda_device_synchronize();
            HighResolutionTimer timer;
            timer.start();
            auto [ra, rb] = circledast_ct(uac, ubc, vac, vbc, ksk_pair.first, ksk_pair.second, ctx);
            CudaInterf::cuda_device_synchronize();
            total_time += timer.stop();
        }
    }

    printf("[n=%d, p=%d, cuda=%d, repeat=%d]:avgTime=%lf\n", n, p, int(cuda), rpt, total_time/rpt);
}

void bench_for_paper()
{
    std::vector<std::pair<int, int>> nplist = {
        // {16, 5},
        // {16, 17},
        // {16, 257},
        // {32, 5},
        // {32, 17},
        // {32, 257},
        // {64, 5},
        // {64, 17},
        // {64, 257},
        // {128, 5},
        // {128, 17},
        // {128, 257},
        // {256, 5},
        // {256, 17},
        {256, 257},
    };
    int rpt = 1;

    for(auto [n,p] : nplist)
    {
        bench_circledast_param(n, p, false, rpt);
        bench_circledast_param(n, p, true , rpt);
    }

}

void bench_circledast()
{
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };

    GentryPolyCtx ctx(n, p, qrp);
    GentryPoly u = GentryPoly::randint(n, p, mods, -1000, 1000);
    GentryPoly v = GentryPoly::randint(n, p, mods, -1000, 1000);
    GentryPoly sk = GentryPoly::sk(n, p, mods);

    auto [ua, ub] = encrypt_gp(u, sk, ctx);
    auto [va, vb] = encrypt_gp(v, sk, ctx);

    // 生成含有qo的sk1, sk2
    // 生成KSK
    auto ksk_pair = create_ksks_for_circledast_ct(sk, qo, ctx);
    auto u2 = decrypt_gp(ua, ub, sk, ctx);
    auto v2 = decrypt_gp(va, vb, sk, ctx);

    ProfilerStart("cd.prof");
    auto [ra, rb] = circledast_ct(ua, ub, va, vb, ksk_pair.first, ksk_pair.second, ctx);
    ProfilerStop();

    GentryPoly r = decrypt_gp(ra, rb, sk, ctx);
    // 直接计算
    GentryPoly w = GentryPoly::zeros_like(u);
    u2.iw_ntt(ctx);
    v2.iw_ntt(ctx);
    GentryPoly::circledast(w, u2, v2);  // 使用u2 @ v2作为对照组，因为我们忽略加密本身的噪声
    w.iw_intt(ctx);
    GentryPoly::sub(r,r,w);
    printf("bench_circledast: error_abs=%ld\n", r.abs());

}


void bench_circledast_cuda()
{
    printf("进入bench_circledast_cuda\n");
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };

    GentryPolyCtx ctx(n, p, qrp);
    printf("生成明文\n");
    GentryPoly u = GentryPoly::randint(n, p, mods, -1000, 1000);
    GentryPoly v = GentryPoly::randint(n, p, mods, -1000, 1000);
    printf("生成私钥\n");
    GentryPoly sk = GentryPoly::sk(n, p, mods);
    printf("加密\n");
    auto [ua, ub] = encrypt_gp(u, sk, ctx);
    auto [va, vb] = encrypt_gp(v, sk, ctx);
    printf("生成KSK\n");
    // 生成含有qo的sk1, sk2
    // 生成KSK
    auto ksk_pair = create_ksks_for_circledast_ct(sk, qo, ctx);
    printf("将明文转为cuda\n");
    GentryPoly uac = ua.to_cuda();
    GentryPoly ubc = ub.to_cuda();
    GentryPoly vac = va.to_cuda();
    GentryPoly vbc = vb.to_cuda();
    printf("开始CCMM\n");
    HighResolutionTimer timer;
    // cudaProfilerStart();
    CudaInterf::cuda_profiler_start();
    timer.start();

    auto [ra, rb] = circledast_ct(uac, ubc, vac, vbc, ksk_pair.first, ksk_pair.second, ctx);
    CudaInterf::cuda_device_synchronize();
    double time = timer.stop();
    CudaInterf::cuda_profiler_stop();
    std::cout << "time(ms):" << time << std::endl;

    GentryPoly r = decrypt_gp(ra.to_cpu(), rb.to_cpu(), sk, ctx);
    // 直接计算
    GentryPoly w = GentryPoly::zeros_like(u);
    auto u2 = decrypt_gp(ua, ub, sk, ctx);
    auto v2 = decrypt_gp(va, vb, sk, ctx);
    u2.iw_ntt(ctx);
    v2.iw_ntt(ctx);
    GentryPoly::circledast(w, u2, v2);  // 使用u2 @ v2作为对照组，因为我们忽略加密本身的噪声
    w.iw_intt(ctx);
    GentryPoly::sub(r,r,w);
    printf("bench_circledast_cuda: error_abs=%ld\n", r.abs());

}
