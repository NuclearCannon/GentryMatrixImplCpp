#include "FHE/encrypt_gp.hpp"
#include "FHE/key_switch_gp.hpp"
#include "FHE/circledast.hpp"
#include "CRT.hpp"
#include <gperftools/profiler.h>
#include "timer.hpp"
#include <cuda_profiler_api.h>
#include <cuda_runtime_api.h>
#include "GPU/cuda_check.hpp"
#include <iostream>



void bench_matmul_cuda()
{
    printf("bench_matmul_cuda\n");
    // 准备参数
    int n = 256;
    int p = 17;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};

    GentryPoly a = GentryPoly::uniform(n, p, mods).to_cuda();
    GentryPoly b = GentryPoly::uniform(n, p, mods).to_cuda();
    GentryPoly c = GentryPoly::uniform(n, p, mods).to_cuda();
    GentryPoly d = GentryPoly::uniform(n, p, mods).to_cuda();
    GentryPoly ac = GentryPoly::zeros(n, p, mods, GPDevice::CUDA);
    GentryPoly ad = GentryPoly::zeros(n, p, mods, GPDevice::CUDA);
    GentryPoly bc = GentryPoly::zeros(n, p, mods, GPDevice::CUDA);
    GentryPoly bd = GentryPoly::zeros(n, p, mods, GPDevice::CUDA);
    

    HighResolutionTimer timer;
    cudaProfilerStart();
    timer.start();
    GentryPoly::circledast(ac, a, c);
    GentryPoly::circledast(ad, a, d);
    GentryPoly::circledast(bc, b, c);
    GentryPoly::circledast(bd, b, d);
    CUDA_CHECK(cudaDeviceSynchronize());
    double time = timer.stop();
    cudaProfilerStop();
    std::cout << "time(ms):" << time << std::endl;

}
