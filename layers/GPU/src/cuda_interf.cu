#include "GPU/cuda_interf.hpp"
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include "GPU/cuda_check.hpp"

void CudaInterf::cuda_profiler_start()
{
    cudaProfilerStart();
}
void CudaInterf::cuda_profiler_stop()
{
    cudaProfilerStop();
}

void CudaInterf::cuda_device_synchronize()
{
    CUDA_CHECK(cudaDeviceSynchronize());
}