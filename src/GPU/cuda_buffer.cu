
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "GPU/cuda_buffer.hpp"
#include "GPU/cuda_check.hpp"

CudaBuffer::CudaBuffer(
    size_t size_bytes
)
{
    CUDA_CHECK(cudaMalloc(&ptr_, size_bytes));
    if(ptr_ == nullptr)
    {
        std::cerr << "CudaBuffer" <<std::endl;
    }
    len_ = size_bytes;
}

CudaBuffer::~CudaBuffer()
{
    if(ptr_)
    {
        cudaFree(ptr_);
        ptr_ = 0;
    }
}

CudaBuffer::CudaBuffer(CudaBuffer&& other) noexcept
{
    ptr_ = other.ptr_;
    other.ptr_ = 0;
    len_ = other.len_;
    other.len_ = 0;
}
CudaBuffer& CudaBuffer::operator=(CudaBuffer&& other) noexcept
{
    if(ptr_)
    {
        cudaFree(ptr_);
        ptr_ = 0;
    }
    ptr_ = other.ptr_;
    other.ptr_ = 0;
    len_ = other.len_;
    other.len_ = 0;
    return *this;
}

void CudaBuffer::copy_to_host(void* host_ptr) const
{
    CUDA_CHECK(cudaMemcpy(host_ptr, ptr_, len_, cudaMemcpyDeviceToHost));

}
void CudaBuffer::copy_from_host(const void* host_ptr)
{
    CUDA_CHECK(cudaMemcpy(ptr_, host_ptr, len_, cudaMemcpyHostToDevice));
}
