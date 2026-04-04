#pragma once
#include "flints.hpp"
#include <cstdint>
#include "montgomery.hpp"
#include "GPU/cuda_buffer.hpp"



void circledast_u64_cpu(uint64_t* C, const uint64_t* A, const uint64_t* B, size_t n, size_t p, uint64_t q);

void circledast_u64_gpu(const CudaBuffer& C, const CudaBuffer& A, const CudaBuffer& B, size_t n, size_t p, uint64_t q);