#pragma once
#include "GPU/cuda_buffer.hpp"
#include <cstdint>
#include "montgomery.hpp"
#include <vector>

class CudaMatmulTaskSet
{
private:
    int n_;
    uint64_t M_;
    std::vector<uint64_t*> Atasks, Btasks, Ctasks;

public:
    CudaMatmulTaskSet(int n, uint64_t q): n_(n), M_(q) {};
    void run();
    void append(const CudaBuffer& C,const CudaBuffer& A,const CudaBuffer& B);
    
};