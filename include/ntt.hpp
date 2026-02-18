#pragma once
#include <vector>
#include "flints.hpp"
#include "vec64.hpp"
#include "montgomery.hpp"
#include "GPU/cuda_ntt.hpp"



class StandardNTTer {
private:
    size_t n_;
    size_t logn_;
    uint64_t q_;

    MontgomeryMultiplier mm_;
    vec64 roots_mont_;
    vec64 iroots_mont_;

    CudaBuffer roots_cuda_;
    CudaBuffer iroots_cuda_;
public:
    StandardNTTer(size_t n, uint64_t q, uint64_t qroot);
    ~StandardNTTer();

    void ntt(uint64_t* dst) const;
    void intt(uint64_t* dst) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;

    void ntt_batch_cuda(const CudaBuffer& dst, size_t batch_size) const;
    void intt_batch_cuda(const CudaBuffer& dst, size_t batch_size) const;

};

class TwistedNtterXY
{
private:
    int n_;
    uint64_t q_;
    vec64 zeta_pos_pows_mont_;
    vec64 zeta_neg_pows_mont_;
    StandardNTTer std_ntter_;
    MontgomeryMultiplier mm_;

    CudaBuffer zeta_pos_pows_mont_cuda_;
    CudaBuffer zeta_neg_pows_mont_cuda_;

public:
    TwistedNtterXY(int n, uint64_t q, uint64_t qroot);
    ~TwistedNtterXY();
    void ntt_mont(vec64& dst, const vec64& src) const;
    void intt_mont(vec64& dst, const vec64& src) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;

    void ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
    void intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
};




class TwistedNtterW
{
private:
    uint64_t p_;
    uint64_t q_;
    StandardNTTer std_ntter_;
    MontgomeryMultiplier mm_;

    // b_[i] = eta^{gamma^i}（蒙哥马利形式）
    vec64 b_;

    // binv_[i] = inv(etas[i])（蒙哥马利形式）
    vec64 binv_;

    CudaBuffer b_cuda_, b_inv_cuda_;
public:
    TwistedNtterW(int p , uint64_t q, uint64_t qroot);
    ~TwistedNtterW();
    void ntt_mont(vec64& dst, const vec64& src) const;
    void intt_mont(vec64& dst, const vec64& src) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;

    void ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
    void intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
};

class NTTerI
{
    MontgomeryMultiplier mm_;
    uint64_t I_mont_, I2_inv_mont_, inv2_mont_;
public:
    NTTerI(uint64_t q, uint64_t qroot);
    void ntt_batch(
        uint64_t* low, 
        uint64_t* hig, 
        size_t batch_size
    ) const;
    void intt_batch(
        uint64_t* low, 
        uint64_t* hig, 
        size_t batch_size
    ) const;
};