#pragma once
#include <vector>
#include "flints.hpp"

#include "montgomery.hpp"
#include "GPU/cuda_ntt.hpp"



class StandardNTTer {
private:
    size_t n_;
    size_t logn_;
    uint64_t q_;      // 模数q
    vec64 roots_;    // nroot的[0,n/2)次幂
    vec64 iroots_;    // nroot^{-1}的[0,n/2)次幂
    uint64_t ninv_;       // n_在模q意义下的乘法逆元

    // 蒙哥马利约简部分
    MontgomeryMultiplier mm;
    vec64 roots_mont_;
    vec64 iroots_mont_;
    uint64_t ninv_mont_;

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

class TwistedNtterXY64
{
private:
    int n_;
    uint64_t q_;
    vec64 zeta_pos_pows_;
    vec64 zeta_neg_pows_;
    vec64 zeta_pos_pows_mont_;
    vec64 zeta_neg_pows_mont_;
    StandardNTTer std_ntter;
    MontgomeryMultiplier mm;

    CudaBuffer zeta_pos_pows_mont_cuda_;
    CudaBuffer zeta_neg_pows_mont_cuda_;

public:
    TwistedNtterXY64(int n, uint64_t q, uint64_t qroot);
    ~TwistedNtterXY64();
    void ntt_mont(vec64& dst, const vec64& src) const;
    void intt_mont(vec64& dst, const vec64& src) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;

    void ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
    void intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
};




class TwistedNtterW64
{
private:
    uint64_t p_;
    uint64_t q_;
    StandardNTTer subntter;
    MontgomeryMultiplier mm;

    // b_[i] = eta^{gamma^i}（蒙哥马利形式）
    vec64 b_;

    // binv_[i] = inv(etas[i])（蒙哥马利形式）
    vec64 binv_;

    CudaBuffer b_cuda_, b_inv_cuda_;
    
    mutable vec64 buf1;
public:
    TwistedNtterW64(int p , uint64_t q, uint64_t qroot);
    ~TwistedNtterW64();
    void ntt_mont(vec64& dst, const vec64& src) const;
    void intt_mont(vec64& dst, const vec64& src) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;

    void ntt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
    void intt_batch_cuda(const CudaBuffer& a, size_t batch_size) const;
};