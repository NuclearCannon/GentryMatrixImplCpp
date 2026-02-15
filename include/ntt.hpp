#pragma once
#include <vector>
#include "flints.hpp"

#include "montgomery.hpp"




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

public:
    StandardNTTer(size_t n, uint64_t q, uint64_t qroot);
    ~StandardNTTer();

    void ntt_mont(uint64_t* dst, const uint64_t* src) const;
    void intt_mont(uint64_t* dst, const uint64_t* src) const;

};

