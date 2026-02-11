#pragma once
#include <vector>
#include <unordered_map>
#include "flints.hpp"
#include "uint64.hpp"
#include "montgomery.hpp"

const std::vector<size_t>& get_bit_reverse_table(size_t n);
const std::vector<size_t>& get_bit_reverse_table_by_logn(size_t log2n);

constexpr int log2(int x)
{
    assert((x & (x - 1)) == 0 && x >= 1);
    int i=-1;
    while(x)
    {
        x>>=1;
        i++;
    }    
    return i;
}


class StandardNTTer {
private:
    size_t n_;
    size_t logn_;
    u64 q_;      // 模数q
    vec64 roots_;    // nroot的[0,n/2)次幂
    vec64 iroots_;    // nroot^{-1}的[0,n/2)次幂
    u64 ninv_;       // n_在模q意义下的乘法逆元

    // 蒙哥马利约简部分
    MontgomeryMultiplier mm;
    vec64 roots_mont_;
    vec64 iroots_mont_;
    u64 ninv_mont_;

    void _ntt_standard_inner_mont(u64* dst, const u64* src, const u64* roots) const;
public:
    StandardNTTer(size_t n, u64 q, u64 qroot);
    ~StandardNTTer();

    void ntt_mont(u64* dst, const u64* src) const;
    void intt_mont(u64* dst, const u64* src) const;

};


// Rader NTT
class RaderNTTer64 {
private:
    u64 p_, eta_, q_, pinv_;
    vec64 gpp, gnp;
    vec64 b1ntt, b2ntt;
    StandardNTTer subntter;
    

    // 蒙哥马利约简部分
    MontgomeryMultiplier mm;
    vec64 b1ntt_mont, b2ntt_mont;

    mutable vec64 buf_a_, buf_c_, buf_a_ntt_, buf_c_ntt_;
    void _rader_inner_mont(u64* dst, const u64* src, const vec64& bntt) const;

public:
    RaderNTTer64(u64 p, u64 q, u64 qroot);
    ~RaderNTTer64();
    
    void rader_mont(u64* dst, const u64* src) const;
    void irader_mont(u64* dst, const u64* src) const;

};