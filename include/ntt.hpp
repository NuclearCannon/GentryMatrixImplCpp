#pragma once
#include <vector>
#include <unordered_map>
#include "flints.hpp"
#include "uint64.hpp"

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
    u64 qroot_;  // q的一个原根
    u64 nroot_;   // = qroot ^ ((q-1)/n)

    vec64 roots_;    // nroot的[0,n/2)次幂
    vec64 iroots_;    // nroot^{-1}的[0,n/2)次幂
    u64 ninv_;       // n_在模q意义下的乘法逆元

    
    void _ntt_standard_inner(u64* dst, const u64* src, const u64* roots) const;
public:
    StandardNTTer(size_t n, u64 q, u64 qroot);
    ~StandardNTTer();
    void ntt(u64* dst, const u64* src) const;
    void intt(u64* dst, const u64* src) const;

};


// Rader NTT
class RaderNTTer64 {
private:
    u64 p_, g_, eta_, q_, pinv_;
    vec64 gpp, gnp;
    vec64 b1ntt, b2ntt;
    StandardNTTer subntter;
    void _rader_inner(u64* dst, const u64* src, const vec64& bntt) const;

public:
    RaderNTTer64(u64 p, u64 q, u64 qroot);
    ~RaderNTTer64();
    
    void rader(u64* dst, const u64* src) const;
    void irader(u64* dst, const u64* src) const;

};