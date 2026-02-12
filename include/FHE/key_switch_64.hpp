#pragma once
#include "FHE/encrypt64.hpp"

// KSK: B基分解法
// 这种方法依赖于fmpz对密文中的a成分进行取模，基数B比较小，只需要一个额外模数qo

class KeySwitchKey64Base
{
private:
    std::shared_ptr<const U64CtxChain> cc_low_, cc_hig_;
    uint64_t B_,L_;
    std::vector<std::pair<CRTArray, CRTArray>> cts_;
public:
    KeySwitchKey64Base(const CRTArray& sk_from, const CRTArray& sk_to, uint64_t B, uint64_t L, uint64_t qo, uint64_t qor);
    std::pair<CRTArray, CRTArray> key_switch_big_1(const CRTArray& a) const;
    std::pair<CRTArray, CRTArray> key_switch_big_2(const CRTArray& a, const CRTArray& b) const;
};

// KSK: 从CRT表示中提取分量
// 这种方式计算量极小，但是可能不太好控制噪声
// 暂且地，也只加入一个额外模数
class KeySwitchKey64CRT
{
private:
    std::shared_ptr<const U64CtxChain> cc_low_, cc_hig_;
    std::vector<std::pair<CRTArray, CRTArray>> cts_;
    KeySwitchKey64CRT(std::shared_ptr<const U64CtxChain> cc_low,std::shared_ptr<const U64CtxChain> cc_hig, std::vector<std::pair<CRTArray, CRTArray>> cts);

    // mutable
    mutable CRTArray buf_a_sum_;
    mutable CRTArray buf_b_sum_;
    mutable CRTArray buf_mul_result_;
    mutable CRTArray buf_ntted_;
    mutable CRTArray buf_a_split_i_;
    mutable CRTArray buf_a_res_;
    mutable CRTArray buf_b_res_;

public:
    

    static KeySwitchKey64CRT ksk_gen(const CRTArray& sk_from, const CRTArray& sk_to, uint64_t qo, uint64_t qor);


    std::pair<CRTArray, CRTArray> key_switch_big_1(const CRTArray &a) const;
    std::pair<CRTArray, CRTArray> key_switch_big_2(const CRTArray& a, const CRTArray& b) const;
};

std::pair<KeySwitchKey64CRT, KeySwitchKey64CRT> create_ksks_for_circledast_ct(const CRTArray& sk, uint64_t qo, uint64_t qor);