#pragma once
#include "FHE/encrypt64.hpp"

// KSK: B基分解法
// 这种方法依赖于fmpz对密文中的a成分进行取模，基数B比较小，只需要一个额外模数qo

class KeySwitchKey64Base
{
private:
    std::shared_ptr<const U64CtxChain> cc_low_, cc_hig_;
    u64 B_,L_;
    std::vector<std::pair<CRTArray, CRTArray>> cts_;
public:
    KeySwitchKey64Base(const CRTArray& sk_from, const CRTArray& sk_to, u64 B, u64 L, u64 qo, u64 qor);
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
public:
    KeySwitchKey64CRT(const CRTArray& sk_from, const CRTArray& sk_to, u64 qo, u64 qor);
    std::pair<CRTArray, CRTArray> key_switch_big_1(const CRTArray& a) const;
    std::pair<CRTArray, CRTArray> key_switch_big_2(const CRTArray& a, const CRTArray& b) const;
};