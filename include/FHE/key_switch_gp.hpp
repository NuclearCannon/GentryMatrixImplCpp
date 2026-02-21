#pragma once
#include "FHE/encrypt_gp.hpp"


// KSK: 从CRT表示中提取分量
// 这种方式计算量极小，但是可能不太好控制噪声
// 暂且地，也只加入一个额外模数
class KeySwitchKeyGP
{
private:
    std::vector<std::pair<GentryPoly, GentryPoly>> cts_;
    uint64_t qo_;
    KeySwitchKeyGP(std::vector<std::pair<GentryPoly, GentryPoly>> cts, uint64_t qo);
public:
    static KeySwitchKeyGP ksk_gen(const GentryPoly& sk_from, const GentryPoly& sk_to, uint64_t qo);
    std::pair<GentryPoly, GentryPoly> key_switch_big_1(const GentryPoly &a) const;
    std::pair<GentryPoly, GentryPoly> key_switch_big_2(const GentryPoly& a, const GentryPoly& b) const;
};

