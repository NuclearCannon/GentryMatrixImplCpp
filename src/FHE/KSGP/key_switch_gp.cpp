#include "FHE/key_switch_gp.hpp"




KeySwitchKeyGP::KeySwitchKeyGP(std::vector<std::pair<GentryPoly, GentryPoly>> cts, uint64_t qo):
    cts_(cts), qo_(qo)
{

}

KeySwitchKeyGP KeySwitchKeyGP::ksk_gen(const GentryPoly& sk_from, const GentryPoly& sk_to, uint64_t qo, const GentryPolyCtx& ctx)
{
    // 这里假设sk_from, sk_to都是含qo模数的，且sk_from已经乘过qo
    auto mods = sk_to.moduli();
    GentryPoly sk_from_cp = sk_from;
    std::vector<std::pair<GentryPoly, GentryPoly>> cts;
    for(auto mod : mods)
    {
        if(mod == qo)continue;
        cts.push_back(encrypt_gp(sk_from_cp, sk_to, ctx));
        GentryPoly::mul_scalar(sk_from_cp, sk_from_cp, mod);
    }
    return KeySwitchKeyGP(std::move(cts), qo);
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::key_switch_big_1(const GentryPoly &a) const
{
    throw std::runtime_error("key_switch_big_1尚未实现");
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::key_switch_big_2(const GentryPoly& a, const GentryPoly& b) const
{
    throw std::runtime_error("key_switch_big_2尚未实现");
}


