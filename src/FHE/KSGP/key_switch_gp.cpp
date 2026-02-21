#include "FHE/key_switch_gp.hpp"




KeySwitchKeyGP::KeySwitchKeyGP(std::vector<std::pair<GentryPoly, GentryPoly>> cts, uint64_t qo):
    cts_(cts), qo_(qo)
{

}

KeySwitchKeyGP KeySwitchKeyGP::ksk_gen(const GentryPoly& sk_from, const GentryPoly& sk_to, uint64_t qo)
{
    throw std::runtime_error("ksk_gen尚未实现");
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::key_switch_big_1(const GentryPoly &a) const
{
    throw std::runtime_error("key_switch_big_1尚未实现");
}
std::pair<GentryPoly, GentryPoly> KeySwitchKeyGP::key_switch_big_2(const GentryPoly& a, const GentryPoly& b) const
{
    throw std::runtime_error("key_switch_big_2尚未实现");
}


