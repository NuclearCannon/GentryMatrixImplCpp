#include "FHE/circledast.hpp"

std::pair<CRTArray, CRTArray> circledast_ct(
    const CRTArray& ua,
    const CRTArray& ub,
    const CRTArray& va,
    const CRTArray& vb,
    const KeySwitchKey64CRT& ksk1,
    const KeySwitchKey64CRT& ksk2
)
{
    // 需要先转化为half形式才能circledast
    CRTArray auh = ua.iw_ntt();
    CRTArray avh = va.iw_ntt();
    CRTArray buh = ub.iw_ntt();
    CRTArray bvh = vb.iw_ntt();

    CRTArray auav = auh.circledast(avh).iw_intt();
    CRTArray aubv = auh.circledast(bvh).iw_intt();
    CRTArray buav = buh.circledast(avh).iw_intt();
    CRTArray bubv = buh.circledast(bvh).iw_intt();
    auto [buav0, buav1] = ksk1.key_switch_big_1(buav);
    auto [auav0, auav1] = ksk2.key_switch_big_1(auav);
    return std::make_pair(
        aubv.add(buav0).add(auav0),
        bubv.add(buav1).add(auav1)
    );
}