#include "FHE/circledast.hpp"

std::pair<GentryPoly, GentryPoly> circledast_ct(
    const GentryPoly& ua,
    const GentryPoly& ub,
    const GentryPoly& va,
    const GentryPoly& vb,
    const KeySwitchKeyGP& ksk1,
    const KeySwitchKeyGP& ksk2,
    const GentryPolyCtx& ctx
)
{
    // 需要先转化为half形式才能circledast
    GentryPoly auh = ua; auh.iw_ntt(ctx);
    GentryPoly avh = va; avh.iw_ntt(ctx);
    GentryPoly buh = ub; buh.iw_ntt(ctx);
    GentryPoly bvh = vb; bvh.iw_ntt(ctx);

    GentryPoly auav = GentryPoly::zeros_like(ua, ua.device());
    GentryPoly aubv = GentryPoly::zeros_like(ua, ua.device());
    GentryPoly buav = GentryPoly::zeros_like(ua, ua.device());
    GentryPoly bubv = GentryPoly::zeros_like(ua, ua.device());

    GentryPoly::circledast(auav, auh, avh);
    GentryPoly::circledast(aubv, auh, bvh);
    GentryPoly::circledast(buav, buh, avh);
    GentryPoly::circledast(bubv, buh, bvh);

    auav.iw_intt(ctx);
    aubv.iw_intt(ctx);
    buav.iw_intt(ctx);
    bubv.iw_intt(ctx);

    auto [buav0, buav1] = ksk1.key_switch_big_1(buav, ctx);
    auto [auav0, auav1] = ksk2.key_switch_big_1(auav, ctx);
    GentryPoly::add(aubv, aubv, auav0);
    GentryPoly::add(aubv, aubv, buav0);

    GentryPoly::add(bubv, bubv, auav1);
    GentryPoly::add(bubv, bubv, buav1);

    return std::make_pair(
        std::move(aubv), std::move(bubv)
    );
}