#pragma once
#include "FHE/key_switch_gp.hpp"

std::pair<KeySwitchKeyGP, KeySwitchKeyGP> create_ksks_for_circledast_ct(
    const GentryPoly& sk, uint64_t qo, const GentryPolyCtx &ctx);

std::pair<GentryPoly, GentryPoly> circledast_ct(
    const GentryPoly& ua,
    const GentryPoly& ub,
    const GentryPoly& va,
    const GentryPoly& vb,
    const KeySwitchKeyGP& ksk1,
    const KeySwitchKeyGP& ksk2,
    const GentryPolyCtx& ctx
);