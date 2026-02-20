#pragma once
#include "GentryPoly.hpp"

std::pair<GentryPoly, GentryPoly> encrypt_gp(
    const GentryPoly& msg,
    const GentryPoly& sk,
    const GentryPolyCtx& ctx
);

GentryPoly decrypt_gp(
    const GentryPoly& cta,
    const GentryPoly& ctb,
    const GentryPoly& sk,
    const GentryPolyCtx& ctx
);