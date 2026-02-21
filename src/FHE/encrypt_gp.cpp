#include "FHE/encrypt_gp.hpp"

std::pair<GentryPoly, GentryPoly> encrypt_gp(
    const GentryPoly& msg,
    const GentryPoly& sk,
    const GentryPolyCtx& ctx
)
{
    assert(msg.like(sk));
    // 获取sk和
    GentryPoly a = GentryPoly::uniform(msg.n(), msg.p(), msg.moduli());
    GentryPoly e = GentryPoly::dg(msg.n(), msg.p(), msg.moduli());
    // b = m+e-as
    
    GentryPoly a_ntt = a;
    a_ntt.ntt(ctx);
    GentryPoly b = sk;
    b.ntt(ctx);    // b = NTT(sk)
    GentryPoly::mul(b, a_ntt, sk);      // b = 
    b.intt(ctx);
    GentryPoly::sub(b, msg, b);     // b = msg - b
    GentryPoly::add(b, b, e);       // b = b+e


    return {
        std::move(a),
        std::move(b)
    };
}

GentryPoly decrypt_gp(
    const GentryPoly& cta,
    const GentryPoly& ctb,
    const GentryPoly& sk,
    const GentryPolyCtx& ctx
)
{
    assert(cta.like(ctb));
    assert(cta.like(sk));
    // m = as+b
    GentryPoly m = cta;
    GentryPoly sk_ntt = sk;
    m.ntt(ctx);
    sk_ntt.ntt(ctx);
    GentryPoly::mul(m, m, sk_ntt);
    m.intt(ctx);
    GentryPoly::add(m, m, ctb);
    return m;
}