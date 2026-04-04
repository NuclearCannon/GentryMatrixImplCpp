#include "FHE/key_switch_gp.hpp"
#include "FHE/circledast.hpp"


std::pair<KeySwitchKeyGP, KeySwitchKeyGP> create_ksks_for_circledast_ct(
    const GentryPoly& sk, uint64_t qo, const GentryPolyCtx &ctx)
{
    


    GentryPoly sk1 = sk.transpose().conj().w_inv();
    // sk2 = sk1 多项式乘 sk
    GentryPoly sk_ntt = sk, sk1_ntt = sk1, sk2 = GentryPoly::zeros_like(sk1);
    sk_ntt.ntt(ctx);
    sk1_ntt.ntt(ctx);
    GentryPoly::mul(sk2, sk1_ntt, sk_ntt);
    sk2.intt(ctx);

    sk1.moduli_extend_mult(qo);
    sk2.moduli_extend_mult(qo);

    GentryPoly sk_qo = sk;
    sk_qo.moduli_extend_unsafe(qo);

    return {
        KeySwitchKeyGP::ksk_gen(sk1, sk_qo, qo, ctx),
        KeySwitchKeyGP::ksk_gen(sk2, sk_qo, qo, ctx)
    };
}

