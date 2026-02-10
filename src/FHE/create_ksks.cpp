#include "FHE/key_switch_64.hpp"

std::pair<KeySwitchKey64CRT, KeySwitchKey64CRT> create_ksks_for_circledast_ct(const CRTArray& sk, u64 qo, u64 qor)
{
    CRTArray sk1 = sk.transpose().conj().w_inv();
    CRTArray sk2 = sk1.mul_poly(sk);
    return {
        KeySwitchKey64CRT::ksk_gen(sk1, sk, qo, qor),
        KeySwitchKey64CRT::ksk_gen(sk2, sk, qo, qor)
    };
}

