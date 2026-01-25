#include "FHE/key_switch.hpp"
#include <memory>

std::pair<std::unique_ptr<KeySwitchingKey>, std::unique_ptr<KeySwitchingKey>> create_ksks_for_circledast(
    const ZiqArray& sk,
    const ZiqArrayContext* ctx_q,
    const ZiqArrayContext* ctx_Q,
    const fmpz_t B,
    int L
);

std::pair<ZiqArray, ZiqArray> circledast_ct(
    const ZiqArray& au,
    const ZiqArray& bu,
    const ZiqArray& av,
    const ZiqArray& bv,
    const KeySwitchingKey& ksk1,
    const KeySwitchingKey& ksk2
);

std::pair<ZiqArray, ZiqArray> circledast_ct_timer(
    const ZiqArray& au,
    const ZiqArray& bu,
    const ZiqArray& av,
    const ZiqArray& bv,
    const KeySwitchingKey& ksk1,
    const KeySwitchingKey& ksk2
);