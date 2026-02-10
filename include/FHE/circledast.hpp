#include "key_switch_64.hpp"

std::pair<CRTArray, CRTArray> circledast_ct(
    const CRTArray& ua,
    const CRTArray& ub,
    const CRTArray& va,
    const CRTArray& vb,
    const KeySwitchKey64CRT& ksk1,
    const KeySwitchKey64CRT& ksk2
);