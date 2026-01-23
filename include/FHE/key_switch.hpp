#pragma once
#include "ziq_array.hpp"

class KeySwitchingKey
{
private:
    fmpz_scalar B_, Q_, q_, q_half_, Q_half_;
    int L_;
    const ZiqArrayContext *ctx_q_, *ctx_Q_;
    // sk_from在sk_to下的各段密文
    std::vector<std::pair<ZiqArray, ZiqArray>> cts;


public:
    // 构造时，两个sk可以被以任意的方法传入
    KeySwitchingKey(
        const ZiqArray& sk_from,
        const ZiqArray& sk_to,
        const ZiqArrayContext* ctx_q,
        const ZiqArrayContext* ctx_Q,
        const fmpz_t B,
        int L
    );

    std::pair<ZiqArray, ZiqArray> key_switch_big_1(
        const ZiqArray& ct_a_from
    ) const;

    std::pair<ZiqArray, ZiqArray> key_switch_big_2(
        const ZiqArray& ct_a_from,
        const ZiqArray& ct_b_from
    ) const;

};