#include "FHE/key_switch.hpp"
#include "FHE/encrypt.hpp"

KeySwitchingKey::KeySwitchingKey(
    const ZiqArray& sk_from,
    const ZiqArray& sk_to,
    const ZiqArrayContext* ctx_q,
    const ZiqArrayContext* ctx_Q,
    const fmpz_t B,
    int L
):
    B_(B),
    Q_(ctx_Q->q()),
    q_(ctx_q->q()),
    L_(L),
    ctx_q_(ctx_q),
    ctx_Q_(ctx_Q)
{

    const fmpz* q = q_.raw();
    const fmpz* Q = Q_.raw();
    int length = sk_from.data().len();
    assert(sk_to.data().len() == length);
    assert(ctx_q->get_size() == length);
    assert(ctx_Q->get_size() == length);

    fmpz_vector sk_from_Q = sk_from.data();
    ZiqArray sk_to_Q = sk_to.ctx_switch(ctx_Q);

    fmpz_fdiv_q_si(q_half_.raw(), q, 2);
    fmpz_fdiv_q_si(Q_half_.raw(), Q, 2);


    for(int i=0;i<L;i++)
    {
        // 获取B^i * sk_from * Q // q，存入bis中
        fmpz_vector bis(length);
        fmpz_t B_pow_i, temp_prod;
        fmpz_init(B_pow_i);
        fmpz_init(temp_prod);
        // 计算 B^i
        if (i == 0) {
            fmpz_one(B_pow_i);
        } else {
            fmpz_pow_ui(B_pow_i, B, (ulong)i);
        }
        for (slong j = 0; j < length; ++j) {
            // temp_prod = B^i * sk_from[j]
            fmpz_mul(temp_prod, B_pow_i, sk_from_Q[j]);
            // temp_prod %= q
            fmpz_mod(temp_prod, temp_prod, q);
            // temp_prod = temp_prod * Q
            fmpz_mul(temp_prod, temp_prod, Q);
            // + q//2
            fmpz_add(temp_prod, temp_prod, q_half_.raw());
            // bis[j] = trunc( temp_prod / q ) —— 忽略余数
            fmpz_fdiv_q(bis[j], temp_prod, q);
        }
        fmpz_clear(B_pow_i);
        fmpz_clear(temp_prod);
        // 将这一份结果加密成密文
        auto [a, b] = encrypt(ZiqArray(std::move(bis), ctx_Q), sk_to_Q);
        auto a_ntt = a.iw_ntt().xy_ntt();
        auto b_ntt = b.iw_ntt().xy_ntt();
        cts.push_back(std::make_pair(std::move(a_ntt), std::move(b_ntt)));
    }
}

std::pair<ZiqArray, ZiqArray> KeySwitchingKey::key_switch_big_2(
    const ZiqArray& ct_a_from,
    const ZiqArray& ct_b_from
) const
{
    assert(ct_a_from.ctx() == ctx_q_);
    assert(ct_b_from.ctx() == ctx_q_);

    auto [a1, b1] = key_switch_big_1(ct_a_from);
    auto b2 = b1.add(ct_b_from);
    return std::make_pair(std::move(a1), std::move(b2));
}

std::pair<ZiqArray, ZiqArray> KeySwitchingKey::key_switch_big_1(
    const ZiqArray& ct_a_from
) const
{
    assert(ct_a_from.ctx() == ctx_q_);

    // 对ct_a_from沿B分解
    fmpz_vector a = ct_a_from.data();
    int len = a.len();
    assert(len == ctx_Q_->get_size());

    std::vector<ZiqArray> a_split;
    while(1)
    {
        // 检查a是不是还有0？
        if (_fmpz_vec_is_zero(a.raw(), len)) break;
        // 进行一次取模
        fmpz_vector t_mod(len);
        // let t_mod = a % B
        _fmpz_vec_scalar_mod_fmpz(t_mod.raw(), a.raw(), len, B_.raw());
        // let a = a // B
        _fmpz_vec_scalar_tdiv_q_fmpz(a.raw(), a.raw(), len, B_.raw());
        a_split.push_back(ZiqArray(std::move(t_mod), ctx_Q_));
    }
    // 检查长度
    if (a_split.size() > L_)
    {
        throw std::invalid_argument("ct_a_from不能被完整分解！\n");
    }
    // 逐个相乘
    fmpz_vector a_sum_ntt = fmpz_vector::zeros(len);
    fmpz_vector b_sum_ntt = fmpz_vector::zeros(len);

    for(int i=0;i<a_split.size();i++)
    {
        const auto& [sk_a_ntt, sk_b_ntt] = cts[i];
        ZiqArray a_piece_ntt = a_split[i].iw_ntt().xy_ntt();
        ZiqArray skaa = a_piece_ntt.mul(sk_a_ntt);
        ZiqArray skba = a_piece_ntt.mul(sk_b_ntt);
        _fmpz_vec_add(a_sum_ntt.raw(), a_sum_ntt.raw(), skaa.data().raw(), len);
        _fmpz_vec_add(b_sum_ntt.raw(), b_sum_ntt.raw(), skba.data().raw(), len);
    }
    // 取模Q
    _fmpz_vec_scalar_mod_fmpz(a_sum_ntt.raw(), a_sum_ntt.raw(), len, Q_.raw());
    _fmpz_vec_scalar_mod_fmpz(b_sum_ntt.raw(), b_sum_ntt.raw(), len, Q_.raw());
    // intt
    fmpz_vector a_sum = ctx_Q_->iw_intt(ctx_Q_->xy_intt(a_sum_ntt));
    fmpz_vector b_sum = ctx_Q_->iw_intt(ctx_Q_->xy_intt(b_sum_ntt));
    // 取模Q
    _fmpz_vec_scalar_mod_fmpz(a_sum.raw(), a_sum.raw(), len, Q_.raw());
    _fmpz_vec_scalar_mod_fmpz(b_sum.raw(), b_sum.raw(), len, Q_.raw());
    // 降低模数
    _fmpz_vec_scalar_mul_fmpz(a_sum.raw(), a_sum.raw(), len, q_.raw());
    _fmpz_vec_scalar_mul_fmpz(b_sum.raw(), b_sum.raw(), len, q_.raw());
    // +Q/2
    for(int i=0;i<len;i++)
    {
        fmpz_add(a_sum[i], a_sum[i], Q_half_.raw());
        fmpz_add(b_sum[i], b_sum[i], Q_half_.raw());
    }
    _fmpz_vec_scalar_tdiv_q_fmpz(a_sum.raw(), a_sum.raw(), len, Q_.raw());
    _fmpz_vec_scalar_tdiv_q_fmpz(b_sum.raw(), b_sum.raw(), len, Q_.raw());

    return std::make_pair(
        ZiqArray(std::move(a_sum), ctx_q_), ZiqArray(std::move(b_sum), ctx_q_)
    );

}