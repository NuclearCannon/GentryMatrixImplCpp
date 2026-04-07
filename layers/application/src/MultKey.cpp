#include <iostream>
#include "application.hpp"

static GentryPoly polymul(GentryPoly x, GentryPoly y, const GentryPolyCtx& ctx)
{
    x.ntt(ctx);
    y.ntt(ctx);
    GentryPoly::mul(x, x, y);
    x.intt(ctx);
    return x;
}


MultKey::MultKey(std::unique_ptr<KeySwitchKeyGP> ksk):
    ksk_(std::move(ksk))
{
    
}

MultKey MultKey::gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx)
{
    GentryPoly sk_ntt = *sk.data_;
    sk_ntt.ntt(ctx);
    GentryPoly sk2 = GentryPoly::zeros_like(sk_ntt);
    GentryPoly::mul(sk2, sk_ntt, sk_ntt);
    sk2.intt(ctx);
    // 注意，sk在用于生成KSK前需要扩展模数
    sk2.moduli_extend_mult(qo);
    GentryPoly sk_to = (*sk.data_);
    sk_to.moduli_extend_unsafe(qo);
    return MultKey(
        std::make_unique<KeySwitchKeyGP>(
            KeySwitchKeyGP::ksk_gen(
                sk2,
                sk_to,
                qo, ctx
            )
        )
    );
}

Ciphertext MultKey::run(const Ciphertext& src1, const Ciphertext& src2, const GentryPolyCtx& ctx) const
{
    GentryPoly a1 = *src1.a_;
    GentryPoly b1 = *src1.b_;
    GentryPoly a2 = *src2.a_;
    GentryPoly b2 = *src2.b_;
    a1.ntt(ctx);
    b1.ntt(ctx);
    a2.ntt(ctx);
    b2.ntt(ctx);

    GentryPoly x3 = GentryPoly::zeros_like(a1); // s^2
    GentryPoly a3 = GentryPoly::zeros_like(a1); // s
    GentryPoly b3 = GentryPoly::zeros_like(a1); // 1

    // b3 = b1 * b2
    GentryPoly::mul(b3, b1, b2);
    // a3 = a1 * b2 + a2 * b1
    GentryPoly::mul(a3, a1, b2);
    GentryPoly::mul(x3, a2, b1);    // 把x3当成临时变量用
    GentryPoly::add(a3, a3, x3);
    // x3 = a1 * a2
    GentryPoly::mul(x3, a1, a2);
    x3.intt(ctx);
    a3.intt(ctx);
    b3.intt(ctx);

    // ks
    auto [a_ks, b_ks] = ksk_->key_switch_big_1(x3, ctx);

    // a3, b3 += a_ks, b_ks
    GentryPoly::add(a3, a3, a_ks);
    GentryPoly::add(b3, b3, b_ks);

    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a3)), 
        std::make_unique<GentryPoly>(std::move(b3))
    );
}


void MultKey::test_ct_mult()
{
    // 准备各个参数
    int n = 8;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };
    double delta = 10000;

    GentryPolyCtx ctx(n, p, qrp);


    // 准备明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(5, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    ComplexMatrixGroup mat2 = ComplexMatrixGroup::random(5, n, p);
    Plaintext pt2 = Plaintext::from_cmat(mat2, mods, delta);
    // 准备一个私钥
    SecretKey sk(n, p, mods);
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    Ciphertext ct2 = Ciphertext::encrypt(pt2, sk, ctx);
    MultKey key = MultKey::gen(sk, qo, ctx);
    Ciphertext ct3 = key.run(ct1, ct2, ctx);
    Plaintext pt3 = ct3.decrypt(sk, ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat3 = pt3.to_cmat(delta * delta);
    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat3.at(w,x,y) - mat1.at(w,x,y) * mat2.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_mult: max_diff="<<max_diff<<std::endl;
}