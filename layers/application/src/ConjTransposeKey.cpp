#include <iostream>
#include "application.hpp"

ConjTransposeKey::ConjTransposeKey(std::unique_ptr<KeySwitchKeyGP> ksk):
    ksk_(std::move(ksk))
{
    
}

ConjTransposeKey ConjTransposeKey::gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods)
{
    // 注意，sk在用于生成KSK前需要扩展模数

    GentryPoly sk_to = sk.as_poly(mods);
    GentryPoly sk_from = sk_to.w_inv().transpose().conj();
    sk_from.moduli_extend_mult(qo);
    sk_to.moduli_extend_unsafe(qo);

    return ConjTransposeKey(
        std::make_unique<KeySwitchKeyGP>(
            KeySwitchKeyGP::ksk_gen(
                sk_from,
                sk_to,
                qo, ctx
            )
        )
    );
}

Ciphertext ConjTransposeKey::run(const Ciphertext& src, const GentryPolyCtx& ctx) const
{
    GentryPoly a = (*src.a_).w_inv().transpose().conj();
    GentryPoly b = (*src.b_).w_inv().transpose().conj();
    // ks
    auto [a2, b2] = ksk_->key_switch_big_2(a, b, ctx);

    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a2)), 
        std::make_unique<GentryPoly>(std::move(b2))
    );
}

void ConjTransposeKey::test_pt_transpose()
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
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++) mat1.at(w, x, y) = w*x;
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    Plaintext pt2 = pt1.conj_transpose(ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - std::conj(mat1.at(w,y,x)); // 在这里进行编码前层面上的转置操作
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_pt_transpose: max_diff="<<max_diff<<std::endl;
}

void ConjTransposeKey::test_ct_transpose()
{
    // 准备各个参数
    int n = 8;
    int p = 17;
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
    // 准备一个私钥
    SecretKey sk(n, p);
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    ConjTransposeKey key = ConjTransposeKey::gen(sk, qo, ctx, mods);
    Ciphertext ct2 = key.run(ct1, ctx);
    Plaintext pt2 = ct2.decrypt(sk, ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - std::conj(mat1.at(w,y,x)); // 在这里进行编码前层面上的转置操作
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_transpose: max_diff="<<max_diff<<std::endl;
}