#include "application.hpp"
#include <iostream>

RotateKey::RotateKey(std::unique_ptr<KeySwitchKeyGP> ksk, int x_bias, int y_bias, int w_bias):
    ksk_(std::move(ksk)), x_bias_(x_bias), y_bias_(y_bias), w_bias_(w_bias)
{

}

RotateKey RotateKey::gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods, int x_bias, int y_bias, int w_bias)
{
    // 注意，sk在用于生成KSK前需要扩展模数
    GentryPoly sk_to = sk.as_poly(mods);
    GentryPoly sk_from = sk_to.rotate(x_bias, y_bias, w_bias);
    sk_from.moduli_extend_mult(qo);
    sk_to.moduli_extend_unsafe(qo);

    return RotateKey(
        std::make_unique<KeySwitchKeyGP>(
            KeySwitchKeyGP::ksk_gen(
                sk_from,
                sk_to,
                qo, ctx
            )
        ), 
        x_bias, y_bias, w_bias
    );
}

Ciphertext RotateKey::run(const Ciphertext& src, const GentryPolyCtx& ctx) const
{
    GentryPoly a = (*src.a_).rotate(x_bias_, y_bias_, w_bias_);
    GentryPoly b = (*src.b_).rotate(x_bias_, y_bias_, w_bias_);
    // ks
    auto [a2, b2] = ksk_->key_switch_big_2(a, b, ctx);

    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a2)), 
        std::make_unique<GentryPoly>(std::move(b2))
    );
}


void RotateKey::test_pt_rotate()
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
    // 确定rotate步长
    int x_bias = 5;
    int y_bias = 3;
    int w_bias = 2;

    GentryPolyCtx ctx(n, p, qrp);

    // 准备明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(5, n, p);
    // for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++) mat1.at(w, x, y) = w*x;
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    Plaintext pt2 = pt1.rotate(x_bias, y_bias, w_bias);
    // 恢复成矩阵
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    // mat1.print("mat1");
    // mat2.print("mat2");

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - mat1.at((w+w_bias)%(p-1),(x+x_bias)%n,(y+y_bias)%n); // 在这里进行编码前层面上的转置操作
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_pt_rotate: max_diff="<<max_diff<<std::endl;
}

void RotateKey::test_ct_rotate()
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
    // 确定rotate步长
    int x_bias = 5;
    int y_bias = 3;
    int w_bias = 1;

    GentryPolyCtx ctx(n, p, qrp);

    // 准备明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(5, n, p);
    // for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++) mat1.at(w, x, y) = w*x;
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    SecretKey sk(n, p);
    RotateKey key = RotateKey::gen(sk, qo, ctx, mods, x_bias, y_bias, w_bias);
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    Ciphertext ct2 = key.run(ct1, ctx);
    Plaintext pt2 = ct2.decrypt(sk, ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    // mat1.print("mat1");
    // mat2.print("mat2");

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - mat1.at((w+w_bias)%(p-1),(x+x_bias)%n,(y+y_bias)%n); // 在这里进行编码前层面上的转置操作
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_rotate: max_diff="<<max_diff<<std::endl;
}