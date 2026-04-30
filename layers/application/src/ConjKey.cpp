#include "application.hpp"
#include <iostream>


ConjKey::ConjKey(std::unique_ptr<KeySwitchKeyGP> ksk):
    ksk_(std::move(ksk))
{

}



ConjKey ConjKey::gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods)
{
    GentryPoly sk_to = sk.as_poly(mods);
    throw std::runtime_error("未实现");
}


Ciphertext ConjKey::run(const Ciphertext& src, const GentryPolyCtx& ctx) const
{
    throw std::runtime_error("未实现");
}


void ConjKey::test_pt_conj()
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
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++) mat1.at(w, x, y) = complex(w*x, w*y);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    Plaintext pt2 = pt1.conj(ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - std::conj(mat1.at(w,x,y)); // 在这里进行编码前层面上的转置操作
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    // mat1.print("mat1");
    // mat2.print("mat2");
    std::cout << "test_pt_conj: max_diff="<<max_diff<<std::endl;
}