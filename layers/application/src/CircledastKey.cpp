#include "application.hpp"
#include <iostream>
#include "FHE/circledast.hpp"

CircledastKey::CircledastKey(std::unique_ptr<KeySwitchKeyGP> ksk1, std::unique_ptr<KeySwitchKeyGP> ksk2):
    ksk1_(std::move(ksk1)),
    ksk2_(std::move(ksk2))
{

}

CircledastKey CircledastKey::gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx)
{
    auto [k1, k2] = create_ksks_for_circledast_ct(*sk.data_, qo, ctx);
    std::unique_ptr<KeySwitchKeyGP> k1p = std::make_unique<KeySwitchKeyGP>(std::move(k1));
    std::unique_ptr<KeySwitchKeyGP> k2p = std::make_unique<KeySwitchKeyGP>(std::move(k2));
    return CircledastKey(std::move(k1p), std::move(k2p));

}

Ciphertext CircledastKey::run(const Ciphertext& u, const Ciphertext& v, const GentryPolyCtx& ctx) const
{
    auto [wa, wb] = circledast_ct(
        *u.a_, *u.b_, *v.a_, *v.b_, *ksk1_, *ksk2_, ctx 
    );
    std::unique_ptr<GentryPoly> wap = std::make_unique<GentryPoly>(std::move(wa));
    std::unique_ptr<GentryPoly> wbp = std::make_unique<GentryPoly>(std::move(wb));
    return Ciphertext(std::move(wap), std::move(wbp));

}

Ciphertext CircledastKey::run_cp(const Ciphertext& u, const Plaintext& v, const GentryPolyCtx& ctx) const
{
    auto [wa, wb] = circledast_cp(
        *u.a_, *u.b_, *v.data_, ctx 
    );
    std::unique_ptr<GentryPoly> wap = std::make_unique<GentryPoly>(std::move(wa));
    std::unique_ptr<GentryPoly> wbp = std::make_unique<GentryPoly>(std::move(wb));
    return Ciphertext(std::move(wap), std::move(wbp));
}
Ciphertext CircledastKey::run_pc(const Plaintext& u, const Ciphertext& v, const GentryPolyCtx& ctx) const
{
    auto [wa, wb] = circledast_pc(
        *u.data_, *v.a_, *v.b_, *ksk1_, ctx 
    );
    std::unique_ptr<GentryPoly> wap = std::make_unique<GentryPoly>(std::move(wa));
    std::unique_ptr<GentryPoly> wbp = std::make_unique<GentryPoly>(std::move(wb));
    return Ciphertext(std::move(wap), std::move(wbp));
}


void CircledastKey::test_pt_circledast_end2end()
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

    // 准备两个明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(5, n, p);

    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);

    ComplexMatrixGroup mat2 = ComplexMatrixGroup::random(5, n, p);

    Plaintext pt2 = Plaintext::from_cmat(mat2, mods, delta);
    Plaintext pt3 = pt1.circledast(pt2, ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat3 = pt3.to_cmat(delta * delta);
    // 对mat3乘以n以弥补circledast运算的副作用
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        mat3.at(w, x, y) *= n;
    }
    // 进行一次更标准的AB*乘法
    ComplexMatrixGroup mat4 = mat1.matmul_ABT(mat2);
    // 比较mat3, mat4的内容
    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat3.at(w,x,y) - mat4.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_pt_circledast_end2end: max_diff="<<max_diff<<std::endl;

    
}

void CircledastKey::test_ct_circledast_end2end(bool cuda)
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


    // 准备两个明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(5, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    ComplexMatrixGroup mat2 = ComplexMatrixGroup::random(5, n, p);
    Plaintext pt2 = Plaintext::from_cmat(mat2, mods, delta);
    // 准备一个私钥
    SecretKey sk(n, p, mods);
    // 分别加密成密文
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    Ciphertext ct2 = Ciphertext::encrypt(pt2, sk, ctx);
    // 生成circledast ksk
    CircledastKey key = CircledastKey::gen(sk, qo, ctx);
    Ciphertext ct3 = (cuda) ? key.run(ct1.to_cuda(), ct2.to_cuda(), ctx) : key.run(ct1, ct2, ctx);
    Plaintext pt3 = ct3.decrypt(sk, ctx);
    // 恢复成矩阵
    ComplexMatrixGroup mat3 = pt3.to_cmat(delta * delta);
    // 对mat3乘以n以弥补circledast运算的副作用
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        mat3.at(w, x, y) *= n;
    }

    // 进行一次更标准的AB*乘法
    ComplexMatrixGroup mat4 = mat1.matmul_ABT(mat2);
    // 比较mat3, mat4的内容
    double max_diff = 0;

    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat3.at(w,x,y) - mat4.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_circledast_end2end(CCMM): max_diff="<<max_diff<<std::endl;

    {
        Ciphertext ct_cp = key.run_cp(ct1, pt2, ctx);
        ComplexMatrixGroup mat_cp = ct_cp.decrypt(sk, ctx).to_cmat(delta * delta);
        // 对mat_cp乘以n以弥补circledast运算的副作用
        for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
        {
            mat_cp.at(w, x, y) *= n;
        }
        max_diff = 0;
        for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
        {
            complex diff = mat_cp.at(w,x,y) - mat4.at(w,x,y);
            double diff_abs = std::abs(diff);
            if(diff_abs>max_diff)max_diff = diff_abs;
        }
        std::cout << "test_ct_circledast_end2end(CPMM): max_diff="<<max_diff<<std::endl;
    }
    
    {
        Ciphertext ct_pc = key.run_pc(pt1, ct2, ctx);
        ComplexMatrixGroup mat_pc = ct_pc.decrypt(sk, ctx).to_cmat(delta * delta);
        // 对mat_pc乘以n以弥补circledast运算的副作用
        for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
        {
            mat_pc.at(w, x, y) *= n;
        }
        max_diff = 0;
        for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
        {
            complex diff = mat_pc.at(w,x,y) - mat4.at(w,x,y);
            double diff_abs = std::abs(diff);
            if(diff_abs>max_diff)max_diff = diff_abs;
        }
        std::cout << "test_ct_circledast_end2end(PCMM): max_diff="<<max_diff<<std::endl;
    }

    


    
}

