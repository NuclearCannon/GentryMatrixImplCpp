#include "application.hpp"
#include <iostream>
#include "FHE/encrypt_gp.hpp"

Ciphertext::Ciphertext(std::unique_ptr<GentryPoly> a, std::unique_ptr<GentryPoly> b):
    a_(std::move(a)), b_(std::move(b))
{

}

Ciphertext Ciphertext::encrypt(const Plaintext& pt, const SecretKey& sk, const GentryPolyCtx& ctx)
{
    auto [a, b] = encrypt_gp(
        *(pt.data_),
        sk.as_poly(pt.data_->moduli()),
        ctx
    );
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}

Ciphertext Ciphertext::zeros(int n, int p, std::vector<uint64_t> mods)
{
    return Ciphertext(
        std::make_unique<GentryPoly>(GentryPoly::zeros(n, p, mods)),
        std::make_unique<GentryPoly>(GentryPoly::zeros(n, p, mods))
    );
}

Plaintext Ciphertext::decrypt(const SecretKey& sk, const GentryPolyCtx& ctx)
{
    if(a_->is_cuda())
    {
        return Plaintext(std::make_unique<GentryPoly>(decrypt_gp(a_->to_cpu(), b_->to_cpu(), sk.as_poly(a_->moduli()), ctx)));
    }
    else
    {
        return Plaintext(std::make_unique<GentryPoly>(decrypt_gp(*a_, *b_, sk.as_poly(a_->moduli()), ctx)));
    }
    

}

Ciphertext Ciphertext::to_cuda() const
{
    return Ciphertext(
        std::make_unique<GentryPoly>(a_->to_cuda()),
        std::make_unique<GentryPoly>(b_->to_cuda())
    );
}

void Ciphertext::test_ct_encrypt_and_decrypt()
{
    int n = 16;
    int p = 17;
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    vec64 roots = {6, 11, 6};
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
    };
    GentryPolyCtx ctx(n, p, qrp);

    double delta = 10000;

    ComplexMatrixGroup mat = ComplexMatrixGroup::random(5, n, p);
    Plaintext pt = Plaintext::from_cmat(mat, mods, delta);
    SecretKey sk(n, p);
    Ciphertext ct = Ciphertext::encrypt(pt, sk, ctx);
    Plaintext pt2 = ct.decrypt(sk, ctx);
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - mat.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_encrypt_and_decrypt: max_diff="<<max_diff<<std::endl;
}

Ciphertext Ciphertext::add_pt(const Plaintext& pt, const GentryPolyCtx& ctx) const
{
    GentryPoly a = *a_;
    GentryPoly b = GentryPoly::zeros_like(*b_);
    GentryPoly::add(b, *b_, *pt.data_);
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}

Ciphertext Ciphertext::add(const Ciphertext& other, const GentryPolyCtx& ctx) const
{
    GentryPoly a = GentryPoly::zeros_like(*a_);
    GentryPoly b = GentryPoly::zeros_like(*b_);
    GentryPoly::add(a, *a_, *other.a_);
    GentryPoly::add(b, *b_, *other.b_);
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}

void Ciphertext::add_(const Ciphertext& other, const GentryPolyCtx& ctx)
{
    GentryPoly::add(*a_, *a_, *other.a_);
    GentryPoly::add(*b_, *b_, *other.b_);
}
Ciphertext Ciphertext::mul_pt(const Plaintext& pt, const GentryPolyCtx& ctx) const
{
    GentryPoly a = *a_;
    GentryPoly b = *b_;
    GentryPoly m = *pt.data_;
    a.ntt(ctx);
    b.ntt(ctx);
    m.ntt(ctx);
    GentryPoly a_result = GentryPoly::zeros_like(a);
    GentryPoly b_result = GentryPoly::zeros_like(b);
    GentryPoly::mul(a_result, m, a);
    GentryPoly::mul(b_result, m, b);
    a_result.intt(ctx);
    b_result.intt(ctx);
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a_result)),
        std::make_unique<GentryPoly>(std::move(b_result))
    );
}

void Ciphertext::test_ct_add_pt()
{
    int n = 8, p = 5;
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    vec64 roots = {6, 11, 6};
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
    };
    GentryPolyCtx ctx(n, p, qrp);

    double delta = 1000;

    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(3, n, p);
    ComplexMatrixGroup mat2 = ComplexMatrixGroup::random(3, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    Plaintext pt2 = Plaintext::from_cmat(mat2, mods, delta);
    SecretKey sk(n, p);
    Ciphertext ct = Ciphertext::encrypt(pt1, sk, ctx);
    Ciphertext ct_result = ct.add_pt(pt2, ctx);
    Plaintext pt_result = ct_result.decrypt(sk, ctx);
    ComplexMatrixGroup mat_result = pt_result.to_cmat(delta);

    ComplexMatrixGroup expected(n, p);
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        expected.at(w,x,y) = mat1.at(w,x,y) + mat2.at(w,x,y);
    }

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat_result.at(w,x,y) - expected.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_add_pt: max_diff="<<max_diff<<std::endl;
}

void Ciphertext::test_ct_add()
{
    int n = 8, p = 5;
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    vec64 roots = {6, 11, 6};
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
    };
    GentryPolyCtx ctx(n, p, qrp);

    double delta = 1000;

    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(3, n, p);
    ComplexMatrixGroup mat2 = ComplexMatrixGroup::random(3, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    Plaintext pt2 = Plaintext::from_cmat(mat2, mods, delta);
    SecretKey sk(n, p);
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    Ciphertext ct2 = Ciphertext::encrypt(pt2, sk, ctx);
    Ciphertext ct_result = ct1.add(ct2, ctx);
    Plaintext pt_result = ct_result.decrypt(sk, ctx);
    ComplexMatrixGroup mat_result = pt_result.to_cmat(delta);

    ComplexMatrixGroup expected(n, p);
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        expected.at(w,x,y) = mat1.at(w,x,y) + mat2.at(w,x,y);
    }

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat_result.at(w,x,y) - expected.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_add: max_diff="<<max_diff<<std::endl;
}

void Ciphertext::test_ct_mul_pt()
{
    int n = 8, p = 5;
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    vec64 roots = {6, 11, 6};
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
    };
    GentryPolyCtx ctx(n, p, qrp);

    double delta = 1000;

    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(3, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    SecretKey sk(n, p);
    Ciphertext ct = Ciphertext::encrypt(pt1, sk, ctx);
    Ciphertext ct_result = ct.mul_pt(pt1, ctx);
    Plaintext pt_result = ct_result.decrypt(sk, ctx);
    ComplexMatrixGroup mat_result = pt_result.to_cmat(delta * delta);

    ComplexMatrixGroup expected(n, p);
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        expected.at(w,x,y) = mat1.at(w,x,y) * mat1.at(w,x,y);
    }

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat_result.at(w,x,y) - expected.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_ct_mul_pt: max_diff="<<max_diff<<std::endl;
}
