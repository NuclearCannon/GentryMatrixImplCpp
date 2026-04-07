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
        *(sk.data_),
        ctx
    );
    return Ciphertext(
        std::make_unique<GentryPoly>(std::move(a)),
        std::make_unique<GentryPoly>(std::move(b))
    );
}


Plaintext Ciphertext::decrypt(const SecretKey& sk, const GentryPolyCtx& ctx)
{
    if(a_->is_cuda())
    {
        return Plaintext(std::make_unique<GentryPoly>(decrypt_gp(a_->to_cpu(), b_->to_cpu(), *(sk.data_), ctx)));
    }
    else
    {
        return Plaintext(std::make_unique<GentryPoly>(decrypt_gp(*a_, *b_, *(sk.data_), ctx)));
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
    SecretKey sk(n, p, mods);
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
