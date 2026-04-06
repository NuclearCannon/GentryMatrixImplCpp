#include "application.hpp"
#include <iostream>


Plaintext::Plaintext(std::unique_ptr<GentryPoly> data):
    data_(std::move(data))
{

}

Plaintext Plaintext::from_cmat(const ComplexMatrixGroup& cmg, const std::vector<uint64_t>& mods, double delta)
{
    // ComplexMatrixGroup encoded = cmg.encode();
    // fmpz_vector vec = cmg.encode().to_fmpz_vector(delta);
    fmpz_vector vec = cmg.to_fmpz_vector(delta);
    auto res = crt(vec, mods);  // 执行CRT分解
    size_t n = cmg.get_n(), p = cmg.get_p();
    std::unique_ptr<GentryPoly> data = std::make_unique<GentryPoly>(GentryPoly::from_coeffs(n, p, mods, res));
    return Plaintext(std::move(data));
}

ComplexMatrixGroup Plaintext::to_cmat(double delta) const 
{
    auto r = data_->to_coeffs();
    // length = 2*p*n*n
    auto len = data_->n() * data_->n() * (data_->p()-1) * 2;
    fmpz_vector vec(len);
    auto mods = data_->moduli();
    icrt(vec, r, mods);
    // 对icrt的结果进行向中心的取整
    fmpz_scalar Q = fmpz_scalar::from_ui(1);
    for(auto i:mods)
    {
        // Q *= i
        fmpz_mul_ui(Q.raw(), Q.raw(), i);
    }
    fmpz_vector vec_centered = vec.mod_centered(Q.raw());
    ComplexMatrixGroup encoded = ComplexMatrixGroup::from_fmpz_vector(vec_centered, delta, data_->n(), data_->p());
    // return encoded.decode();
    return encoded;
}

void Plaintext::test_pt_encode_and_decode()
{
    int n=8, p=5;
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    double delta = 1000;

    ComplexMatrixGroup mat(n,p);
    // 混乱赋值给它
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        mat.at(w,x,y) = complex(w+x, w-y) * 3.1415926;
    }
    Plaintext pt = Plaintext::from_cmat(mat, mods, delta);
    
    ComplexMatrixGroup mat2 = pt.to_cmat(delta);

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - mat.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_pt_encode_and_decode: max_diff="<<max_diff<<std::endl;
}