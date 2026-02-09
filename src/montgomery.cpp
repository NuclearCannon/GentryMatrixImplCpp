

#include "montgomery.hpp"

MontgomeryMultiplier::MontgomeryMultiplier(uint64_t M)
{

    assert(M % 2 == 1);
    M_ = M;
    R_ = mod_pow(2, 64, M);
    Rinv_ = mod_inv(R_, M);
    R2_ = mod_mul(R_, R_, M);
    N1_ = - modinv64(M);

}


uint64_t MontgomeryMultiplier::pow(uint64_t base, size_t e) const
{
    u64 result = 1;
    while (e > 0) {
        if (e & 1) {
            result = mul(result, base);
        }
        base = mul(base, base);
        e >>= 1;
    }
    return result;
}

vec64 MontgomeryMultiplier::get_powers(uint64_t x, size_t len) const
{
    vec64 res(len);
    res[0] = encode(1);
    for(int i=1; i<len; i++)res[i] = mul(res[i-1], x);
    return res;
}

vec64 MontgomeryMultiplier::batch_encode(const vec64& src) const
{
    vec64 res(src);
    batch_encode_inplace(res);
    return res;
}

void MontgomeryMultiplier::batch_encode_inplace(vec64& v) const
{
    for(auto& i:v)i=encode(i);
}
void MontgomeryMultiplier::batch_decode_inplace(vec64& v) const
{
    for(auto& i:v)i=decode(i);
}


void MontgomeryMultiplier::vec_mul_mont(vec64& dst, const vec64& src1, const vec64& src2) const
{
    size_t size = dst.size();
    assert(src1.size() == size);
    assert(src2.size() == size);
    for(size_t i=0; i<size; i++)dst[i] = mul(src1[i], src2[i]);
}