#include "ntt.hpp"
#include "modops.hpp"


NTTerI::NTTerI(uint64_t q, uint64_t qroot):
    mm_(q)
{
    uint64_t I = mod_pow(qroot, (q-1)/4, q);
    uint64_t Iinv = mod_pow(I, 3, q);
    I_mont_ = mm_.encode(I);
    
    inv2_mont_ = mm_.encode(mod_inv(2, q));

    I2_inv_mont_ = mm_.mul(mm_.encode(Iinv), inv2_mont_);
}
void NTTerI::ntt_batch(
    uint64_t* low, 
    uint64_t* hig, 
    size_t batch_size
) const
{
    for(size_t i=0; i<batch_size; i++)
    {
        uint64_t R = low[i], I = mm_.mul(hig[i], I_mont_);
        uint64_t P = mod_add(R, I, mm_.M);
        uint64_t N = mod_sub(R, I, mm_.M);
        low[i] = P;
        hig[i] = N;
    }
}
void NTTerI::intt_batch(
    uint64_t* low, 
    uint64_t* hig, 
    size_t batch_size
) const
{
    for(size_t i=0; i<batch_size; i++)
    {
        uint64_t P = low[i], N = hig[i];
        uint64_t R = mm_.mul(mod_add(P, N, mm_.M), inv2_mont_);
        uint64_t I = mm_.mul(mod_sub(P, N, mm_.M), I2_inv_mont_);
        low[i] = R;
        hig[i] = I;
    }
}