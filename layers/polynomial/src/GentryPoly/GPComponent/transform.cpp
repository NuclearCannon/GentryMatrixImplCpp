#include "GentryPoly.hpp"
#include <cstring>
#include "modops.hpp"

// void transpose(GPComponent& dst);
// void conj(GPComponent& dst);
// void w_inv(GPComponent& dst);

void GPComponent::transpose(GPComponent& dst_) const
{
    const vec64& src = this->data_;
    vec64& dst = dst_.data_;
    assert(&src != &dst);
    const size_t n = n_;
    const size_t nn = n*n;
    const size_t pnn = (p_-1)*nn;

    for(int i=0; i<2; i++)
    {
        for(int w=0; w<p_-1; w++)
        {
            for(int x=0; x<n; x++)
            {
                for(int y=0; y<n; y++)
                {
                    dst[i*pnn + w*nn + x*n + y] = src[i*pnn + w*nn + y*n + x];
                }
            }
        }
    }
}
void GPComponent::conj(GPComponent& dst_) const
{
    const vec64& src = this->data_;
    vec64& dst = dst_.data_;
    assert(&src != &dst);
    const size_t n = n_;
    const size_t nn = n*n;
    const size_t pnn = (p_-1)*nn;

    for(int w=0; w<p_-1; w++)
    {
        for(int x=0; x<n_; x++)
        {
            for(int y=0; y<n_; y++)
            {
                // i==0
                dst[w*nn + x*n_ + y] = src[w*nn + x*n_ + y];
                // i==1
                dst[pnn + w*nn + x*n_ + y] = (q_ - src[pnn + w*nn + x*n_ + y]) % q_;
            }
        }
    }
    
}
void GPComponent::w_inv(GPComponent& dst_) const
{
    const vec64& src = this->data_;
    vec64& dst = dst_.data_;
    assert(&src != &dst);
    const size_t n = n_;
    const size_t nn = n*n;
    const size_t pnn = (p_-1)*nn;

    memset(dst.data(), 0, dst.size()*sizeof(uint64_t));

    vec64 gpp = get_powers(3,p_-1,p_), gpp_backward(p_);
    for(int i=0; i<p_-1; i++)gpp_backward[gpp[i]] = i;
    // gpp_backward[0]是无意义的

    for(int w=0; w<p_-1; w++)
    {
        size_t w2 = gpp_backward[p_-gpp[w]];
        for(int i=0; i<2; i++)
        {
            for(int xy=0; xy<nn; xy++)
            {        
                // 迁移src[:,w,:,:]到dst[:,(p-w)%p,:,:]
                // src[:,w,:,:] 是关于 W^{gamma^w}的系数。现在需要变成W^{-{gamma^w}} = W^{p-{gamma^w}}
                const uint64_t& s = src[i*pnn + w*nn + xy];
                uint64_t& d = dst[i*pnn + w2*nn + xy];
                d = mod_add(d, s, q_);
            }
        }
    }
}