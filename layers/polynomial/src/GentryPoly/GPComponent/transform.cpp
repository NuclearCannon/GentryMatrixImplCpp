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

void GPComponent::automorphism_XY(GPComponent& dst_, int x, int y) const
{
    assert(x > 0);
    assert(y > 0);
    assert(x % 2 == 1);
    assert(y % 2 == 1);
    const vec64& src = this->data_;
    vec64& dst = dst_.data_;
    assert(&src != &dst);
    const size_t n = n_;
    const size_t nn = n*n;
    const size_t pnn = (p_-1)*nn;

    memset(dst.data(), 0, dst.size()*sizeof(uint64_t));

    for(int w=0; w<p_-1; w++)
    {
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                // 将src[w, i, j]移动到dst[w, i2, j2]
                int i2 = (i * x) % (4*n);   // zeta^4n=1，因此可以取模4n
                int j2 = (j * y) % (4*n);
                int iexp = 0;   // 虚数单位i的指数
                if(i2>=n){iexp++; i2-=n;}
                if(i2>=n){iexp++; i2-=n;}
                if(i2>=n){iexp++; i2-=n;}

                if(j2>=n){iexp--; j2-=n;}
                if(j2>=n){iexp--; j2-=n;}
                if(j2>=n){iexp--; j2-=n;}
                // iexp in [-3,3]
                iexp = (iexp + 4) % 4;
                // iexp in [0,3]


                uint64_t s1 = src[w*nn + i*n + j];
                uint64_t s2 = src[pnn + w*nn + i*n + j];
                uint64_t& d1 = dst[w*nn + i2*n + j2];
                uint64_t& d2 = dst[pnn + w*nn + i2*n + j2];
                if(iexp == 0)
                {
                    d1 = mod_add(d1, s1, q_);
                    d2 = mod_add(d2, s2, q_);
                }
                else if (iexp == 1)
                {
                    d1 = mod_sub(d1, s2, q_);
                    d2 = mod_add(d2, s1, q_);
                }
                else if (iexp == 2)
                {
                    d1 = mod_sub(d1, s1, q_);
                    d2 = mod_sub(d2, s2, q_);
                }
                else if (iexp == 3)
                {
                    d1 = mod_add(d1, s2, q_);
                    d2 = mod_sub(d2, s1, q_);
                }

            }
        }
    }
}

void GPComponent::automorphism_W(GPComponent& dst_, int delta_w) const
{

    const vec64& src = this->data_;
    vec64& dst = dst_.data_;
    assert(&src != &dst);
    const size_t n = n_;
    const size_t nn = n*n;
    const size_t pnn = (p_-1)*nn;

    memset(dst.data(), 0, dst.size()*sizeof(uint64_t));

    for(int w=0; w<p_-1; w++)
    {
        int w2 = (w-delta_w+p_-1)%(p_-1);
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                


                uint64_t s1 = src[w*nn + i*n + j];
                uint64_t s2 = src[pnn + w*nn + i*n + j];
                uint64_t& d1 = dst[w2*nn + i*n + j];
                uint64_t& d2 = dst[pnn + w2*nn + i*n + j];
                d1 = s1;
                d2 = s2;

            }
        }
    }
}