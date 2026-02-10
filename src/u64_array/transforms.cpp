#include "u64_array.hpp"
#include <cstring>


void U64Context::transpose(vec64& dst, const vec64& src) const
{
    assert(&src != &dst);
    for(int i=0; i<2; i++)
    {
        for(int w=0; w<p_-1; w++)
        {
            for(int x=0; x<n_; x++)
            {
                for(int y=0; y<n_; y++)
                {
                    dst[i*pnn_ + w*nn_ + x*n_ + y] = src[i*pnn_ + w*nn_ + y*n_ + x];
                }
            }
        }
    }
}
void U64Context::conj(vec64& dst, const vec64& src) const
{
    assert(&src != &dst);
    for(int i=0; i<2; i++)
    {
        for(int w=0; w<p_-1; w++)
        {
            for(int x=0; x<n_; x++)
            {
                for(int y=0; y<n_; y++)
                {
                    dst[i*pnn_ + w*nn_ + x*n_ + y] = src[(1-i)*pnn_ + w*nn_ + x*n_ + y];
                }
            }
        }
    }
}
void U64Context::w_inv(vec64& dst, const vec64& src) const
{
    assert(&src != &dst);
    memset(dst.data(), 0, dst.size()*sizeof(u64));
    for(int i=0; i<2; i++)
    {
        for(int x=0; x<n_; x++)
        {
            for(int y=0; y<n_; y++)
            {
                for(int w=0; w<p_-1; w++)
                {
                    // 迁移src[:,w,:,:]到dst[:,(p-w)%p,:,:]
                    const u64& s = src[i*pnn_ + w*nn_ + x*n_ + y];

                    int w2 = (p_-w)%p_;
                    if (w2 == p_-1)
                    {
                        for(int w3 = 0; w3<p_-1; w3++)
                        {
                            u64& d = dst[i*pnn_ + w3*nn_ + x*n_ + y];
                            d = mod_sub(d, s, q_);
                        }
                    }
                    else
                    {
                        u64& d = dst[i*pnn_ + w2*nn_ + x*n_ + y];
                        d = mod_add(d, s, q_);
                    }
                }
            }
        }
        
    }
}


void U64CtxChain::transpose(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->transpose(dst[i], src[i]);
    }
}

void U64CtxChain::conj(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->conj(dst[i], src[i]);
    }
}

void U64CtxChain::w_inv(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->w_inv(dst[i], src[i]);
    }
}


CRTArray CRTArray::transpose() const{
    CRTArray res(cc_);
    cc_->transpose(res.data_, data_);
    return res;
}

CRTArray CRTArray::conj() const{
    CRTArray res(cc_);
    cc_->conj(res.data_, data_);
    return res;
}

CRTArray CRTArray::w_inv() const{
    CRTArray res(cc_);
    cc_->w_inv(res.data_, data_);
    return res;
}
