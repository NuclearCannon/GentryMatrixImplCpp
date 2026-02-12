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

    for(int w=0; w<p_-1; w++)
    {
        for(int x=0; x<n_; x++)
        {
            for(int y=0; y<n_; y++)
            {
                // i==0
                dst[w*nn_ + x*n_ + y] = src[w*nn_ + x*n_ + y];
                // i==1
                dst[pnn_ + w*nn_ + x*n_ + y] = (q_ - src[pnn_ + w*nn_ + x*n_ + y]) % q_;
            }
        }
    }
    
}
void U64Context::w_inv(vec64& dst, const vec64& src) const
{
    assert(&src != &dst);
    memset(dst.data(), 0, dst.size()*sizeof(uint64_t));

    vec64 gpp = get_powers(3,p_-1,p_), gpp_backward(p_);
    for(int i=0; i<p_-1; i++)gpp_backward[gpp[i]] = i;
    // gpp_backward[0]是无意义的
    size_t nn = nn_;

    for(int w=0; w<p_-1; w++)
    {
        size_t w2 = gpp_backward[p_-gpp[w]];
        for(int i=0; i<2; i++)
        {
            for(int xy=0; xy<nn; xy++)
            {        
                // 迁移src[:,w,:,:]到dst[:,(p-w)%p,:,:]
                // src[:,w,:,:] 是关于 W^{gamma^w}的系数。现在需要变成W^{-{gamma^w}} = W^{p-{gamma^w}}
                const uint64_t& s = src[i*pnn_ + w*nn + xy];
                uint64_t& d = dst[i*pnn_ + w2*nn + xy];
                d = mod_add(d, s, q_);
                    
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
