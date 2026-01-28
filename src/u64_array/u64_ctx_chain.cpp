#include "u64_array.hpp"
#include <cassert>


U64CtxChain::U64CtxChain(int n, int p, const vec64& mods, const vec64& roots):
    n_(n), p_(p), mods_(mods), roots_(roots), mod_prod_(1)
{
    chain_len_ = mods.size();
    assert(chain_len_ == roots.size());
    for(size_t i=0; i<chain_len_; i++)
    {
        ctxs_.push_back(
            std::make_shared<U64Context>(
                n, p, mods_[i], roots_[i]
            )
        );
        fmpz_mul_ui(mod_prod_.raw(), mod_prod_.raw(), mods_[i]);
    }
    
}
U64CtxChain::~U64CtxChain()
{
    // do nothing
}

void U64CtxChain::iw_ntt(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->iw_ntt(dst[i], src[i]);
    }
}
void U64CtxChain::iw_intt(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->iw_intt(dst[i], src[i]);
    }
}
void U64CtxChain::xy_ntt(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->xy_ntt(dst[i], src[i]);
    }
}
void U64CtxChain::xy_intt(vv64& dst, const vv64& src) const
{
    assert(dst.size() == chain_len_);
    assert(src.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->xy_intt(dst[i], src[i]);
    }
}

// 逐位加法
void U64CtxChain::add(vv64& dst, const vv64& src1, const vv64& src2) const
{
    assert(dst.size() == chain_len_);
    assert(src1.size() == chain_len_);
    assert(src2.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->add(dst[i], src1[i], src2[i]);
    }
}
// 逐位减法
void U64CtxChain::sub(vv64& dst, const vv64& src1, const vv64& src2) const
{
    assert(dst.size() == chain_len_);
    assert(src1.size() == chain_len_);
    assert(src2.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->sub(dst[i], src1[i], src2[i]);
    }
}
// 逐位乘法
void U64CtxChain::mul(vv64& dst, const vv64& src1, const vv64& src2) const
{
    assert(dst.size() == chain_len_);
    assert(src1.size() == chain_len_);
    assert(src2.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->mul(dst[i], src1[i], src2[i]);
    }
}
// 逐位负
void U64CtxChain::neg(vv64& dst, const vv64& src1) const
{
    assert(dst.size() == chain_len_);
    assert(src1.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->neg(dst[i], src1[i]);
    }
}
// 标量乘
void U64CtxChain::mul_scalar(vv64& dst, const vv64& src_vec, u64 src_scalar) const
{
    assert(dst.size() == chain_len_);
    assert(src_vec.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        ctxs_[i]->mul_scalar(dst[i], src_vec[i], src_scalar);
    }
}
// 比较
bool U64CtxChain::eq(const vv64& src1, const vv64& src2) const
{
    assert(src1.size() == chain_len_);
    assert(src2.size() == chain_len_);
    for(size_t i=0;i<chain_len_;i++)
    {
        if (!(ctxs_[i]->eq(src1[i], src2[i])))return false;
    }
    return true;
}