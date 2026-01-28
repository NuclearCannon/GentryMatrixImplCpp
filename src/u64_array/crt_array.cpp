#include "u64_array.hpp"
#include "CRT.hpp"
#include <cassert>
// 构造函数: 全零数组
CRTArray::CRTArray(std::shared_ptr<const U64CtxChain> cc):
    cc_(cc),
    data_(cc_->get_chain_length(), std::vector<uint64_t>(cc_->get_size(), 0))
{
    
}

// 构造函数: 从 vector 构造（需正确大小）
CRTArray::CRTArray(vv64 data, std::shared_ptr<const U64CtxChain> cc):
    cc_(cc),
    data_(std::move(data))
{
    assert(data_.size() == cc_->get_chain_length());
    auto size = cc_->get_size();
    for(auto& t : data_)
    {
        assert(t.size() == size);
    }
}

CRTArray CRTArray::from_fmpz_vector(const fmpz_vector& data, std::shared_ptr<const U64CtxChain> cc)
{
    assert(data.len() == cc->get_size());
    auto data2 = crt(data, cc->get_mods());
    return CRTArray(data2, cc);
}
fmpz_vector CRTArray::to_fmpz_vector() const
{
    int size = cc_->get_size();
    fmpz_vector dst(size);
    icrt(dst, data_, cc_->get_mods());
    return dst;
}

fmpz_vector CRTArray::to_fmpz_vector_centered() const
{
    int size = cc_->get_size();
    fmpz_vector buf1(size);
    icrt(buf1, data_, cc_->get_mods());
    fmpz_vector buf2 = buf1.mod_centered(cc_->get_mod_prod());
    return buf2;
}

CRTArray::~CRTArray() = default;
CRTArray::CRTArray(const CRTArray&) = default;
CRTArray::CRTArray(CRTArray&&) = default;


// 运算符重载：逐元素操作
CRTArray CRTArray::add(const CRTArray& other) const
{
    CRTArray res(cc_);
    cc_->add(res.data_, data_, other.data_);
    return res;
}
CRTArray CRTArray::neg() const
{
    CRTArray res(cc_);
    cc_->neg(res.data_, data_);
    return res;
}
CRTArray CRTArray::sub(const CRTArray& other) const
{
    CRTArray res(cc_);
    cc_->sub(res.data_, data_, other.data_);
    return res;
}
CRTArray CRTArray::mul(const CRTArray& other) const
{
    CRTArray res(cc_);
    cc_->mul(res.data_, data_, other.data_);
    return res;
}
CRTArray CRTArray::mul_scalar(u64 other) const
{
    CRTArray res(cc_);
    cc_->mul_scalar(res.data_, data_, other);
    return res;
}
CRTArray CRTArray::mul_poly(const CRTArray& other) const
{
    vv64 ntt1(cc_->get_chain_length(), std::vector<uint64_t>(cc_->get_size()));
    vv64 ntt2(cc_->get_chain_length(), std::vector<uint64_t>(cc_->get_size()));
    vv64 ntt3(cc_->get_chain_length(), std::vector<uint64_t>(cc_->get_size()));
    cc_->iw_ntt(ntt1, data_);
    cc_->xy_ntt(ntt1, ntt1);
    cc_->iw_ntt(ntt2, other.data_);
    cc_->xy_ntt(ntt2, ntt2);
    cc_->mul(ntt3, ntt1, ntt2);
    cc_->xy_intt(ntt3, ntt3);
    cc_->iw_intt(ntt3, ntt3);
    return CRTArray(std::move(ntt3), cc_);
}

bool CRTArray::eq(const CRTArray& other) const
{
    return cc_->eq(data_, other.data_);
}

CRTArray CRTArray::iw_ntt() const
{
    CRTArray res(cc_);
    cc_->iw_ntt(res.data_, data_);
    return res;
}
CRTArray CRTArray::iw_intt() const
{
    CRTArray res(cc_);
    cc_->iw_intt(res.data_, data_);
    return res;
}
CRTArray CRTArray::xy_ntt() const
{
    CRTArray res(cc_);
    cc_->xy_ntt(res.data_, data_);
    return res;
}
CRTArray CRTArray::xy_intt() const
{
    CRTArray res(cc_);
    cc_->xy_intt(res.data_, data_);
    return res;
}
CRTArray CRTArray::all_ntt() const
{
    return iw_ntt().xy_ntt();
}
CRTArray CRTArray::all_intt() const
{
    return xy_intt().iw_intt();
}
