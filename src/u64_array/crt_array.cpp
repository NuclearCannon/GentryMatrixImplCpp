#include "u64_array.hpp"
#include "CRT.hpp"
#include <cassert>
#include <cstring>






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

CRTArray::CRTArray(const vec64& data, std::shared_ptr<const U64CtxChain> cc):
    cc_(cc),
    data_(cc_->get_chain_length(), data)
{
    int size = cc_->get_size();
    for(int i=0; i<cc_->get_chain_length(); i++)
    {
        u64 mod = cc_->get_mods()[i];
        vec64& row = data_[i];
        for(int j=0; j<size; j++)row[j] %= mod;
    }
}

void CRTArray::set_from_raw(const vec64& data)
{
    int size = cc_->get_size();
    for(int i=0; i<cc_->get_chain_length(); i++)
    {
        u64 mod = cc_->get_mods()[i];
        vec64& row = data_[i];
        for(int j=0; j<size; j++)row[j] = data[j] % mod;
    }
}

void CRTArray::set_to_zero(const vec64& data)
{
    int size = cc_->get_size();
    for(int i=0; i<cc_->get_chain_length(); i++)
    {
        vec64& row = data_[i];
        memset(row.data(), 0, row.size()*sizeof(u64));
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
void CRTArray::adde(const CRTArray& other)
{
    cc_->add(data_, data_, other.data_);
}
void CRTArray::mul_scalar_e(u64 other)
{
    cc_->mul_scalar(data_, data_, other);
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
CRTArray CRTArray::mul_mont(const CRTArray& other) const
{
    CRTArray res(cc_);
    cc_->mul_mont(res.data_, data_, other.data_);
    return res;
}
void CRTArray::mul_mont3(CRTArray& dst, const CRTArray& src1, const CRTArray& src2)
{
    assert(dst.get_cc().get() == src1.get_cc().get());
    assert(dst.get_cc().get() == src2.get_cc().get());
    dst.get_cc()->mul_mont(dst.data_, src1.data_, src2.data_);
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
    CRTArray res(cc_);
    cc_->iw_ntt(res.data_, data_);
    cc_->xy_ntt(res.data_, res.data_);
    return res;
}
void CRTArray::all_ntt_to(CRTArray& dst) const
{
    cc_->iw_ntt(dst.data_, data_);
    cc_->xy_ntt(dst.data_, dst.data_);
}
CRTArray CRTArray::all_intt() const
{
    CRTArray res(cc_);
    cc_->xy_intt(res.data_, data_);
    cc_->iw_intt(res.data_, res.data_);
    return res;
}
void CRTArray::all_intt_to(CRTArray& dst) const
{
    cc_->xy_intt(dst.data_, data_);
    cc_->iw_intt(dst.data_, dst.data_);
}
CRTArray CRTArray::mont_encode() const
{
    CRTArray res(cc_);
    cc_->mont_encode(res.data_, data_);
    return res;
}
CRTArray CRTArray::mont_decode() const
{
    CRTArray res(cc_);
    cc_->mont_decode(res.data_, data_);
    return res;
}
void CRTArray::mont_encode_inplace()
{
    cc_->mont_encode(data_, data_);
}
void CRTArray::mont_decode_inplace()
{
    cc_->mont_decode(data_, data_);
}
