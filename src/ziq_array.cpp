#include "ziq_array.hpp"



ZiqArray::ZiqArray(int n, int p, const ZiqArrayContext* ctx):
    ctx_(ctx), 
    size_(2*n*n*(p-1)),
    data_(fmpz_vector::zeros(n))
{

}

// 构造函数: 从 vector 构造（需正确大小）
ZiqArray::ZiqArray(fmpz_vector data, const ZiqArrayContext* ctx):
    ctx_(ctx), 
    size_(data.len()),
    data_(std::move(data))
{
    _fmpz_vec_scalar_mod_fmpz(data_.raw(), data_.raw(), size_, ctx->q());
}

ZiqArray::~ZiqArray() = default;


ZiqArray::ZiqArray(ZiqArray&&) = default;

// 运算符重载：逐元素操作
ZiqArray ZiqArray::add(const ZiqArray& other) const
{
    assert(ctx_ == other.ctx_);
    fmpz_vector result(size_);
    ctx_->add(result, data_, other.data_);
    return ZiqArray(result, ctx_);
}

ZiqArray ZiqArray::neg() const
{
    fmpz_vector result(size_);
    ctx_->neg(result, data_);
    return ZiqArray(result, ctx_);
}

ZiqArray ZiqArray::sub(const ZiqArray& other) const
{
    assert(ctx_ == other.ctx_);
    fmpz_vector result(size_);
    ctx_->sub(result, data_, other.data_);
    return ZiqArray(result, ctx_);
}

ZiqArray ZiqArray::mul(const ZiqArray& other) const
{
    assert(ctx_ == other.ctx_);
    fmpz_vector result(size_);
    ctx_->mul(result, data_, other.data_);
    return ZiqArray(result, ctx_);
}

ZiqArray ZiqArray::mul_scalar(const fmpz_t other) const
{
    fmpz_vector result(size_);
    ctx_->mul_scalar(result, data_, other);
    return ZiqArray(result, ctx_);
}

bool ZiqArray::eq(const ZiqArray& other) const
{
    assert(ctx_ == other.ctx_);
    return ctx_->eq(data_, other.data_);
}


ZiqArray ZiqArray::iw_ntt() const
{
    return ZiqArray(ctx_->iw_ntt(data_), ctx_);
}
ZiqArray ZiqArray::iw_intt() const
{
    return ZiqArray(ctx_->iw_intt(data_), ctx_);
}
ZiqArray ZiqArray::xy_ntt () const
{
    return ZiqArray(ctx_->xy_ntt(data_), ctx_);
}
ZiqArray ZiqArray::xy_intt() const
{
    return ZiqArray(ctx_->xy_intt(data_), ctx_);
}

ZiqArray ZiqArray::circledast(const ZiqArray& other) const
{
    assert(ctx_ == other.ctx_);
    return ZiqArray(ctx_->circledast(data_, other.data_), ctx_);
}

ZiqArray ZiqArray::ctx_switch(const ZiqArrayContext* new_ctx) const
{
    fmpz_vector data = data_.mod_centered(ctx_->q());
    return ZiqArray(data, new_ctx);
}