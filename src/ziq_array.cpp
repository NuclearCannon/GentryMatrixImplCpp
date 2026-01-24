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

ZiqArray::ZiqArray(const ZiqArray&) = default;
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

ZiqArray ZiqArray::mul_poly(const ZiqArray& other) const
{
    ZiqArray ntt1 = this->iw_ntt().xy_ntt();
    ZiqArray ntt2 = other.iw_ntt().xy_ntt();
    ZiqArray ntt3 = ntt1.mul(ntt2);
    ZiqArray poly3 = ntt3.xy_intt().iw_intt();
    return poly3;
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



ZiqArray ZiqArray::transpose() const
{
    int p1 = ctx_->get_p()-1, n = ctx_->get_n(), size = ctx_->get_size();
    int nn = n*n;
    int pnn = p1*nn;
    
    fmpz_vector result(size);
    for(int i=0;i<2;i++)
    {
        for(int w=0;w<p1;w++)
        {
            for(int x=0;x<n;x++)
            {
                for(int y=0;y<n;y++)
                {
                    fmpz_set(
                        result[i*pnn + w*nn + x*n + y],
                        data_[i*pnn + w*nn + y*n + x]
                    );
                }
            }
        }
    }
    return ZiqArray(std::move(result), ctx_);
}
ZiqArray ZiqArray::conj() const
{
    // TODO: 可能都有很大的优化空间
    int p1 = ctx_->get_p()-1, n = ctx_->get_n(), size = ctx_->get_size();
    int nn = n*n;
    int pnn = p1*nn;
    
    fmpz_vector result(size);
    for(int i=0;i<2;i++)
    {
        for(int w=0;w<p1;w++)
        {
            for(int x=0;x<n;x++)
            {
                for(int y=0;y<n;y++)
                {
                    if (i == 0)
                    {
                        fmpz_set(
                            result[i*pnn + w*nn + x*n + y],
                            data_[i*pnn + w*nn + x*n + y]
                        );
                    }
                    else
                    {
                        fmpz_neg(
                            result[i*pnn + w*nn + x*n + y],
                            data_[i*pnn + w*nn + x*n + y]
                        );
                    }
                    
                }
            }
        }
    }
    return ZiqArray(std::move(result), ctx_);
}
ZiqArray ZiqArray::w_inversion() const
{
    // TODO: 可能都有很大的优化空间
    int p = ctx_->get_p();
    int p1 = p-1, n = ctx_->get_n(), size = ctx_->get_size();
    int nn = n*n;
    int pnn = p1*nn;
    
    fmpz_vector result = fmpz_vector::zeros(size);
    for(int i=0;i<2;i++)
    {
        for(int w=0;w<p1;w++)
        {
            for(int x=0;x<n;x++)
            {
                for(int y=0;y<n;y++)
                {
                    const fmpz* src = data_[i*pnn + w*nn + x*n + y];
                    // 对于data[i,w,x,y]
                    // 对应的W次数为w
                    // 新的w次数为(p-w)%p
                    // p-w有可能为p-1，此时需要特殊考虑
                    int w2 = (p-w)%p;
                    if (w2 == p1)
                    {
                        // 特殊考虑
                        for(int k=0;k<p1;k++)
                        {
                            fmpz* dst = result[i*pnn + k*nn + x*n + y];
                            fmpz_sub(dst, dst, src);
                        }
                    }
                    else
                    {
                        fmpz* dst = result[i*pnn + w2*nn + x*n + y];
                        fmpz_add(dst, dst, src);
                    }

                }
            }
        }
    }
    return ZiqArray(std::move(result), ctx_);
}