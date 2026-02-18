#include "GentryPoly.hpp"
#include "math_utils.hpp"

void GPComponent::i_ntt(const GPCCtx& ctx)
{
    assert(like(ctx));
    size_t pnn = (p_-1)*n_*n_;
    assert(data_.size() == 2*pnn);
    ctx.ntter_i_.ntt_batch(data_.data(), data_.data()+pnn, pnn);
}
void GPComponent::i_intt(const GPCCtx& ctx)
{
    assert(like(ctx));
    size_t pnn = (p_-1)*n_*n_;
    assert(data_.size() == 2*pnn);
    ctx.ntter_i_.intt_batch(data_.data(), data_.data()+pnn, pnn);
}
void GPComponent::w_ntt(const GPCCtx& ctx)
{
    std::vector<uint64_t> buf(data_.size());
    size_t nn = n_*n_;
    transpose_rect_restrict(buf.data(), data_.data(), 2*p_-2, nn);
    ctx.ntter_w_.ntt_batch(buf.data(), 2*nn);
    transpose_rect_restrict(data_.data(), buf.data(), nn, 2*p_-2);
}
void GPComponent::w_intt(const GPCCtx& ctx)
{
    std::vector<uint64_t> buf(data_.size());
    size_t nn = n_*n_;
    transpose_rect_restrict(buf.data(), data_.data(), 2*p_-2, nn);
    ctx.ntter_w_.intt_batch(buf.data(), 2*nn);
    transpose_rect_restrict(data_.data(), buf.data(), nn, 2*p_-2);
}
void GPComponent::x_ntt(const GPCCtx& ctx)
{
    size_t n = n_;
    size_t nn = n*n;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    for(int i=0;i<2*pnn;i+=nn)transpose_inplace(data_.data()+i, n);
    ctx.ntter_p_.ntt_batch(data_.data(), pn);
    ctx.ntter_n_.ntt_batch(data_.data() + pnn, pn);
    for(int i=0;i<2*pnn;i+=nn)transpose_inplace(data_.data()+i, n);
}
void GPComponent::x_intt(const GPCCtx& ctx)
{
    size_t n = n_;
    size_t nn = n*n;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    

    for(int i=0;i<2*pnn;i+=nn)transpose_inplace(data_.data()+i, n);
    ctx.ntter_p_.intt_batch(data_.data(), pn);
    ctx.ntter_n_.intt_batch(data_.data() + pnn, pn);
    for(int i=0;i<2*pnn;i+=nn)transpose_inplace(data_.data()+i, n);
}
void GPComponent::y_ntt(const GPCCtx& ctx)
{
    size_t n = n_;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    
    ctx.ntter_n_.ntt_batch(data_.data(), pn);
    ctx.ntter_p_.ntt_batch(data_.data() + pnn, pn);
}
void GPComponent::y_intt(const GPCCtx& ctx)
{
    size_t n = n_;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    
    ctx.ntter_n_.intt_batch(data_.data(), pn);
    ctx.ntter_p_.intt_batch(data_.data() + pnn, pn);
}