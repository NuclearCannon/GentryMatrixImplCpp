#include "GentryPoly.hpp"
#include "math_utils.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"


void GPComponentCuda::i_ntt(const GPCCtx& ctx)
{
    assert(like(ctx));
    size_t pnn = (p_-1)*n_*n_;
    assert(data_.size() == 2*pnn*sizeof(uint64_t));
    ctx.ntter_i_->ntt_batch_cuda(data_, pnn);
}
void GPComponentCuda::i_intt(const GPCCtx& ctx)
{
    assert(like(ctx));
    size_t pnn = (p_-1)*n_*n_;
    assert(data_.size() == 2*pnn*sizeof(uint64_t));
    ctx.ntter_i_->intt_batch_cuda(data_, pnn);

}
void GPComponentCuda::w_ntt(const GPCCtx& ctx)
{
    CudaBuffer buf(data_.size());
    size_t nn = n_*n_;
    // 然后做转置(2,p-1,n,n) -> (n,n,2,p-1)
    cuda_transpose_rect_restrict(buf, data_, 2*(p_-1), nn);
    // 然后做W-NTT
    ctx.ntter_w_->ntt_batch_cuda(buf, 2*nn);
    // 然后做转置(n,n,2,p-1) -> (2,p-1,n,n)
    cuda_transpose_rect_restrict(data_, buf, nn, 2*(p_-1));
}
void GPComponentCuda::w_intt(const GPCCtx& ctx)
{
    CudaBuffer buf(data_.size());
    size_t nn = n_*n_;
    // 然后做转置(2,p-1,n,n) -> (n,n,2,p-1)
    cuda_transpose_rect_restrict(buf, data_, 2*(p_-1), nn);
    // 然后做W-NTT
    ctx.ntter_w_->intt_batch_cuda(buf, 2*nn);
    // 然后做转置(n,n,2,p-1) -> (2,p-1,n,n)
    cuda_transpose_rect_restrict(data_, buf, nn, 2*(p_-1));
}
void GPComponentCuda::x_ntt(const GPCCtx& ctx)
{
    CudaBuffer buf(data_.size());
    size_t n = n_;
    size_t nn = n*n;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    for(int i=0; i<2*p_-2; i++)
    {
        cuda_transpose_rect_restrict(
            buf.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            data_.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            n_, n_
        );
    }
    // NTT
    // 转置
    ctx.ntter_p_->ntt_batch_cuda(buf.slice(0, pnn*sizeof(uint64_t)), pn);
    ctx.ntter_n_->ntt_batch_cuda(buf.slice(pnn*sizeof(uint64_t), 2*pnn*sizeof(uint64_t)), pn);
    // 再转置
    for(int i=0; i<2*p_-2; i++)
    {
        cuda_transpose_rect_restrict(
            data_.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            buf.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            n_, n_
        );
    }
}
void GPComponentCuda::x_intt(const GPCCtx& ctx)
{
    CudaBuffer buf(data_.size());
    size_t n = n_;
    size_t nn = n*n;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    for(int i=0; i<2*p_-2; i++)
    {
        cuda_transpose_rect_restrict(
            buf.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            data_.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            n_, n_
        );
    }
    // NTT
    // 转置
    ctx.ntter_p_->intt_batch_cuda(buf.slice(0, pnn*sizeof(uint64_t)), pn);
    ctx.ntter_n_->intt_batch_cuda(buf.slice(pnn*sizeof(uint64_t), 2*pnn*sizeof(uint64_t)), pn);
    // 再转置
    for(int i=0; i<2*p_-2; i++)
    {
        cuda_transpose_rect_restrict(
            data_.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            buf.slice(i*nn*sizeof(uint64_t), (i+1)*nn*sizeof(uint64_t)), 
            n_, n_
        );
    }
}
void GPComponentCuda::y_ntt(const GPCCtx& ctx)
{
    size_t n = n_;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    ctx.ntter_n_->ntt_batch_cuda(data_.slice(0, pnn*sizeof(uint64_t)), pn);
    ctx.ntter_p_->ntt_batch_cuda(data_.slice(pnn*sizeof(uint64_t), 2*pnn*sizeof(uint64_t)), pn);
}
void GPComponentCuda::y_intt(const GPCCtx& ctx)
{
    size_t n = n_;
    size_t pn = (p_-1)*n;
    size_t pnn = pn*n;
    ctx.ntter_n_->intt_batch_cuda(data_.slice(0, pnn*sizeof(uint64_t)), pn);
    ctx.ntter_p_->intt_batch_cuda(data_.slice(pnn*sizeof(uint64_t), 2*pnn*sizeof(uint64_t)), pn);
}