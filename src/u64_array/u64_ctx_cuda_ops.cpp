#include "u64_array.hpp"
#include "GPU/cuda_u64_ctx_ops.hpp"
void U64Context::iw_ntt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const
{
    CudaBuffer buf(size_ * sizeof(uint64_t));
    // 先做i-NTT
    cuda_i_ntt(dst, src, pnn_, I_mont_, mm_);
    // 然后做转置(2,p-1,n,n) -> (n,n,2,p-1)
    cuda_transpose_rect_restrict(buf, dst, 2*(p_-1), nn_);
    // 然后做W-NTT
    ntter_w->ntt_batch_cuda(buf, 2*nn_);
    // 然后做转置(n,n,2,p-1) -> (2,p-1,n,n)
    cuda_transpose_rect_restrict(dst, buf, nn_, 2*(p_-1));

}
void U64Context::iw_intt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const
{
    CudaBuffer buf(size_ * sizeof(uint64_t));
    // 先做转置
    cuda_transpose_rect_restrict(buf, src, 2*(p_-1), nn_);
    // 然后做W-iNTT
    ntter_w->intt_batch_cuda(buf, 2*nn_);
    // 然后做转置(n,n,2,p-1) -> (2,p-1,n,n)
    cuda_transpose_rect_restrict(dst, buf, nn_, 2*(p_-1));
    // 最后做i-intt
    cuda_i_intt(dst, dst, pnn_, I_inv_mont_, mm_);
}
void U64Context::xy_ntt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const
{
    CudaBuffer buf(size_ * sizeof(uint64_t));
    const size_t n = this->n_;
    const size_t nn = this->nn_;
    const size_t pnn = this->pnn_;
    const size_t pn = n*(p_-1);
    // 先转置
    // TODO: 更好的转置方法？
    for(int i=0; i<2*p_-2; i++)
    {
        cuda_transpose_rect_restrict(
            buf.slice(i*nn_*sizeof(uint64_t), (i+1)*nn_*sizeof(uint64_t)), 
            src.slice(i*nn_*sizeof(uint64_t), (i+1)*nn_*sizeof(uint64_t)), 
            n_, n_
        );
    }
    // NTT
    // 转置
    ntter_p->ntt_batch_cuda(buf.slice(0, pnn*sizeof(uint64_t)), pn);
    ntter_n->ntt_batch_cuda(buf.slice(pnn*sizeof(uint64_t), 2*pnn*sizeof(uint64_t)), pn);
    // 再转置
    for(int i=0; i<2*p_-2; i++)
    {
        cuda_transpose_rect_restrict(
            dst.slice(i*nn_*sizeof(uint64_t), (i+1)*nn_*sizeof(uint64_t)), 
            buf.slice(i*nn_*sizeof(uint64_t), (i+1)*nn_*sizeof(uint64_t)), 
            n_, n_
        );
    }
    ntter_n->ntt_batch_cuda(dst.slice(0, pnn*sizeof(uint64_t)), pn);
    ntter_p->ntt_batch_cuda(dst.slice(pnn*sizeof(uint64_t), 2*pnn*sizeof(uint64_t)), pn);
}
void U64Context::xy_intt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const
{
    throw std::runtime_error("TODO: 尚未实现\n");
}