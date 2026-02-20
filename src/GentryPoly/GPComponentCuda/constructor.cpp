#include "GentryPoly.hpp"

GPComponentCuda::GPComponentCuda(size_t n, size_t p, uint64_t q):
    n_(n),
    p_(p),
    q_(q),
    mm_(q),
    data_(2*n*n*(p-1)*sizeof(uint64_t))
{

}

GPComponentCuda GPComponentCuda::from_buffer(size_t n, size_t p, uint64_t q, const uint64_t* data)
{
    GPComponentCuda result(n, p, q);
    result.set_from_buffer(data);
    return result;
}


GPComponentCuda::~GPComponentCuda() = default;
GPComponentCuda::GPComponentCuda(const GPComponentCuda& other):
    n_(other.n_),
    p_(other.p_),
    q_(other.q_),
    mm_(other.mm_),
    data_(other.data_.size())
{
    data_.copy_from_other(other.data_);
}
GPComponentCuda::GPComponentCuda(GPComponentCuda&&) = default;