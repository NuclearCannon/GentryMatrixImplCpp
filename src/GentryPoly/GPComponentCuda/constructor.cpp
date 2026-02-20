#include "GentryPoly.hpp"

GPComponentCuda::GPComponentCuda(size_t n, size_t p, uint64_t q):
    n_(n),
    p_(p),
    q_(q),
    mm_(q),
    data_(2*n*n*(p-1)*sizeof(uint64_t))
{

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


GPComponentCuda GPComponentCuda::copy_from_cpu(const GPComponent& other)
{
    GPComponentCuda result(other.n_, other.p_, other.q_);
    result.set_from_cpu(other);
    return result;
}