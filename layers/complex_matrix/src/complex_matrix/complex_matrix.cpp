#include "complecx_matrix.hpp"

ComplexMatrixGroup::ComplexMatrixGroup(size_t n, size_t p, std::vector<complex> data):
    p_(p), n_(n), data_(data)
{

}

ComplexMatrixGroup::ComplexMatrixGroup(size_t n, size_t p):
    p_(p), n_(n), data_((p-1)*n*n, complex(0))
{

}