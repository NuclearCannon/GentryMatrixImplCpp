#include "complecx_matrix.hpp"
#include "random.hpp"


ComplexMatrixGroup ComplexMatrixGroup::random(double abs_max, size_t n, size_t p)
{
    size_t len = (p-1)*n*n;
    std::vector<complex> result(len);
    for(size_t i=0; i<len; i++)
    {
        double abs = ramdom_generators::random_real() * abs_max;
        double angle = ramdom_generators::random_real() * M_PI * 2;
        result[i] = std::polar(abs, angle);
    }
    return ComplexMatrixGroup(n, p, result);

}