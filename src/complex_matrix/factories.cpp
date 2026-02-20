#include "complecx_matrix.hpp"
#include "random.hpp"


ComplexMatrixGroup ComplexMatrixGroup::random(double abs_max, size_t n, size_t p)
{
    size_t len = (p-1)*n*n;
    std::vector<complex> result(len);
    for(size_t i=0; i<len; i++)
    {
        double abs = random_generators::random_real() * abs_max;
        double angle = random_generators::random_real() * M_PI * 2;
        result[i] = std::polar(abs, angle);
    }
    return ComplexMatrixGroup(n, p, result);

}


ComplexMatrixGroup ComplexMatrixGroup::eye(size_t n, size_t p)
{
    size_t len = (p-1)*n*n;
    std::vector<complex> result(len);
    for(size_t i=0; i<n; i++)
    {
        for(size_t w=0; w<p-1; w++)
        {
            result[w*n*n + i*(n+1)] = 1;
        }
    }
    return ComplexMatrixGroup(n, p, result);

}

ComplexMatrixGroup ComplexMatrixGroup::example_matrix(size_t n, size_t p)
{
    ComplexMatrixGroup result(n, p);
    for(size_t w=0; w<p-1; w++)
    {
        for(size_t i=0; i<n; i++)
        {
            for(size_t j=0; j<n; j++)
            {
                result.at(w, i, j) = complex(i,j);
            }
        }
    }

    return result;

}