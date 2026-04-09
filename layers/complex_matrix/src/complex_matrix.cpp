#include "complex_matrix.hpp"

ComplexMatrixGroup::ComplexMatrixGroup(size_t n, size_t p, std::vector<complex> data):
    p_(p), n_(n), data_(data)
{

}

ComplexMatrixGroup::ComplexMatrixGroup(size_t n, size_t p):
    p_(p), n_(n), data_((p-1)*n*n, complex(0))
{

}

static void print_fixed_double(double value, int int_digits, int frac_digits) {
    // 总宽度 = 符号(1) + 整数部分位数 + 小数点(1) + 小数位数
    int total_width = 1 + int_digits + 1 + frac_digits;
    printf("%+0*.*f", total_width, frac_digits, value);
}

void ComplexMatrixGroup::print(const char* name) const
{
    if(name)
        printf("[print ComplexMatrixGroup %s]\n", name);
    else
        printf("[print ComplexMatrixGroup]\n");

    
    for(int w=0; w<p_-1; w++)
    {
        printf("[w=%d]\n", w);
        for(int x=0; x<n_; x++)
        {
            for(int y=0; y<n_; y++)
            {
                print_fixed_double(std::real(at(w, x, w)), 3, 2);
                printf(" ");
                // printf("%.2lf ", std::real(at(w, x, w)));
            }
            printf("\t\t");
            for(int y=0; y<n_; y++)
            {
                print_fixed_double(std::imag(at(w, x, w)), 3, 2);
                printf(" ");
                // printf("%.2lf ", std::imag(at(w, x, w)));
            }
            printf("\n");
        }
        
    }
    printf("[print ComplexMatrixGroup finish]\n");
}