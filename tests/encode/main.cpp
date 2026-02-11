#include "complecx_matrix.hpp"
#include <iostream>


int test_encode_decode()
{
    int n=32, p=17;
    ComplexMatrixGroup mat(n,p);
    // 混乱赋值给它
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        mat.at(w,x,y) = complex(w+x, w-y) * complex(1000000);
    }
    ComplexMatrixGroup mat2 = mat.encode();
    ComplexMatrixGroup mat3 = mat2.decode();
    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat3.at(w,x,y) - mat.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_encode_decode: max_diff="<<max_diff<<std::endl;
    return 0;
}

int test_fmpz_vector()
{
    int n=32, p=17;
    ComplexMatrixGroup mat(n,p);
    // 混乱赋值给它
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        mat.at(w,x,y) = complex(w+x, w-y) * M_PI;   // 引入无理数成分以确保四舍五入误差存在
    }
    double delta = 100000;
    ComplexMatrixGroup mat2 = ComplexMatrixGroup::from_fmpz_vector(mat.to_fmpz_vector(delta), delta, n, p);
    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - mat.at(w,x,y);
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_fmpz_vector: max_diff="<<max_diff<<std::endl;
    return 0;
}

int main()
{
    test_encode_decode();
    test_fmpz_vector();
}