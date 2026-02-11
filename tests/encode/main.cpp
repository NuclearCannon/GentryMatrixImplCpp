#include "complecx_matrix.hpp"
#include <iostream>


int main()
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
    std::cout << "test encode: max_diff="<<max_diff<<std::endl;
    return 0;
}