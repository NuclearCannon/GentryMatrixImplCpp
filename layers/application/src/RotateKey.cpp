#include "application.hpp"
#include <iostream>

void RotateKey::test_pt_rotate()
{
    // 准备各个参数
    int n = 8;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };
    double delta = 10000;
    // 确定rotate步长
    int x_bias = 5;
    int y_bias = 3;

    GentryPolyCtx ctx(n, p, qrp);

    // 准备明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(5, n, p);
    // for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++) mat1.at(w, x, y) = w*x;
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    Plaintext pt2 = pt1.rotate_XY(x_bias, y_bias);
    // 恢复成矩阵
    ComplexMatrixGroup mat2 = pt2.to_cmat(delta);

    // mat1.print("mat1");
    // mat2.print("mat2");

    double max_diff = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex diff = mat2.at(w,x,y) - mat1.at(w,(x+x_bias)%n,(y+y_bias)%n); // 在这里进行编码前层面上的转置操作
        double diff_abs = std::abs(diff);
        if(diff_abs>max_diff)max_diff = diff_abs;
    }
    std::cout << "test_pt_rotate: max_diff="<<max_diff<<std::endl;
}