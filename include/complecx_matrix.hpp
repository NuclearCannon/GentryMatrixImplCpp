#pragma once
#include <cstddef>
#include <complex>
#include <vector>
#include "flints.hpp"


using complex = std::complex<double>;

// 一个(p-1, n, n)数组
// 在朴素形式下，其[w, x, y]位置是第w个矩阵的[x,y]位置
class ComplexMatrixGroup
{
private:
    size_t p_, n_;

    std::vector<complex> data_;
    // 私有构造函数，从内容中直接构造
    ComplexMatrixGroup(size_t n, size_t p, std::vector<complex> data);
public:
    ComplexMatrixGroup(size_t n, size_t p); // 0矩阵

    ComplexMatrixGroup(const ComplexMatrixGroup&) = default;
    ComplexMatrixGroup(ComplexMatrixGroup&&) = default;
    ~ComplexMatrixGroup() = default;

    ComplexMatrixGroup& operator=(const ComplexMatrixGroup&) = delete;
    ComplexMatrixGroup& operator=(ComplexMatrixGroup&&) = delete;

    inline complex& at(size_t w, size_t x, size_t y)
    {
        return data_.at((w*n_ + x)*n_ + y);
    }
    const complex& at(size_t w, size_t x, size_t y) const 
    {
        return data_.at((w*n_ + x)*n_ + y);
    }

    ComplexMatrixGroup encode() const;
    ComplexMatrixGroup decode() const;

    fmpz_vector to_fmpz_vector(double delta) const;

    static ComplexMatrixGroup from_fmpz_vector(const fmpz_vector& src, double delta, size_t n, size_t p);

    ComplexMatrixGroup matmul_ABT(const ComplexMatrixGroup& other) const;

    static ComplexMatrixGroup random(double abs_max, size_t n, size_t p);

    // 单位矩阵（调试用）
    static ComplexMatrixGroup eye(size_t n, size_t p);

    // mat[w,x,y] = x+yi
    static ComplexMatrixGroup example_matrix(size_t n, size_t p);

};


const std::vector<size_t>& get_bit_reverse_table_by_logn(size_t log2n);
