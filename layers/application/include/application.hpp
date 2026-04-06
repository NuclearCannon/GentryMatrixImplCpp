#pragma once
#include "complex_matrix.hpp"
#include <memory>
#include "GentryPoly.hpp"

class Plaintext
{
private:
    std::unique_ptr<GentryPoly> data_;

    Plaintext(std::unique_ptr<GentryPoly> data);

    Plaintext(Plaintext&&) = default;

public:
    ~Plaintext() = default;

    // 从复矩阵组中编码一个明文出来
    // n, p直接从复矩阵组中读取，模数链则需要外界传递进来
    static Plaintext from_cmat(const ComplexMatrixGroup&, const std::vector<uint64_t>& mods, double delta);

    // 从明文对象中恢复出复矩阵
    ComplexMatrixGroup to_cmat(double delta) const;



// ===========================以下是单元测试==========================
    static void test_pt_encode_and_decode();
};