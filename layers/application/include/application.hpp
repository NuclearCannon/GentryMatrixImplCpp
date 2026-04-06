#pragma once
#include "complex_matrix.hpp"
#include <memory>
#include "GentryPoly.hpp"

class Ciphertext;

class Plaintext
{
    friend class Ciphertext;
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

class SecretKey
{
    friend class Ciphertext;
private:
    std::unique_ptr<GentryPoly> data_;
public:
    SecretKey(size_t n, size_t p, const std::vector<uint64_t>& mods);

    ~SecretKey() = default;
    
};

class Ciphertext
{
private:
    std::unique_ptr<GentryPoly> a_, b_;
    Ciphertext(std::unique_ptr<GentryPoly> a, std::unique_ptr<GentryPoly> b);

public:
    static Ciphertext encrypt(const Plaintext& pt, const SecretKey& sk, const GentryPolyCtx& ctx);
    Plaintext decrypt(const SecretKey& sk, const GentryPolyCtx& ctx);

    static void test_ct_encrypt_and_decrypt();
};