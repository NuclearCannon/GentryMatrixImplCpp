#pragma once
#include "complex_matrix.hpp"
#include <memory>
#include "GentryPoly.hpp"
#include "FHE/key_switch_gp.hpp"

class Ciphertext;
class CircledastKey;
class ConjTransposeKey;
class MultKey;
class RotateKey;

class Plaintext
{
    friend class Ciphertext;
    friend class ConjTransposeKey;
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

    Plaintext circledast(const Plaintext& other, const GentryPolyCtx& ctx) const;
    Plaintext conj_transpose(const GentryPolyCtx& ctx) const;
    Plaintext rotate(int x_bias, int y_bias, int w_bias) const;

// ===========================以下是单元测试==========================
    static void test_pt_encode_and_decode();
};



class SecretKey
{
    friend class Ciphertext;
    friend class CircledastKey;
    friend class ConjTransposeKey;
    friend class MultKey;
    friend class RotateKey;
private:
    std::unique_ptr<GentryPoly> data_;
public:
    SecretKey(size_t n, size_t p, const std::vector<uint64_t>& mods);

    ~SecretKey() = default;
    
};




class Ciphertext
{
    friend class CircledastKey;
    friend class ConjTransposeKey;
    friend class MultKey;
    friend class RotateKey;
private:
    std::unique_ptr<GentryPoly> a_, b_;
    Ciphertext(std::unique_ptr<GentryPoly> a, std::unique_ptr<GentryPoly> b);

public:
    static Ciphertext encrypt(const Plaintext& pt, const SecretKey& sk, const GentryPolyCtx& ctx);
    Plaintext decrypt(const SecretKey& sk, const GentryPolyCtx& ctx);

    static void test_ct_encrypt_and_decrypt();

    // 这将会创建一个新的密文对象，它的数据都放在cuda上
    Ciphertext to_cuda() const;
};

class CircledastKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk1_, ksk2_;
    CircledastKey(std::unique_ptr<KeySwitchKeyGP> ksk1, std::unique_ptr<KeySwitchKeyGP> ksk2);
public:
    static CircledastKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx);
    Ciphertext run(const Ciphertext& u, const Ciphertext& v, const GentryPolyCtx& ctx) const;

    static void test_pt_circledast_end2end();
    static void test_ct_circledast_end2end(bool cuda = false);
};

class ConjTransposeKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk_;
    ConjTransposeKey(std::unique_ptr<KeySwitchKeyGP> ksk);


public:
    static ConjTransposeKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx);
    Ciphertext run(const Ciphertext& src, const GentryPolyCtx& ctx) const;
    static void test_pt_transpose();
    static void test_ct_transpose();

};


class MultKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk_;
    MultKey(std::unique_ptr<KeySwitchKeyGP> ksk);


public:
    static MultKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx);
    Ciphertext run(const Ciphertext& src1, const Ciphertext& src2, const GentryPolyCtx& ctx) const;
    // static void test_pt_transpose(); // TODO: 如果有需要再做
    static void test_ct_mult();

};

class RotateKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk_;
    int x_bias_, y_bias_, w_bias_;
    RotateKey(std::unique_ptr<KeySwitchKeyGP> ksk, int x_bias, int y_bias, int w_bias);

public:
    static RotateKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, int x_bias, int y_bias, int w_bias);
    Ciphertext run(const Ciphertext& src, const GentryPolyCtx& ctx) const;
    static void test_pt_rotate();
    static void test_ct_rotate();
};