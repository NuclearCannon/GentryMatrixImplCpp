#pragma once
#include "complex_matrix.hpp"
#include <memory>
#include "GentryPoly.hpp"
#include "FHE/key_switch_gp.hpp"

class Ciphertext;
class CircledastKey;
class ConjTransposeKey;
class ConjKey;
class MultKey;
class RotateKey;

class Plaintext
{
    friend class Ciphertext;
    friend class ConjTransposeKey;
    friend class CircledastKey;
private:
    std::unique_ptr<GentryPoly> data_;

    Plaintext(std::unique_ptr<GentryPoly> data);

    Plaintext(Plaintext&&) = default;

public:
    ~Plaintext() = default;

    // 从复矩阵组中编码一个明文出来
    // n, p直接从复矩阵组中读取，模数链则需要外界传递进来
    static Plaintext from_cmat(const ComplexMatrixGroup&, const std::vector<uint64_t>& mods, double delta);
    static Plaintext from_scalar(int n, int p, complex scalar, const std::vector<uint64_t>& mods, double delta);

    static Plaintext _from_cmat_without_encoding(const ComplexMatrixGroup&, const std::vector<uint64_t>& mods, double delta);

    // 从明文对象中恢复出复矩阵
    ComplexMatrixGroup to_cmat(double delta) const;
    ComplexMatrixGroup _to_cmat_without_decoding(double delta) const;

    Plaintext circledast(const Plaintext& other, const GentryPolyCtx& ctx) const;
    Plaintext conj_transpose(const GentryPolyCtx& ctx) const;
    Plaintext conj(const GentryPolyCtx& ctx) const;
    Plaintext rotate(int x_bias, int y_bias, int w_bias) const;

// ===========================以下是单元测试==========================
    static void test_pt_encode_and_decode();
};



class SecretKey
{
    // friend class Ciphertext;
    // friend class CircledastKey;
    // friend class ConjTransposeKey;
    // friend class MultKey;
    // friend class RotateKey;
private:
    int n_, p_;
    std::vector<int> sk_data_;
public:
    SecretKey(size_t n, size_t p);

    ~SecretKey() = default;

    GentryPoly as_poly(const std::vector<uint64_t>& mods) const;

    
};




class Ciphertext
{
    friend class CircledastKey;
    friend class ConjTransposeKey;
    friend class ConjKey;
    friend class MultKey;
    friend class RotateKey;
private:
    std::unique_ptr<GentryPoly> a_, b_;
    Ciphertext(std::unique_ptr<GentryPoly> a, std::unique_ptr<GentryPoly> b);

public:
    static Ciphertext encrypt(const Plaintext& pt, const SecretKey& sk, const GentryPolyCtx& ctx);
    static Ciphertext zeros(int n, int p, std::vector<uint64_t> mods);
    Plaintext decrypt(const SecretKey& sk, const GentryPolyCtx& ctx);
    double check_abs(const SecretKey& sk, const GentryPolyCtx& ctx, double delta = 1);
    static void test_ct_encrypt_and_decrypt();

    // 这将会创建一个新的密文对象，它的数据都放在cuda上
    Ciphertext to_cuda() const;

    Ciphertext add_pt(const Plaintext& pt, const GentryPolyCtx& ctx) const;
    Ciphertext add(const Ciphertext& other) const;
    Ciphertext neg() const;
    Ciphertext sub(const Ciphertext& other) const;
    Ciphertext mul_int(int other) const;
    Ciphertext mul_i() const;
    void add_(const Ciphertext& other, const GentryPolyCtx& ctx);
    Ciphertext mul_pt(const Plaintext& pt, const GentryPolyCtx& ctx) const;

    // 自举的步骤之一：简单模数提升
    Ciphertext naive_moduli_extend(const std::vector<uint64_t>& extra_moduli) const;

    // 模数约简，用于控制噪声
    Ciphertext moduli_reduce(const std::vector<uint64_t>& moduli) const;

    static void test_ct_add_pt();
    static void test_ct_add();
    static void test_ct_mul_pt();
    static void test_ct_mul_i();
};

class CircledastKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk1_, ksk2_;
    CircledastKey(std::unique_ptr<KeySwitchKeyGP> ksk1, std::unique_ptr<KeySwitchKeyGP> ksk2);
public:
    static CircledastKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods);
    Ciphertext run(const Ciphertext& u, const Ciphertext& v, const GentryPolyCtx& ctx) const;
    Ciphertext run_cp(const Ciphertext& u, const Plaintext& v, const GentryPolyCtx& ctx) const;
    Ciphertext run_pc(const Plaintext& u, const Ciphertext& v, const GentryPolyCtx& ctx) const;

    static void test_pt_circledast_end2end();
    static void test_ct_circledast_end2end(bool cuda = false);
};

class ConjTransposeKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk_;
    ConjTransposeKey(std::unique_ptr<KeySwitchKeyGP> ksk);


public:
    static ConjTransposeKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods);
    Ciphertext run(const Ciphertext& src, const GentryPolyCtx& ctx) const;
    static void test_pt_transpose();
    static void test_ct_transpose();

};

class ConjKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk_;
    ConjKey(std::unique_ptr<KeySwitchKeyGP> ksk);


public:
    static ConjKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods);
    Ciphertext run(const Ciphertext& src, const GentryPolyCtx& ctx) const;
    static void test_pt_conj();
    static void test_ct_conj();

};


class MultKey
{
private:
    std::unique_ptr<KeySwitchKeyGP> ksk_;
    MultKey(std::unique_ptr<KeySwitchKeyGP> ksk);


public:
    static MultKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods);
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
    static RotateKey gen(const SecretKey& sk, uint64_t qo, const GentryPolyCtx& ctx, const std::vector<uint64_t>& mods, int x_bias, int y_bias, int w_bias);
    Ciphertext run(const Ciphertext& src, const GentryPolyCtx& ctx) const;
    static void test_pt_rotate();
    static void test_ct_rotate();
};