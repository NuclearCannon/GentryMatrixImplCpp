#pragma once
#include "flints.hpp"

class TwistedNtterXY
{
private:
    int n_;  // 做模X^n-I或者X^n+I的NTT
    fmpz_t q_; // 模数
    fmpz_t I_; // I^2+1==0 (mod q),  I=zeta^n
    fmpz_t zeta_;  // 4n次本原单位根
    fmpz_t omega_; // omega=zeta^4, n次本原单位根
    fmpz_t omega_inv_; // omega的乘法逆元
    fmpz_mod_ctx_t q_ctx_; // “模上下文”，它的存在有利于快速模运算
    // zeta的正数次幂和负数次幂，质数范围是[0..n]
    // 注意：zeta_neg_pows实际上还乘以了n_inv
    fmpz_vector zeta_pos_pows_, zeta_neg_pows_ninv_;
    fmpz_vector buffer_;    // 私有缓冲区

public:
    TwistedNtterXY(int n, fmpz_t q, fmpz_t zeta);
    ~TwistedNtterXY();

    // 不是const，因为会修改buffer
    void ntt(const fmpz_vector& src, fmpz_vector& dst);
    void intt(const fmpz_vector& src, fmpz_vector& dst);

};




class TwistedNtterW
{
private:
    int p_;      // 小质数长度
    int g_;      // p的生成元
    fmpz_t q_;   // 模数
    fmpz_t p_inv_;   // 模数
    fmpz_mod_ctx_t q_ctx_;   // 模上下文，用于加速模计算
    std::vector<int> gpp_;   // g的正指数序列
    std::vector<int> gnp_;   // g的负指数序列
    fmpz_t eta_;   // p-本原单位根
    fmpz_mod_poly_t a_poly_;   // 序列a（用作缓冲区）
    fmpz_mod_poly_t b_poly_;   // 序列b的多项式形式
    fmpz_mod_poly_t b2_poly_;   // 序列b的多项式形式
    fmpz_mod_poly_t c_poly_;   // 序列c=a*b（用作缓冲区）
    fmpz_mod_poly_t xn1_;   // 常多项式X^{p-1}-1

    fmpz_vector buf1, buf2;

public:
    TwistedNtterW(int p,int g,fmpz_t q,fmpz_t eta);
    ~TwistedNtterW();

    // 不是const，因为会修改buffer
    void ntt(const fmpz_vector& src, fmpz_vector& dst);
    void intt(const fmpz_vector& src, fmpz_vector& dst);

private:
    
};
