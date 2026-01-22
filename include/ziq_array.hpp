#pragma once
#include <vector>
#include <stdexcept>
#include "flints.hpp"
#include "twisted_ntter.hpp"


class MatmulContext {
    int n_;
    fmpz_t q_;
    fmpz_mod_mat_t A_mat, B_mat, C_mat; // 矩阵乘法时使用的缓冲区
public:
    MatmulContext(int n, fmpz_t q);
    ~MatmulContext();

    // C = A @ B.T  (mod q)
    void matmul_transpose(fmpz* C, const fmpz* A, const fmpz* B);
};

class ZiqArray;

class ZiqArrayContext {
private:
    int n_, p_;
    int nn_, pnn_, size_;
    fmpz_t q_, I_, I_inv_;
    fmpz_mod_ctx_t q_ctx_;
    TwistedNtterXY *ntter_p, *ntter_n;    // 这两个ntter会分别负责模X^n-I的和模X^n+I的
    TwistedNtterW *ntter_w;
    MatmulContext *mm_ctx;

    fmpz_vector *buf_p, *buf_n, *buf_half, *buf_size;
public:

    ZiqArrayContext(int n, int p, int g, fmpz_t q, fmpz_t zeta, fmpz_t eta);

    ~ZiqArrayContext();

    fmpz_vector iw_ntt(const fmpz_vector& src) const;
    fmpz_vector iw_intt(const fmpz_vector& src) const;

    fmpz_vector xy_ntt(const fmpz_vector& src) const;
    fmpz_vector xy_intt(const fmpz_vector& src) const;

    fmpz_vector circledast(const fmpz_vector& A, const fmpz_vector& B) const;

    
    // 逐位加法
    void add(fmpz_vector& dst, const fmpz_vector src1, const fmpz_vector& src2) const;
    // 逐位减法
    void sub(fmpz_vector& dst, const fmpz_vector src1, const fmpz_vector& src2) const;
    // 逐位乘法
    void mul(fmpz_vector& dst, const fmpz_vector src1, const fmpz_vector& src2) const;
    // 逐位负
    void neg(fmpz_vector& dst, const fmpz_vector src1) const;
    // 标量乘
    void mul_scalar(fmpz_vector& dst, const fmpz_vector src_vec, const fmpz_t src_scalar) const;
    // 比较
    bool eq(const fmpz_vector src1, const fmpz_vector& src2) const;

    ZiqArray zeros() const;
    ZiqArray uniform() const;
    ZiqArray dg() const;

    inline const fmpz* q() const { return q_; }

};

class ZiqArray {
private:
    const ZiqArrayContext* ctx_;
    int size_;
    fmpz_vector data_;

public:
    
    // 构造函数: 全零数组
    ZiqArray(int n, int p, const ZiqArrayContext* ctx);

    // 构造函数: 从 vector 构造（需正确大小）
    ZiqArray(fmpz_vector data, const ZiqArrayContext* ctx);

    ~ZiqArray();

    // 禁用赋值（保持 immutable）
    ZiqArray& operator=(const ZiqArray&) = delete;
    ZiqArray& operator=(ZiqArray&&) = delete;

    ZiqArray(const ZiqArray&) = delete; // 原则上不应该允许复制
    ZiqArray(ZiqArray&&);   // 可以移动

    // 运算符重载：逐元素操作
    ZiqArray add(const ZiqArray& other) const;
    ZiqArray neg() const;
    ZiqArray sub(const ZiqArray& other) const;
    ZiqArray mul(const ZiqArray& other) const;
    ZiqArray mul_scalar(const fmpz_t other) const;

    bool eq(const ZiqArray& other) const;

    ZiqArray iw_ntt() const;
    ZiqArray iw_intt() const;
    ZiqArray xy_ntt () const;
    ZiqArray xy_intt() const;
    ZiqArray circledast(const ZiqArray& other) const;

    inline const fmpz_vector& data() const { return data_; };
    inline const ZiqArrayContext* ctx() const { return ctx_; };

};