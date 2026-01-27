#pragma once
#include <vector>
#include <stdexcept>
#include "flints.hpp"
#include "twisted_ntter.hpp"
#include "uint64.hpp"
#include <memory>

class TwistedNtterXY64
{
private:
    int n_;
    u64 q_;
    vec64 zeta_pos_pows_;
    vec64 zeta_neg_pows_ninv_;
    vec64 omega_pos_pows_;
    vec64 omega_neg_pows_;

public:
    TwistedNtterXY64(int n, u64 q, u64 zeta);
    ~TwistedNtterXY64();
    void ntt(vec64& dst, const vec64& src) const;
    void intt(vec64& dst, const vec64& src) const;
};




class TwistedNtterW64
{
private:
    int p_;
    u64 q_;
    u64 p_inv_; // p的乘法逆元
    vec64 eta_powers_;

    void _ntt_no_alias(vec64& dst, const vec64& src) const;
    void _intt_no_alias(vec64& dst, const vec64& src) const;
public:
    TwistedNtterW64(int p , u64 q, u64 eta);
    ~TwistedNtterW64();
    void ntt(vec64& dst, const vec64& src) const;
    void intt(vec64& dst, const vec64& src) const;
};


class U64Context {
private:
    int n_, p_;
    int nn_, pnn_, size_;
    u64 q_, I_, I_inv_;
    std::unique_ptr<const TwistedNtterXY64> ntter_p, ntter_n;    // 这两个ntter会分别负责模X^n-I的和模X^n+I的
    std::unique_ptr<const TwistedNtterW64> ntter_w;
public:

    U64Context(int n, int p, u64 q, u64 root_q);

    ~U64Context();

    void iw_ntt(vec64& dst, const vec64& src) const;
    void iw_intt(vec64& dst, const vec64& src) const;
    void xy_ntt(vec64& dst, const vec64& src) const;
    void xy_intt(vec64& dst, const vec64& src) const;

    // 逐位加法
    void add(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位减法
    void sub(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位乘法
    void mul(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位负
    void neg(vec64& dst, const vec64& src1) const;
    // 标量乘
    void mul_scalar(vec64& dst, const vec64& src_vec, u64 src_scalar) const;
    // 比较
    bool eq(const vec64& src1, const vec64& src2) const;


    inline u64 q() const { return q_; }
    inline int get_n() const {return n_; }
    inline int get_p() const {return p_; }
    inline int get_size() const {return size_; }
    

};



class U64CtxChain 
{
private:
    int n_, p_;   // 多项式尺寸
    size_t chain_len_; // 下面几个vector的长度
    // 模数链和对应的单位根
    vec64 mods_, roots_;
    std::vector<std::shared_ptr<const U64Context>> ctxs_;
    
public:
    U64CtxChain(int n, int p, const vec64& mods, const vec64& roots);
    ~U64CtxChain();

    void iw_ntt(vv64& dst, const vv64& src) const;
    void iw_intt(vv64& dst, const vv64& src) const;
    void xy_ntt(vv64& dst, const vv64& src) const;
    void xy_intt(vv64& dst, const vv64& src) const;

    // 逐位加法
    void add(vv64& dst, const vv64& src1, const vv64& src2) const;
    // 逐位减法
    void sub(vv64& dst, const vv64& src1, const vv64& src2) const;
    // 逐位乘法
    void mul(vv64& dst, const vv64& src1, const vv64& src2) const;
    // 逐位负
    void neg(vv64& dst, const vv64& src1) const;
    // 标量乘
    void mul_scalar(vv64& dst, const vv64& src_vec, u64 src_scalar) const;
    // 比较
    bool eq(const vv64& src1, const vv64& src2) const;

    inline int get_chain_length() const {return chain_len_; }
    inline int get_n() const {return n_; }
    inline int get_p() const {return p_; }
    inline int get_size() const {return 2*n_*n_*(p_-1); }
    inline const vec64& get_mods() const {return mods_; }

};

// ZiqArray替代品
class CRTArray
{
private:
    std::shared_ptr<const U64CtxChain> cc_;
    vv64 data_;

public:
    
    // 构造函数: 全零数组
    CRTArray(std::shared_ptr<const U64CtxChain> cc);

    // 构造函数: 从 vector 构造（需正确大小）
    CRTArray(vv64 data, std::shared_ptr<const U64CtxChain> cc);

    static CRTArray from_fmpz_vector(const fmpz_vector& data, std::shared_ptr<const U64CtxChain> cc);
    void to_fmpz_vector(fmpz_vector& dst);

    ~CRTArray();

    // 禁用赋值（保持 immutable）
    CRTArray& operator=(const CRTArray&) = delete;
    CRTArray& operator=(CRTArray&&) = delete;

    CRTArray(const CRTArray&); // 可以复制
    CRTArray(CRTArray&&);   // 可以移动

    // 运算符重载：逐元素操作
    CRTArray add(const CRTArray& other) const;
    CRTArray neg() const;
    CRTArray sub(const CRTArray& other) const;
    CRTArray mul(const CRTArray& other) const;
    CRTArray mul_scalar(u64 other) const;
    CRTArray mul_poly(const CRTArray& other) const;

    bool eq(const CRTArray& other) const;

    CRTArray iw_ntt() const;
    CRTArray iw_intt() const;
    CRTArray xy_ntt() const;
    CRTArray xy_intt() const;
    CRTArray all_ntt() const;
    CRTArray all_intt() const;
    

};