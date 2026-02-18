#pragma once
#include "montgomery.hpp"
#include "GPU/cuda_buffer.hpp"
#include <cstdint>
#include <vector>
#include "ntt.hpp"

class GPComponent;

class GPCCtx
{
private:
    size_t n_, p_;
    uint64_t q_;
    // 用于I轴NTT
    NTTerI ntter_i_;
    TwistedNtterXY ntter_p_, ntter_n_;
    TwistedNtterW ntter_w_;
public:
    GPCCtx(size_t n, size_t p, uint64_t q, uint64_t qroot);

    friend class GPComponent;
};

// GentryPoly在某个CRT上的分量
class GPComponent
{
private:
    struct tag_no_data {};
    
    size_t n_,p_; // 描述尺寸
    uint64_t q_;    // 模数
    MontgomeryMultiplier mm_;   // 私有蒙哥马利乘法器
    std::vector<uint64_t> data_;    // 数据缓冲区

    GPComponent(size_t n, size_t p, uint64_t q, tag_no_data);

public:
    // 构造函数，未初始化数据
    GPComponent(size_t n, size_t p, uint64_t q);

    // 从给定数据中构造
    // TODO: 要不要改成std::option<GPComponent>之类的返回值？
    static GPComponent from_data(size_t n, size_t p, uint64_t q, std::vector<uint64_t> data);
    // 析构函数：默认即可
    ~GPComponent();
    // 允许复制（以便塞进容器），但是用户应该慎用它
    GPComponent(const GPComponent&);
    // 允许移动（以便塞进容器）
    GPComponent(GPComponent&&);
    // 不允许赋值
    GPComponent& operator=(const GPComponent&) = delete;
    GPComponent& operator=(GPComponent&&) = delete;

    // setters: 无

    // getters
    inline size_t get_n() const {return n_;}
    inline size_t get_p() const {return p_;}
    inline size_t get_size() const {return data_.size();}
    inline uint64_t get_q() const {return q_;}
    inline const MontgomeryMultiplier& get_mm() const {return mm_;}
    inline const std::vector<uint64_t>& get_data() const {return data_;}
    inline bool like(const GPComponent& other) const
    {
        return (n_ == other.n_) && (p_ == other.p_) && (q_ == other.q_);
    }
    inline bool like(const GPCCtx& other) const
    {
        return (n_ == other.n_) && (p_ == other.p_) && (q_ == other.q_);
    }

    // 定义一些无需上下文的运算

    static void neg(GPComponent& dst, const GPComponent& src);
    static void add(GPComponent& dst, const GPComponent& src1, const GPComponent& src2);
    static void sub(GPComponent& dst, const GPComponent& src1, const GPComponent& src2);
    static void mul(GPComponent& dst, const GPComponent& src1, const GPComponent& src2);
    static void mul_scalar(GPComponent& dst, const GPComponent& src1, uint64_t src_scalar);
    static void mont_encode(GPComponent& dst, const GPComponent& src);
    static void mont_decode(GPComponent& dst, const GPComponent& src);
    static void mul_mont(GPComponent& dst, const GPComponent& src1, const GPComponent& src2);

    // 定义一些需要上下文的运算
    // 由于NTT天然适合原地操作，只提供原地版本
    void i_ntt(const GPCCtx&);
    void i_intt(const GPCCtx&);
    void w_ntt(const GPCCtx&);
    void w_intt(const GPCCtx&);
    void x_ntt(const GPCCtx&);
    void x_intt(const GPCCtx&);
    void y_ntt(const GPCCtx&);
    void y_intt(const GPCCtx&);

    inline void iw_ntt(const GPCCtx& ctx) {
        i_ntt(ctx); w_ntt(ctx);
    }

    inline void iw_intt(const GPCCtx& ctx) {
        w_intt(ctx); i_intt(ctx);
    }

    inline void xy_ntt(const GPCCtx& ctx) {
        x_ntt(ctx); y_ntt(ctx);
    }

    inline void xy_intt(const GPCCtx& ctx) {
        y_intt(ctx); x_intt(ctx);
    }

    inline void ntt(const GPCCtx& ctx) {
        i_ntt(ctx); w_ntt(ctx); x_ntt(ctx); y_ntt(ctx);
    }

    inline void intt(const GPCCtx& ctx) {
        y_intt(ctx); x_intt(ctx); w_intt(ctx); i_intt(ctx);
    }

    bool eq(const GPComponent&) const;


};

// GentryPoly在某个CRT上的分量（CUDA）
class GPComponentCuda
{
private:
    size_t n_,p_; // 描述尺寸
    uint64_t q_;    // 模数
    MontgomeryMultiplier mm_;   // 私有蒙哥马利乘法器
    CudaBuffer data_;    // 数据缓冲区
public:
};

// Zq[i,X,Y,W]/<i^2+1, X^n-i, Y^n+i, Phi_p(W)>
class GentryPoly;