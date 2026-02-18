#pragma once
#include "montgomery.hpp"
#include "GPU/cuda_buffer.hpp"
#include <cstdint>
#include <vector>


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
    ~GPComponent() = default;
    // 允许复制（以便塞进容器），但是用户应该慎用它
    GPComponent(const GPComponent&) = default;
    GPComponent& operator=(const GPComponent&) = default;
    // 允许移动（以便塞进容器）
    GPComponent(GPComponent&&) = default;
    GPComponent& operator=(GPComponent&&) = default;

    // setters: 无

    // getters
    inline size_t get_n() const {return n_;}
    inline size_t get_p() const {return p_;}
    inline size_t get_size() const {return data_.size();}
    inline uint64_t get_q() const {return q_;}
    inline const MontgomeryMultiplier& get_mm() const {return mm_;}
    inline const std::vector<uint64_t>& get_data() const {return data_;}
    inline bool like(const GPComponent& other)
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