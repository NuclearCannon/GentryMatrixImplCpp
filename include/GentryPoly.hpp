#pragma once
#include "montgomery.hpp"
#include "GPU/cuda_buffer.hpp"
#include <cstdint>
#include <vector>
#include "ntt.hpp"
#include <unordered_map>
#include <variant>
#include <memory>

enum class GPDevice {
    CPU, CUDA
};

class GPComponent;
class GPComponentCuda;

class GPCCtx
{
private:
    size_t n_, p_;
    uint64_t q_;
    // 用于I轴NTT
    std::unique_ptr<NTTerI> ntter_i_;
    std::unique_ptr<TwistedNtterXY> ntter_p_, ntter_n_;
    std::unique_ptr<TwistedNtterW> ntter_w_;
public:
    GPCCtx(size_t n, size_t p, uint64_t q, uint64_t qroot);
    GPCCtx(GPCCtx&&);
    friend class GPComponent;
    friend class GPComponentCuda;
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
    friend class GPComponentCuda;
    friend class GentryPoly;
    // 构造函数，未初始化数据
    GPComponent(size_t n, size_t p, uint64_t q);

    static GPComponent zeros(size_t n, size_t p, uint64_t q) {
        return GPComponent(n, p, q);
    }

    // 从给定数据中构造
    // TODO: 要不要改成std::option<GPComponent>之类的返回值？
    static GPComponent from_data(size_t n, size_t p, uint64_t q, std::vector<uint64_t> data);
    static GPComponent from_signed_data(size_t n, size_t p, uint64_t q, const std::vector<int64_t>& data);
    // 析构函数：默认即可
    ~GPComponent();
    // 允许复制（以便塞进容器），但是用户应该慎用它
    GPComponent(const GPComponent&);
    // 允许移动（以便塞进容器）
    GPComponent(GPComponent&&);
    GPComponent& operator=(GPComponent&&) = default;
    // 不允许赋值
    GPComponent& operator=(const GPComponent&) = delete;
    

    // setters: 
    void set_from_data(const std::vector<uint64_t>& data);

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

    std::vector<int64_t> to_signed() const;

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
    friend class GentryPoly;
    // 构造函数，未初始化数据
    GPComponentCuda(size_t n, size_t p, uint64_t q);

    // 析构函数：默认即可
    ~GPComponentCuda();
    // 允许复制（以便塞进容器），但是用户应该慎用它
    GPComponentCuda(const GPComponentCuda&);
    // 允许移动（以便塞进容器）
    GPComponentCuda(GPComponentCuda&&);
    // 不允许赋值
    GPComponentCuda& operator=(const GPComponentCuda&) = delete;
    GPComponentCuda& operator=(GPComponentCuda&&) = default;

    // setters
    void set_from_other(const GPComponentCuda& other);
    void set_from_cuda_buffer(const CudaBuffer&);
    void set_from_cpu(const GPComponent& other);
    static GPComponentCuda copy_from_cpu(const GPComponent&);

    // getters
    inline size_t get_n() const {return n_;}
    inline size_t get_p() const {return p_;}
    inline size_t get_size() const {return data_.size()/sizeof(uint64_t);}
    inline uint64_t get_q() const {return q_;}
    inline const MontgomeryMultiplier& get_mm() const {return mm_;}

    void to_buffer(uint64_t* data) const;
    void to_cpu(GPComponent& other) const;
    
    inline bool like(const GPComponentCuda& other) const
    {
        return (n_ == other.n_) && (p_ == other.p_) && (q_ == other.q_);
    }
    inline bool like(const GPComponent& other) const
    {
        return (n_ == other.n_) && (p_ == other.p_) && (q_ == other.q_);
    }
    inline bool like(const GPCCtx& other) const
    {
        return (n_ == other.n_) && (p_ == other.p_) && (q_ == other.q_);
    }

    // 定义一些无需上下文的运算

    static void neg(GPComponentCuda& dst, const GPComponentCuda& src);
    static void add(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2);
    static void sub(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2);
    static void mul(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2);
    static void mul_scalar(GPComponentCuda& dst, const GPComponentCuda& src1, uint64_t src_scalar);
    static void mont_encode(GPComponentCuda& dst, const GPComponentCuda& src);
    static void mont_decode(GPComponentCuda& dst, const GPComponentCuda& src);
    static void mul_mont(GPComponentCuda& dst, const GPComponentCuda& src1, const GPComponentCuda& src2);

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


};

class GentryPolyCtx {
public:
    using QRootPair = std::pair<uint64_t, uint64_t>;

    GentryPolyCtx(size_t n, size_t p, const std::vector<QRootPair>& q_and_qroots);

    const class GPCCtx& get_ctx(uint64_t q) const;

    size_t n() const { return n_; }
    size_t p() const { return p_; }

private:
    size_t n_, p_;
    std::unordered_map<uint64_t, GPCCtx> ctx_map_;
};

// Zq[i,X,Y,W]/<i^2+1, X^n-i, Y^n+i, Phi_p(W)>
class GentryPoly {
public:

    // 从现有 components 构造（内部使用）
    // TODO: 这俩构造函数的is_cuda字段真是多余
    // TODO: 为什么不能从n,p,moduli中构造呢？
    GentryPoly(std::vector<uint64_t> moduli,
               std::vector<GPComponent> cpu_comps);
    GentryPoly(std::vector<uint64_t> moduli,
               std::vector<GPComponentCuda> cuda_comps);

    // 用户指定值构造（仅 CPU）
    static GentryPoly from_coeffs(
        size_t n, size_t p,
        const std::vector<uint64_t>& moduli,
        const std::vector<std::vector<uint64_t>>& coeffs_mod_q
    );

    // 设备查询
    bool is_cuda() const { return device_ == GPDevice::CUDA; }
    bool is_cpu() const { return device_ == GPDevice::CPU; }

    // 属性
    size_t n() const;
    size_t p() const;
    const std::vector<uint64_t>& moduli() const { return moduli_; }

    // like 判断：same device, n, p, moduli (order-sensitive)
    bool like(const GentryPoly& other) const;

    // ===== 静态运算 =====
private:
    using CpuOp2 = void (*)(GPComponent& , const GPComponent& );
    using CpuOp3 = void (*)(GPComponent& , const GPComponent& , const GPComponent& );
    using CudaOp2 = void (*)(GPComponentCuda& , const GPComponentCuda& );
    using CudaOp3 = void (*)(GPComponentCuda& , const GPComponentCuda& , const GPComponentCuda& );

    template<CpuOp2 cpu_op, CudaOp2 cuda_op>
    static void _op2_tmpl(GentryPoly& dst, const GentryPoly& src);
    template<CpuOp3 cpu_op, CudaOp3 cuda_op>
    static void _op3_tmpl(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2);
public:
    static void neg(GentryPoly& dst, const GentryPoly& src);
    static void add(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2);
    static void sub(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2);
    static void mul(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2);
    static void mul_scalar(GentryPoly& dst, const GentryPoly& src1, uint64_t src_scalar);
    static void mont_encode(GentryPoly& dst, const GentryPoly& src);
    static void mont_decode(GentryPoly& dst, const GentryPoly& src);
    static void mul_mont(GentryPoly& dst, const GentryPoly& src1, const GentryPoly& src2);

    // NTT / INTT：原地，需传入 ctx_set
private:
    using CpuNttOp = void (GPComponent::*)(const GPCCtx &ctx);
    using CudaNttOp = void (GPComponentCuda::*)(const GPCCtx &ctx);

    template<CpuNttOp cpu_op, CudaNttOp cuda_op>
    void _ntt_tmpl(const GentryPolyCtx&);

public:

    void i_ntt(const GentryPolyCtx&);
    void i_intt(const GentryPolyCtx&);
    void w_ntt(const GentryPolyCtx&);
    void w_intt(const GentryPolyCtx&);
    void x_ntt(const GentryPolyCtx&);
    void x_intt(const GentryPolyCtx&);
    void y_ntt(const GentryPolyCtx&);
    void y_intt(const GentryPolyCtx&);
    inline void iw_ntt(const GentryPolyCtx& ctx) {
        i_ntt(ctx); w_ntt(ctx);
    }

    inline void iw_intt(const GentryPolyCtx& ctx) {
        w_intt(ctx); i_intt(ctx);
    }

    inline void xy_ntt(const GentryPolyCtx& ctx) {
        x_ntt(ctx); y_ntt(ctx);
    }

    inline void xy_intt(const GentryPolyCtx& ctx) {
        y_intt(ctx); x_intt(ctx);
    }

    inline void ntt(const GentryPolyCtx& ctx) {
        i_ntt(ctx); w_ntt(ctx); x_ntt(ctx); y_ntt(ctx);
    }

    inline void intt(const GentryPolyCtx& ctx) {
        y_intt(ctx); x_intt(ctx); w_intt(ctx); i_intt(ctx);
    }

private:
    GPDevice device_;
    std::vector<uint64_t> moduli_;

    using CpuStorage = std::vector<GPComponent>;
    using CudaStorage = std::vector<GPComponentCuda>;
    std::variant<CpuStorage, CudaStorage> storage_;

    // 辅助访问器
    const CpuStorage& cpu_components() const;
    CpuStorage& cpu_components();
    const CudaStorage& cuda_components() const;
    CudaStorage& cuda_components();

    // 内部 like 检查（不检查 device）
    bool like_no_device(const GentryPoly& other) const;

public:
    // 下面定义一些工厂函数。所有的工厂函数都返回CPU版本的结果。

    
    static GentryPoly zeros(size_t n, size_t p, const std::vector<uint64_t>& moduli, GPDevice dev = GPDevice::CPU);
    static GentryPoly zeros_like(const GentryPoly& other, GPDevice dev = GPDevice::CPU) {
        return zeros(other.n(), other.p(), other.moduli(), dev);
    }
    static GentryPoly dg(size_t n, size_t p, const std::vector<uint64_t>& moduli);
    static GentryPoly sk(size_t n, size_t p, const std::vector<uint64_t>& moduli);
    static GentryPoly uniform(size_t n, size_t p, const std::vector<uint64_t>& moduli);
    static GentryPoly randint(size_t n, size_t p, const std::vector<uint64_t>& moduli, int lb, int ub);

    // 
    int64_t abs() const;

    void moduli_extend_mult(uint64_t mod);
    void moduli_extend_unsafe(uint64_t mod);

    // 这会将自己按照KS的要求切分为多个分量
    // 运行结束后，自己的取值会被摧毁
    std::vector<std::vector<uint64_t>> split_by_moduli();
    std::vector<CudaBuffer> split_by_moduli_cuda();

    // 移除模数链中的一个模数，并且使得自己的值除以它
    void moduli_reduce(uint64_t modulus);

    void set_from_vec64(const std::vector<uint64_t>&);
    void set_from_cuda_buffer(const CudaBuffer&);

    GentryPoly to_cpu() const;
    GentryPoly to_cuda() const;
    bool eq(const GentryPoly&) const;

    // 其实就是所有元素的简单加和
    uint64_t hash() const;
};