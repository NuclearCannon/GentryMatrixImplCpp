#pragma once
#include <vector>
#include <stdexcept>
#include "flints.hpp"

#include "ntt.hpp"
#include <memory>
#include "montgomery.hpp"

class TwistedNtterXY64
{
private:
    int n_;
    uint64_t q_;
    vec64 zeta_pos_pows_;
    vec64 zeta_neg_pows_;
    vec64 zeta_pos_pows_mont_;
    vec64 zeta_neg_pows_mont_;
    StandardNTTer std_ntter;
    MontgomeryMultiplier mm;

public:
    TwistedNtterXY64(int n, uint64_t q, uint64_t qroot);
    ~TwistedNtterXY64();
    void ntt_mont(vec64& dst, const vec64& src) const;
    void intt_mont(vec64& dst, const vec64& src) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;
};




class TwistedNtterW64
{
private:
    uint64_t p_;
    uint64_t q_;
    StandardNTTer subntter;
    MontgomeryMultiplier mm;

    // b_[i] = eta^{gamma^i}（蒙哥马利形式）
    vec64 b_;

    // binv_[i] = inv(etas[i])（蒙哥马利形式）
    vec64 binv_;
    
    mutable vec64 buf1;
public:
    TwistedNtterW64(int p , uint64_t q, uint64_t qroot);
    ~TwistedNtterW64();
    void ntt_mont(vec64& dst, const vec64& src) const;
    void intt_mont(vec64& dst, const vec64& src) const;

    void ntt_batch(uint64_t* dst, size_t batch_size) const;
    void intt_batch(uint64_t* dst, size_t batch_size) const;
};


class U64Context {
private:
    int n_, p_;
    int nn_, pnn_, size_;
    uint64_t q_, I_, I_inv_;
    uint64_t I_mont_, I_inv_mont_;
    uint64_t inv2_mont;
    std::unique_ptr<const TwistedNtterXY64> ntter_p, ntter_n;    // 这两个ntter会分别负责模X^n-I的和模X^n+I的
    std::unique_ptr<const TwistedNtterW64> ntter_w;
    MontgomeryMultiplier mm_;
    mutable vec64 buf_size_, bufn_;
public:

    U64Context(int n, int p, uint64_t q, uint64_t root_q);

    ~U64Context();

    void iw_ntt(vec64& dst, const vec64& src) const;
    void iw_intt(vec64& dst, const vec64& src) const;
    void xy_ntt(vec64& dst, const vec64& src) const;
    void xy_intt(vec64& dst, const vec64& src) const;

    void transpose(vec64& dst, const vec64& src) const;
    void conj(vec64& dst, const vec64& src) const;
    void w_inv(vec64& dst, const vec64& src) const;

    // 逐位加法
    void add(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位减法
    void sub(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位减法（不安全版本，这个版本下不假设src<mod）
    void sub_unsafe(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位乘法
    void mul(vec64& dst, const vec64& src1, const vec64& src2) const;
    // 逐位负
    void neg(vec64& dst, const vec64& src1) const;
    // 标量乘
    void mul_scalar(vec64& dst, const vec64& src_vec, uint64_t src_scalar) const;
    // 比较
    bool eq(const vec64& src1, const vec64& src2) const;


    inline uint64_t q() const { return q_; }
    inline int get_n() const {return n_; }
    inline int get_p() const {return p_; }
    inline int get_size() const {return size_; }
    inline const MontgomeryMultiplier& get_multiplier() const {return mm_; }
    

    void iw_ntt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const;
    void iw_intt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const;
    void xy_ntt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const;
    void xy_intt_cuda(const CudaBuffer& dst, const CudaBuffer& src) const;
};



class U64CtxChain 
{
private:
    int n_, p_;   // 多项式尺寸
    size_t chain_len_; // 下面几个vector的长度
    // 模数链和对应的单位根
    vec64 mods_, roots_;
    std::vector<std::shared_ptr<const U64Context>> ctxs_;
    fmpz_scalar mod_prod_;
    
public:
    U64CtxChain(int n, int p, const vec64& mods, const vec64& roots);
    ~U64CtxChain();

    void iw_ntt(vv64& dst, const vv64& src) const;
    void iw_intt(vv64& dst, const vv64& src) const;
    void xy_ntt(vv64& dst, const vv64& src) const;
    void xy_intt(vv64& dst, const vv64& src) const;

    void transpose(vv64& dst, const vv64& src) const;
    void conj(vv64& dst, const vv64& src) const;
    void w_inv(vv64& dst, const vv64& src) const;

    // 逐位加法
    void add(vv64& dst, const vv64& src1, const vv64& src2) const;
    // 逐位减法
    void sub(vv64& dst, const vv64& src1, const vv64& src2) const;
    // 逐位乘法
    void mul(vv64& dst, const vv64& src1, const vv64& src2) const;
    void mul_mont(vv64& dst, const vv64& src1, const vv64& src2) const;
    // 逐位负
    void neg(vv64& dst, const vv64& src1) const;
    // 标量乘
    void mul_scalar(vv64& dst, const vv64& src_vec, uint64_t src_scalar) const;
    // 比较
    bool eq(const vv64& src1, const vv64& src2) const;

    void mont_encode(vv64& dst, const vv64& src) const;
    void mont_decode(vv64& dst, const vv64& src) const;

    inline int get_chain_length() const {return chain_len_; }
    inline int get_n() const {return n_; }
    inline int get_p() const {return p_; }
    inline int get_size() const {return 2*n_*n_*(p_-1); }
    inline const vec64& get_mods() const {return mods_; }
    inline const vec64& get_roots() const {return roots_; }
    inline const fmpz* get_mod_prod() const {return mod_prod_.raw(); }
    inline const std::vector<std::shared_ptr<const U64Context>>& get_ctx() const {return ctxs_; }

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

    // 构造函数: 从raw构造
    CRTArray(const vec64& data, std::shared_ptr<const U64CtxChain> cc);

    void set_from_raw(const vec64& data);
    void set_to_zero();

    static CRTArray from_fmpz_vector(const fmpz_vector& data, std::shared_ptr<const U64CtxChain> cc);
    fmpz_vector to_fmpz_vector() const;
    fmpz_vector to_fmpz_vector_centered() const;

    ~CRTArray();

    // 禁用赋值（保持 immutable）
    CRTArray& operator=(const CRTArray&) = delete;
    CRTArray& operator=(CRTArray&&) = delete;

    CRTArray(const CRTArray&); // 可以复制
    CRTArray(CRTArray&&);   // 可以移动

    // 运算符重载：逐元素操作
    CRTArray add(const CRTArray& other) const;
    void adde(const CRTArray& other) ;
    void mul_scalar_e(uint64_t other) ;
    CRTArray neg() const;
    CRTArray sub(const CRTArray& other) const;
    CRTArray mul(const CRTArray& other) const;
    CRTArray mul_mont(const CRTArray& other) const;
    static void mul_mont3(CRTArray& dst, const CRTArray& src1, const CRTArray& src2);
    CRTArray mul_scalar(uint64_t other) const;
    CRTArray mul_poly(const CRTArray& other) const;

    bool eq(const CRTArray& other) const;

    CRTArray iw_ntt() const;
    CRTArray iw_intt() const;
    CRTArray xy_ntt() const;
    CRTArray xy_intt() const;
    CRTArray all_ntt() const;
    void all_ntt_to(CRTArray& dst) const;
    CRTArray all_intt() const;
    void all_intt_to(CRTArray& dst) const;

    static CRTArray zeros(std::shared_ptr<const U64CtxChain> cc);
    static CRTArray uniform(std::shared_ptr<const U64CtxChain> cc);
    static CRTArray dg(std::shared_ptr<const U64CtxChain> cc);

    // Y^0分量上均匀三元，其余为0
    static CRTArray sk(std::shared_ptr<const U64CtxChain> cc);
    static CRTArray randint(std::shared_ptr<const U64CtxChain> cc, int64_t start, int64_t end);
    inline std::shared_ptr<const U64CtxChain> get_cc() const { return cc_; }

    // 缩减模数链长度，主要用于Key Switch
    CRTArray mod_reduce(std::shared_ptr<const U64CtxChain> cc2) const;

    // 分解为raw形式
    vv64 mod_by_modulo() const;

    CRTArray mont_encode() const;
    CRTArray mont_decode() const;
    void mont_encode_inplace();
    void mont_decode_inplace();

    CRTArray circledast(const CRTArray& other) const;

    CRTArray transpose() const;
    CRTArray conj() const;
    CRTArray w_inv() const;

    inline const std::vector<std::vector<uint64_t>>& get_data() const {
        return data_;
    }



};


// CRTArray替代品
// 本类的接口更加着重原地操作，以this作为输出位置，而不是返回新对象

class CRTArrayGPU
{
private:
    std::shared_ptr<const U64CtxChain> cc_;
    std::vector<std::unique_ptr<CudaBuffer>> cuda_data_;
public:
    // 构造函数: 未初始化状态
    CRTArrayGPU(std::shared_ptr<const U64CtxChain> cc);
    ~CRTArrayGPU();
    // 不可以拷贝构造，请用copy_from
    CRTArrayGPU(const CRTArrayGPU&) = delete;
    // 不可以拷贝赋值，请用copy_from
    CRTArrayGPU& operator=(const CRTArrayGPU&) = delete;
    // 可以移动构造
    CRTArrayGPU(CRTArrayGPU&&);
    // 不可以移动赋值，请用copy_from
    CRTArrayGPU& operator=(CRTArrayGPU&&) = delete;

    // 从内存的vv64中载入数据
    void set_from_vv64(const vv64& data);
    // 将数据写入内存（vv64格式）
    vv64 export_to_vv64() const;
    // 设为全0
    void set_to_zero();
    // 从另一个对象处拷贝（取代拷贝赋值）
    void copy_from(const CRTArrayGPU& other);

    

    // 一元运算
    void neg_inplace();
    void mul_scalar_inplace(uint64_t scalar);   // 注：scalar是未经encode的
    void mont_encode_inplace();
    void mont_decode_inplace();

    void iw_ntt_inplace();
    void iw_intt_inplace();
    void xy_ntt_inplace();
    void xy_intt_inplace();

    // 二元运算：加、减、乘
    // 乘法仅提供蒙哥马利形式
    // 别名是被允许的
    void add(const CRTArrayGPU& src1, const CRTArrayGPU& src2);
    void sub(const CRTArrayGPU& src1, const CRTArrayGPU& src2);
    void mul_mont(const CRTArrayGPU& src1, const CRTArrayGPU& src2);

    // 比较相等
    bool eq(const CRTArrayGPU& other) const;

    inline std::shared_ptr<const U64CtxChain> get_cc() const { return cc_; }

};