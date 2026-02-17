#pragma once
#include "u64_array.hpp"
#include <optional>
#include <variant>
#include <stdexcept>


/// 面向用户的Zq[i][X,Y,W]/(...)类型
/// 可能表示其
class PlaintextPoly
{
public:
    enum class Device {
        CPU, CUDA
    };
private:
    std::variant<std::monostate, CRTArray, CRTArrayGPU> value_;
    std::shared_ptr<U64CtxChain> cc_;

public:
    // 构造函数：构造内容未初始化的对象
    PlaintextPoly(
        std::shared_ptr<U64CtxChain> cc,
        Device device
    );
    // 自动析构
    ~PlaintextPoly() = default;
    // 自动移动构造
    PlaintextPoly(PlaintextPoly&&) = default;
    // 禁止移动赋值
    PlaintextPoly& operator=(PlaintextPoly&&) = delete;
    // 禁止复制
    PlaintextPoly(const PlaintextPoly&) = delete;
    PlaintextPoly& operator=(const PlaintextPoly&) = delete;


    Device device() const {
        if (std::holds_alternative<CRTArray>(value_))Device::CPU;
        if (std::holds_alternative<CRTArrayGPU>(value_))Device::CUDA;
        throw std::runtime_error("PlaintextPoly::device: 没有检查到device!\n");
    }

    std::shared_ptr<U64CtxChain> cc() const {
        return cc_;
    }

    // 切换设备
    PlaintextPoly to(Device) const;


    bool like(const PlaintextPoly& other) const {
        return (device() == other.device()) && (cc_ == other.cc_);
    }

    // 单目运算
    static void neg(PlaintextPoly& dst, const PlaintextPoly& src);
    static void mont_encode(PlaintextPoly& dst, const PlaintextPoly& src);
    static void mont_decode(PlaintextPoly& dst, const PlaintextPoly& src);

    static void iw_ntt(PlaintextPoly& dst, const PlaintextPoly& src);
    static void iw_intt(PlaintextPoly& dst, const PlaintextPoly& src);
    static void xy_ntt(PlaintextPoly& dst, const PlaintextPoly& src);
    static void xy_intt(PlaintextPoly& dst, const PlaintextPoly& src);

    static void transpose(PlaintextPoly& dst, const PlaintextPoly& src);
    static void conj(PlaintextPoly& dst, const PlaintextPoly& src);
    static void w_inv(PlaintextPoly& dst, const PlaintextPoly& src);

    // 双目运算

    // 加法
    static void add(PlaintextPoly& dst, const PlaintextPoly& src1, const PlaintextPoly& src2);
    // 减法
    static void sub(PlaintextPoly& dst, const PlaintextPoly& src1, const PlaintextPoly& src2);
    // 乘法（普通语义）
    static void mul(PlaintextPoly& dst, const PlaintextPoly& src1, const PlaintextPoly& src2);
    // 乘法（蒙哥马利语义）
    static void mul_mont(PlaintextPoly& dst, const PlaintextPoly& src1, const PlaintextPoly& src2);
    // 标量乘
    static void mul_scalar(PlaintextPoly& dst, const PlaintextPoly& src1, uint64_t src2);

    // setter
    
    // 清零
    void set_zero();

    // 从另一个对象中拷贝数据。禁止跨设备拷贝
    void set_from_other(const PlaintextPoly& other);

    void set_from_fmpz(const fmpz_vector&);

    fmpz_vector to_fmpz_vector() const;
    fmpz_vector to_fmpz_vector_centered() const;

    // 工厂函数：
    static PlaintextPoly zeros(std::shared_ptr<U64CtxChain> cc, Device device);
    static PlaintextPoly uniform(std::shared_ptr<U64CtxChain> cc, Device device);
    static PlaintextPoly dg(std::shared_ptr<U64CtxChain> cc, Device device);
    static PlaintextPoly sk(std::shared_ptr<U64CtxChain> cc, Device device);
    static PlaintextPoly randint(std::shared_ptr<U64CtxChain> cc, Device device, int lb, int ub);

    // for KS
    // PlaintextPoly mod_reduce() const;
    // mod_by_modulo

};