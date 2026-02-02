#pragma once
#include <cstdint>
#include <iostream>
#include <type_traits>

/*

在本头文件中，我们定义可以在编译时确定模数和尺寸的整数类和array类

*/

template<uint64_t M>
class ModInt {
public:
    using value_type = uint64_t;

private:
    value_type val_;
    static_assert(M>0);
    // 辅助函数：标准化到 [0, M)
    static constexpr value_type normalize(value_type x) noexcept {
        return x % M;
    }

public:
    // 构造函数
    constexpr ModInt() noexcept : val_(0) {}
    
    template<typename T>
    constexpr ModInt(T x) noexcept : val_(normalize(static_cast<value_type>(x % static_cast<T>(M)))) {
        static_assert(std::is_integral_v<T>);
        if constexpr (std::is_signed_v<T>) {
            if (x < 0) {
                val_ = (val_ == 0 ? 0 : M - val_);
            }
        }
    }

    // 拷贝与移动
    constexpr ModInt(const ModInt&) noexcept = default;
    constexpr ModInt(ModInt&&) noexcept = default;
    constexpr ModInt& operator=(const ModInt&) noexcept = default;
    constexpr ModInt& operator=(ModInt&&) noexcept = default;

    // 转换为原始值
    constexpr value_type val() const noexcept { return val_; }

    // 运算符重载：一元
    constexpr ModInt operator+() const noexcept { return *this; }
    constexpr ModInt operator-() const noexcept { return val_ == 0 ? ModInt(0) : ModInt(M - val_); }

    // 自增自减
    constexpr ModInt& operator++() noexcept {
        val_ = (val_ + 1) % M;
        return *this;
    }
    constexpr ModInt operator++(int) noexcept {
        ModInt tmp(*this);
        ++(*this);
        return tmp;
    }
    constexpr ModInt& operator--() noexcept {
        val_ = (val_ + M - 1) % M;
        return *this;
    }
    constexpr ModInt operator--(int) noexcept {
        ModInt tmp(*this);
        --(*this);
        return tmp;
    }

    // 复合赋值运算符
    constexpr ModInt& operator+=(const ModInt& rhs) noexcept {
        val_ = (val_ + rhs.val_) % M;
        return *this;
    }
    constexpr ModInt& operator-=(const ModInt& rhs) noexcept {
        val_ = (val_ + M - rhs.val_) % M;
        return *this;
    }
    constexpr ModInt& operator*=(const ModInt& rhs) noexcept {
        __uint128_t product = (__uint128_t)val_ * rhs.val_;
        val_ = (u64)(product % M);
        return *this;
    }

    // 二元运算符（非成员，定义在类内 friend）
    friend constexpr ModInt operator+(const ModInt& lhs, const ModInt& rhs) noexcept {
        return ModInt(lhs) += rhs;
    }
    friend constexpr ModInt operator-(const ModInt& lhs, const ModInt& rhs) noexcept {
        return ModInt(lhs) -= rhs;
    }
    friend constexpr ModInt operator*(const ModInt& lhs, const ModInt& rhs) noexcept {
        return ModInt(lhs) *= rhs;
    }

    // 比较运算符
    friend constexpr bool operator==(const ModInt& lhs, const ModInt& rhs) noexcept {
        return lhs.val_ == rhs.val_;
    }
    friend constexpr bool operator!=(const ModInt& lhs, const ModInt& rhs) noexcept {
        return !(lhs == rhs);
    }
    friend constexpr bool operator<(const ModInt& lhs, const ModInt& rhs) noexcept {
        return lhs.val_ < rhs.val_;
    }
    friend constexpr bool operator<=(const ModInt& lhs, const ModInt& rhs) noexcept {
        return lhs.val_ <= rhs.val_;
    }
    friend constexpr bool operator>(const ModInt& lhs, const ModInt& rhs) noexcept {
        return lhs.val_ > rhs.val_;
    }
    friend constexpr bool operator>=(const ModInt& lhs, const ModInt& rhs) noexcept {
        return lhs.val_ >= rhs.val_;
    }

    // 幂运算（快速幂）
    constexpr ModInt pow(uint64_t exp) const noexcept {
        ModInt base = *this;
        ModInt result(1);
        while (exp > 0) {
            if (exp & 1) result *= base;
            base *= base;
            exp >>= 1;
        }
        return result;
    }

    // 求乘法逆元（仅当 M 为质数时对任意非零元有效；否则需 gcd(a, M) == 1）
    constexpr ModInt inv() const noexcept {
        return pow(M - 2);
    }

    // 流输出
    friend std::ostream& operator<<(std::ostream& os, const ModInt& x) {
        return os << x.val_;
    }
};

// 为方便推导，可提供字面量或辅助函数（可选）