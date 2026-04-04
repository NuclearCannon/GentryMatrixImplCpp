#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_vec.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/ulong_extras.h>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <ctype.h>
#include <cstdint>

// 各类转化函数

std::string fmpz_to_string(fmpz_t x);
void string_to_fmpz(const std::string& s, fmpz_t x);

// C++风格的fmpz_t
// 这是一个纯header的class，以便调用者内联它
// 它的工作是帮助调用者管理fmpz_t的生命周期
// 它适合作为临时对象、类成员
class fmpz_scalar
{
private:
    fmpz_t data_;
public:
    inline fmpz_scalar() {
        fmpz_init(data_);
    }
    inline ~fmpz_scalar() {
        fmpz_clear(data_);
    }
    inline fmpz_scalar(const fmpz_scalar& other) {
        fmpz_init_set(data_, other.data_);
    }
    inline fmpz_scalar(const fmpz* other) {
        fmpz_init_set(data_, other);
    }
    // 安全的移动构造（使用 swap）
    inline fmpz_scalar(fmpz_scalar&& other) noexcept {
        fmpz_init(data_);           // 初始化为 0
        fmpz_swap(data_, other.data_); // 窃取数据
    }
    

    inline static fmpz_scalar from_si(ssize_t x) {
        static_assert(std::is_same_v<mp_limb_signed_t, ssize_t>);
        fmpz_scalar result;
        fmpz_init_set_si(result.data_, x);
        return result;
    }

    inline static fmpz_scalar from_ui(size_t x) {
        static_assert(std::is_same_v<mp_limb_t, size_t>);
        fmpz_scalar result;
        fmpz_init_set_ui(result.data_, x);
        return result;
    }

    // 拷贝赋值
    inline fmpz_scalar& operator=(const fmpz_scalar& other) {
        if (this != &other) {
            // 先清空自己
            fmpz_clear(data_);
            // 再初始化并复制
            fmpz_init_set(data_, other.data_);
        }
        return *this;
    }

    inline __attribute__((always_inline))
    fmpz* raw() {
        return data_;
    }

    inline __attribute__((always_inline))
    const fmpz* raw() const {
        return data_;
    }

    
};



class fmpz_vector
{
private:
    fmpz* data_;
    int len_;
public:

    fmpz_vector(int len);
    fmpz_vector(const std::vector<std::string>&);
    fmpz_vector(const std::vector<uint64_t>&);

    ~fmpz_vector();

    fmpz_vector(const fmpz_vector&); // 拷贝构造

    fmpz_vector(fmpz_vector&&); // 移动构造

    void print() const;

    inline __attribute__((always_inline)) 
    fmpz* raw(){
        return data_;
    }

    inline __attribute__((always_inline)) 
    const fmpz* raw() const {
        return data_;
    }

    inline __attribute__((always_inline)) 
    int len() const { return len_; }

    inline __attribute__((always_inline))
    fmpz* at(int idx) {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline __attribute__((always_inline))
    fmpz* operator[](int idx) {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline __attribute__((always_inline))
    const fmpz* at(int idx) const {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline __attribute__((always_inline))
    const fmpz* operator[](int idx) const {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }


    std::vector<std::string> export_to_vec_str() const;

    static fmpz_vector zeros(int len);
    static fmpz_vector uniform(int len, const fmpz_t q);
    static fmpz_vector dg(int len);

    // 返回按q模到[-q/2, q/2]的新向量
    fmpz_vector mod_centered(const fmpz_t q) const;

    // 绝对值的最大值
    int64_t max_abs() const;
    std::vector<uint64_t> to_uint64() const;

};

