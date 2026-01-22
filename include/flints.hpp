#pragma once
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_vec.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mod_poly.h>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>

// 各类转化函数

std::string fmpz_to_string(fmpz_t x);
void string_to_fmpz(const std::string s, fmpz_t x);

// C++风格的fmpz_t
// 这是一个纯header的class，以便调用者内联它
// 它的工作是帮助调用者管理fmpz_t的生命周期
// 它适合作为临时对象、类成员
// 但是它没有移动语义，因此它不应该作为函数返回值！
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
    // 不能安全移动！
    inline fmpz_scalar(fmpz_scalar&& other) = delete;
    

    // 本函数以string_view为形参，可以传入const char*, std::string等实参，很灵活
    inline fmpz_scalar(std::string_view str) {
        fmpz_init(data_);
        int error = fmpz_set_str(data_, str.data(), 10);
        if (error)
        {
            throw std::invalid_argument("Invalid integer string: " + std::string(str));
        }
    }
    inline fmpz_scalar(long x) {
        // fmpz_init_set_si以mp_limb_signed_t 为第二个形参
        // 我们最好静态断言一下，这是一回事
        static_assert(std::is_same_v<mp_limb_signed_t, long>);
        fmpz_init_set_si(data_, x);
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

    inline fmpz* raw() {
        return data_;
    }
    inline const fmpz* raw() const {
        return data_;
    }

    inline std::string to_string() const {
        char* str = fmpz_get_str(nullptr, 10, data_); // 获取十进制字符串
        std::string result(str);
        flint_free(str); // 释放内存
        return result;
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

    ~fmpz_vector();

    fmpz_vector(const fmpz_vector&); // 拷贝构造

    fmpz_vector(fmpz_vector&&); // 移动构造

    void print() const;

    inline fmpz* raw(){
        return data_;
    }
    inline const fmpz* raw() const {
        return data_;
    }

    inline int len() const { return len_; }

    inline fmpz* at(int idx) {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline fmpz* operator[](int idx) {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline const fmpz* at(int idx) const {
        assert(idx < len_);
        assert(idx >= 0);
        return data_ + idx;
    }

    inline const fmpz* operator[](int idx) const {
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
    long max_abs() const;

};