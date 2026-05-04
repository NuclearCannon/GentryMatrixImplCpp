#include "application.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>

static std::vector<uint64_t> get_mods_from_qrp(const std::vector<std::pair<uint64_t, uint64_t>>& src, int l, int r)
{
    std::vector<uint64_t> result;
    for(int i=l; i<r; i++)result.push_back(src[i].first);
    return result;
}

static std::vector<uint64_t> slice(const std::vector<uint64_t>& src, int l, int r)
{
    std::vector<uint64_t> result;
    for(int i=l; i<r; i++)result.push_back(src[i]);
    return result;
}

bool veccmp(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b)
{
    if(a.size() != b.size())return false;
    for(int i=0; i<a.size(); i++)
    {
        if(a[i]!=b[i])return false;
    }
    return true;
}


long double findMinAbs(long double a, long double b) {
    long double bestValue = a + (-10) * b;  // 初始值设为 k = -10
    long double bestAbs = std::abs(bestValue);
    
    for (int k = -9; k <= 10; ++k) {
        long double current = a + k * b;
        long double currentAbs = std::abs(current);
        
        if (currentAbs < bestAbs) {
            bestAbs = currentAbs;
            bestValue = current;
        }
    }
    
    return bestValue;
}

/*

模数链分配：
C2S: 3层
乘以2\pi/q: 4层。一次pcmul(*2pi)，三次模约简。多出去的模约简是为了除以q。毕竟，2pi/q这么小的数值无法被编码。
sin: 14层
    除以2^10: 1层
    泰勒展开：2层
    除以24：1层
    10次自平方：10层
S2C: 3层
共计20层

*/
class BootstrapKey
{
private:
    int n_, p_;
    const GentryPolyCtx* ctx_;
    std::vector<uint64_t> mods_all_;
    static constexpr int low_ = 1;   // 我们假设输入密文的模数就是mods_all_[:low_] 实际上，low必须等于2
    static constexpr int top_ = 24;   // 我们假设最大密文的模数就是mods_all_[:top_]
    static constexpr int r_ = 10;
    static constexpr long double DELTA = (long double)(1ULL<<40);
    // C2S
    CircledastKey mmkey_top_;
    ConjTransposeKey ctkey_top_;
    std::vector<RotateKey> rtkeys_top_;

    // 分离虚实
    ConjKey cjkey_t3_;
    ConjKey cjkey_sin_;

    // S2C
    CircledastKey mmkey_s2c_;
    ConjTransposeKey ctkey_s2c_;
    std::vector<RotateKey> rtkeys_s2c_;


    // sin
    std::vector<MultKey> mtkeys_bs_;
    MultKey mtkey_c1_, mtkey_c2_;

public:
    BootstrapKey(
        int n,
        int p,
        const SecretKey& sk,
        uint64_t qo, 
        const GentryPolyCtx& ctx,
        const std::vector<uint64_t>& mods
        
    ):
        mmkey_top_(   CircledastKey::gen(sk, qo, ctx, slice(mods, 0, top_))),
        ctkey_top_(ConjTransposeKey::gen(sk, qo, ctx, slice(mods, 0, top_))),
        cjkey_t3_(ConjKey::gen(sk, qo, ctx, slice(mods, 0, top_-4))),
        cjkey_sin_(ConjKey::gen(sk, qo, ctx, slice(mods, 0, top_-18))),
        mmkey_s2c_(   CircledastKey::gen(sk, qo, ctx, slice(mods, 0, top_-18))),
        ctkey_s2c_(ConjTransposeKey::gen(sk, qo, ctx, slice(mods, 0, top_-18))),
        mtkey_c1_ (MultKey::gen(sk, qo, ctx, slice(mods, 0, top_-5))),
        mtkey_c2_ (MultKey::gen(sk, qo, ctx, slice(mods, 0, top_-6)))
    {
        n_ = n;
        p_ = p;
        ctx_ = &ctx;
        mods_all_ = mods;
        // 生成C2S的一系列RotateKey
        for(int l=0; l<p-1; l++)
        {
            rtkeys_top_.push_back(RotateKey::gen(sk, qo, ctx, slice(mods, 0, top_), 0,0,l));
            rtkeys_s2c_.push_back(RotateKey::gen(sk, qo, ctx, slice(mods, 0, top_-18), 0,0,(p-1-l)%(p-1)));
        }
        // 生成bs所需的一系列mtkey。共计10个，第一个的模数链长度是
        for(int i=0; i<10; i++)
        {
            int j = top_-8-i;
            mtkeys_bs_.push_back(MultKey::gen(sk, qo, ctx, slice(mods, 0, j)));
        }
    }

    /// _c2s的输入密文的模数链是top, 缩放因子是delta
    // 模数链长度：top->top-3
    Ciphertext _c2s(const Ciphertext& input, long double delta)
    {
        int n=n_, p=p_;
        std::vector<uint64_t> mods = input.get_moduli();
        int modslen = mods.size();
        assert(modslen == top_);

        
        // 构造XY-DFT辅助矩阵C
        ComplexMatrixGroup C = ComplexMatrixGroup(n, p);    // zero
        {
            
            for(int i=0, i5=1; i<n; i++, i5=i5*5%(4*n))
            {
                // i5 = 5^i
                for(int j=0; j<n; j++)
                {
                    // 计算zeta^{-5^i * j}
                    complex x = std::polar<long double>(1, - M_PI_2 / n * (i5 * j % (4*n)));
                    for(int k=0; k<p-1; k++)
                    {
                        C.at(k, j, i) = x;
                    }
                }
            }
        }
        // 将C化为明文
        Plaintext ptC1 = Plaintext::from_cmat(C, mods, mods[modslen-1]);
        Plaintext ptC2 = Plaintext::from_cmat(C, mods, mods[modslen-2]);
        // 执行XY-DFT
        Ciphertext m02 = mmkey_top_.run_pc(
            ptC1, ctkey_top_.run(
                mmkey_top_.run_cp(
                    input, ptC2, *ctx_
                ), *ctx_
            ), *ctx_
        );// m02的缩放因子是delta^3
        // 执行W-DFT
        Ciphertext m03 = Ciphertext::zeros(n, p, mods); // m03的缩放因子是delta^4
        
        for(int l=0, gamma=1; l<p-1; l++, gamma=gamma*3%p)
        {
            // 构造Rotate l的ksk
            const RotateKey& rkey = rtkeys_top_[l];
            Ciphertext rotated = rkey.run(m02, *ctx_);
            // 构造临时明文
            complex x = std::polar<long double>(1, - (2 * M_PI / ((long double)p)) * (long double)gamma);
            x = (x-1.0L)/((long double)p);
            Plaintext tmppt = Plaintext::from_scalar(n, p, x, mods, mods[modslen-3]);
            
            Ciphertext multed = rotated.mul_pt(tmppt, *ctx_);
            m03.add_(multed, *ctx_);
        }
        return m03.moduli_reduce(slice(mods, modslen-3, modslen));
    }

    Ciphertext _s2c(const Ciphertext& input, long double delta, long double factor)
    {
        int n=n_, p=p_;
        std::vector<uint64_t> mods = input.get_moduli();
        assert(veccmp(mods, slice(mods_all_, 0, top_-18)));
        const GentryPolyCtx& ctx = *ctx_;
        // 构造XY-DFT辅助矩阵C2=C*
        ComplexMatrixGroup C2 = ComplexMatrixGroup(n, p);    // zero
            
        for(int i=0, i5=1; i<n; i++, i5=i5*5%(4*n))
        {
            // i5 = 5^i
            for(int j=0; j<n; j++)
            {
                // 注意在这里存在一个n倍关系
                complex x = std::polar<long double>(n, M_PI_2 / n * (i5 * j % (4*n)));
                for(int k=0; k<p-1; k++)
                {
                    C2.at(k, i, j) = x;
                }
            }
        }
        
        // 将C化为明文
        Plaintext ptC = Plaintext::from_cmat(C2, mods, delta);
        // 执行XY-DFT
        Ciphertext m02 = mmkey_s2c_.run_pc(
            ptC, ctkey_s2c_.run(
                mmkey_s2c_.run_cp(
                    input, ptC, ctx
                ), ctx
            ), ctx
        );// m02的缩放因子是delta^3
        // 执行W-DFT
        Ciphertext m03 = Ciphertext::zeros(n, p, mods); // m03的缩放因子是delta^4
        for(int l=0, gamma=1; l<p-1; l++, gamma=gamma*3%p)
        {
            // 构造Rotate l的ksk
            const RotateKey& rkey = rtkeys_s2c_[l];
            Ciphertext rotated = rkey.run(m02, ctx);
            // 构造临时明文
            complex x = std::polar<long double>(factor, (2 * M_PI / ((long double)p)) * (long double)gamma);
            Plaintext tmppt = Plaintext::from_scalar(n, p, x, mods, delta);
            Ciphertext multed = rotated.mul_pt(tmppt, ctx);
            m03.add_(multed, ctx);
        }
        return m03.moduli_reduce(slice(mods_all_, top_-21, top_-18));

    }


    Ciphertext _sin(const Ciphertext& input, long double delta)
    {
        const int input_layer = top_-4;
        assert(veccmp(slice(mods_all_, 0, input_layer), input.get_moduli()));
        // 准备各个参数
        int n = n_;
        int p = p_;
        const GentryPolyCtx& ctx = *ctx_;
        

        // 准备一个明文
        std::vector<std::vector<uint64_t>> chains;
        // chains[i]是去掉倒数的i个元素后的模数链
        for(int i=input_layer; i>0; i--)
            chains.push_back(slice(mods_all_, 0, i));

        std::vector<std::vector<uint64_t>> to_remove;
        // to_remove[i]只包括倒数第i个模数
        for(int i=input_layer; i>0; i--)
            to_remove.push_back(slice(mods_all_, i-1, i));
        
        // 除以2^r，消耗一层
        constexpr int r = 10;
        constexpr long double exp2r = (long double)(1ULL<<r);
        constexpr long double exp2r_inv = 1.0/exp2r;
        Plaintext pt_2r_inv = Plaintext::from_scalar(n, p, exp2r_inv, chains[0], to_remove[0][0]);
        Ciphertext x_0 = input.mul_pt(pt_2r_inv, ctx).moduli_reduce(to_remove[0]);
        Plaintext one_0 = Plaintext::from_scalar(n, p, 1, chains[1], to_remove[1][0]);
        Plaintext one_1 = Plaintext::from_scalar(n, p, 1, chains[2], to_remove[2][0]);
        Ciphertext x2_1 = mtkey_c1_.run(x_0,  x_0,  ctx).moduli_reduce(to_remove[1]);
        
        Ciphertext x_1 = x_0.mul_pt(one_0, ctx).moduli_reduce(to_remove[1]);
        Ciphertext x4_2 = mtkey_c2_.run(x2_1, x2_1, ctx).moduli_reduce(to_remove[2]);
        Ciphertext x3_2 = mtkey_c2_.run(x2_1, x_1, ctx).moduli_reduce(to_remove[2]).mul_int(4);
        Ciphertext x2_2 = x2_1.mul_pt(one_1, ctx).moduli_reduce(to_remove[2]).mul_int(12);
        Ciphertext x1_2 = x_1.mul_pt(one_1, ctx).moduli_reduce(to_remove[2]).mul_int(24);
        
        x4_2.add_(x3_2, ctx);
        x4_2.add_(x2_2, ctx);
        x4_2.add_(x1_2, ctx);
        
        
        Ciphertext y0 = x4_2
            .add_pt(Plaintext::from_scalar(n, p, 24.0, chains[3], delta), ctx)
            .mul_pt(Plaintext::from_scalar(n, p, 1.0/24.0, chains[3], to_remove[3][0]), ctx)
            .moduli_reduce(to_remove[3]);


        // 对y0进行10次连续平方
        // y0对应的模数链是chains[4]，也即[0:19]
        Ciphertext y1 = mtkeys_bs_[0].run(y0, y0, ctx).moduli_reduce(to_remove[4]);
        Ciphertext y2 = mtkeys_bs_[1].run(y1, y1, ctx).moduli_reduce(to_remove[5]);
        Ciphertext y3 = mtkeys_bs_[2].run(y2, y2, ctx).moduli_reduce(to_remove[6]);
        Ciphertext y4 = mtkeys_bs_[3].run(y3, y3, ctx).moduli_reduce(to_remove[7]);
        Ciphertext y5 = mtkeys_bs_[4].run(y4, y4, ctx).moduli_reduce(to_remove[8]);
        Ciphertext y6 = mtkeys_bs_[5].run(y5, y5, ctx).moduli_reduce(to_remove[9]);
        Ciphertext y7 = mtkeys_bs_[6].run(y6, y6, ctx).moduli_reduce(to_remove[10]);
        Ciphertext y8 = mtkeys_bs_[7].run(y7, y7, ctx).moduli_reduce(to_remove[11]);
        Ciphertext y9 = mtkeys_bs_[8].run(y8, y8, ctx).moduli_reduce(to_remove[12]);
        Ciphertext y10 =mtkeys_bs_[9].run(y9, y9, ctx).moduli_reduce(to_remove[13]);


        // 现在y10的虚部应该就是sin的取值了
        return y10;
    }

    Ciphertext main(const Ciphertext& input)
    {
        // 简单模数提高
        assert(veccmp(input.get_moduli(), slice(mods_all_, 0, low_)));
        Ciphertext t1 = input.naive_moduli_extend(slice(mods_all_, low_, top_));
        assert(veccmp(t1.get_moduli(), slice(mods_all_, 0, top_)));
        // C2S
        Ciphertext t2 = this->_c2s(t1, DELTA);  // 它的模数链为top-3=27级
        assert(veccmp(t2.get_moduli(), slice(mods_all_, 0, top_-3)));
        // t2的槽中元素（缩放后后）的周期是mods[0]/delta
        // 我们希望它的周期是2pi，因此这里需要乘以pi
        // 为什么是pi而不是2pi？因为后面的虚实分离有乘2的副作用
        // 我们希望去除mods[0]，但是模约简却使用mods[26]，因此我们还需乘以mods[26]/mods[0]
        long double scalar = M_PI;
        scalar *= mods_all_[top_-4];
        scalar /= mods_all_[0];

        Plaintext pt_pi = Plaintext::from_scalar(n_, p_, scalar, t2.get_moduli(), DELTA); // 嘿，或许可以省一层下来
        Ciphertext t3 = t2.mul_pt(pt_pi, *ctx_).moduli_reduce(slice(mods_all_, top_-4, top_-3));
        assert(veccmp(t3.get_moduli(), slice(mods_all_, 0, top_-4)));
        // 提取实部虚部
        Ciphertext t3_conj = cjkey_t3_.run(t3, *ctx_);
        Ciphertext real = t3.add(t3_conj).mul_i();
        Ciphertext imag = t3.sub(t3_conj);
        assert(veccmp(real.get_moduli(), slice(mods_all_, 0, top_-4)));
        assert(veccmp(imag.get_moduli(), slice(mods_all_, 0, top_-4)));
        Ciphertext sin_real = this->_sin(real, DELTA);
        Ciphertext sin_imag = this->_sin(imag, DELTA);
        assert(veccmp(sin_real.get_moduli(), slice(mods_all_, 0, top_-18)));
        assert(veccmp(sin_imag.get_moduli(), slice(mods_all_, 0, top_-18)));
        // 提取他们的虚部，这会造成绝对值乘2的副作用
        Ciphertext res_real = sin_real.sub(cjkey_sin_.run(sin_real, *ctx_)).mul_i().neg();
        Ciphertext res_imag = sin_imag.sub(cjkey_sin_.run(sin_imag, *ctx_));
        // 合并
        Ciphertext res_of_sin = res_real.add(res_imag);
        assert(veccmp(res_of_sin.get_moduli(), slice(mods_all_, 0, top_-18)));
        
        // 乘以q/4pi。不是说q/2pi吗？因为我们刚刚造成了副作用，现在要补偿回去。
        // 但是我们什么都不用做。q部分一开始就没乘上去，而剩下来的部分可以放到s2c中作为它的一部分。
        Ciphertext res = this->_s2c(res_of_sin, DELTA, 1.0091/(4 * M_PI));
        assert(veccmp(res.get_moduli(), slice(mods_all_, 0, top_-21)));
        
        return res;
    }

    
};



void test_bsk()
{
    printf("进入test_bsk\n");
    // 准备各个参数
    int n = 8;
    int p = 5;

    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        //268435456
        {1099512345601, 7},
        {1099512693761, 3},
        {1099513216001, 3}, // 1329232072329048371585914525299471361
        {1099513999361, 6}, // -6044647636049595256012800
        {1099514434561, 31},
        {1099514782721, 3},
        {1099515566081, 6},
        {1099516001281, 7},
        {1099516871681, 3},
        {1099518873601, 13},
        {1099519395841, 7},
        {1099519482881, 3},
        {1099520005121, 6},
        {1099520788481, 3},
        {1099521832961, 3},
        {1099522007041, 7},
        {1099523834881, 11},
        {1099523921921, 3},
        {1099524705281, 3},
        {1099525662721, 13},
        {1099526794241, 6},
        {1099526968321, 13},
        {1099527490561, 7},
        {1099528012801, 37},
        {1099528099841, 3},
        {1099532712961, 13},
        {1099532974081, 14},
        {1099533583361, 6},
        {1099533844481, 6},
        {1099534540801, 13},
        {1099535063041, 22},
        {1099536107521, 7},
        {1099536368641, 19},
        {1099538196481, 19},
        {1099538718721, 47},
        {1099545072641, 6},
        {1099545856001, 3},
        {1099546813441, 7},
        {1099548728321, 3},
        {1099549250561, 7},
        {1099550556161, 11},
        {1099550991361, 31},
        {1099551078401, 3},
        {1099552558081, 7},
        {1099553080321, 7},
        {1099554211841, 6},
        {1099557345281, 3},
        {1099557780481, 33},
        {1099558650881, 3},
        {1099559695361, 3},
        {qo, qor}
    };
    long double delta = (long double)(1ULL<<40);    // 缩放因子和每一级模数都是2^40级

    GentryPolyCtx ctx(n, p, qrp);
    std::vector<uint64_t> mods = get_mods_from_qrp(qrp, 0, 50);

    // 准备一个明文
    printf("准备一个明文\n");

    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(0.01, n, p);
    // Plaintext pt1 = Plaintext::from_cmat(mat1, slice(mods, 0, 3), delta);
    Plaintext pt1 = Plaintext::from_cmat(mat1, slice(mods, 0, 1), delta);
    // 准备一个私钥
    printf("准备一个私钥\n");
    SecretKey sk(n, p);
    printf("准备自举所需的公钥\n");
    auto t_bsk_start = std::chrono::steady_clock::now();
    BootstrapKey bsk(n, p, sk, qo, ctx, mods);
    auto t_bsk_end = std::chrono::steady_clock::now();
    double bsk_construct_time = std::chrono::duration<double>(t_bsk_end - t_bsk_start).count();
    printf("bsk构造耗时: %.6f 秒\n", bsk_construct_time);
    // 加密成密文
    printf("加密成密文\n");
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    printf("开始自举\n");
    auto t_main_start = std::chrono::steady_clock::now();
    Ciphertext ct2 = bsk.main(ct1);
    auto t_main_end = std::chrono::steady_clock::now();
    double bsk_main_time = std::chrono::duration<double>(t_main_end - t_main_start).count();
    printf("bsk.main耗时: %.6f 秒\n", bsk_main_time);
    printf("自举结束，检查取值\n");

    Plaintext decrypted = ct2.decrypt(sk, ctx);
    
    // 这里检查最终结果
    ComplexMatrixGroup mat2 = decrypted.to_cmat(delta);
    ComplexMatrixGroup mat0 = mat1.decode();
    long double max_diff = 0;
    long double sum_rate = 0;
    for(int w=0; w<p-1; w++)for(int x=0; x<n; x++)for(int y=0; y<n; y++)
    {
        complex got = mat2.at(w, x, y);
        complex expected = mat1.at(w, x, y);
        // complex diff1 = got / expected;
        // printf("[final %d %d %d] got(%+.6Lf,%+.6Lf), exp(%+.6Lf,%+.6Lf) diff1(%+.6Lf,%+.6Lf) diff2(%+.6Lf,%+.6Lf)\n",
        //     w, x, y, 
        //     got.real(), got.imag(), 
        //     expected.real(),  expected.imag(), 
        //     diff1.real(),  diff1.imag(), 
        //     diff2.real(), diff2.imag()
        // );
        long double diff_abs = std::abs(got - expected);
        long double rate = std::abs(got / expected);
        if(diff_abs > max_diff)max_diff = diff_abs;
        sum_rate += rate;
    }
    sum_rate /= n*n*(p-1);
    std::cout << "最大绝对差" << max_diff << std::endl;
    std::cout << "平均绝对值比" << sum_rate << std::endl;
    
}