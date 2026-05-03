#include "application.hpp"
#include <iostream>

void exp_bootstrapping()
{
    // 准备各个参数
    int n = 8;
    int p = 5;
    // 质数链
    vec64 mods = {70368747120641, 70368747294721, 70368748426241};
    // 原根链
    vec64 roots = {6, 11, 6};
    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {70368747120641, 6},
        {70368747294721, 11},
        {70368748426241, 6},
        {qo, qor}
    };
    double delta = 100000;

    GentryPolyCtx ctx(n, p, qrp);

    // 准备一个明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(100, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    // 准备一个私钥
    SecretKey sk(n, p);
    CircledastKey mmkey = CircledastKey::gen(sk, qo, ctx, mods);
    ConjTransposeKey ctkey = ConjTransposeKey::gen(sk, qo, ctx, mods);
    // 加密成密文
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);

    bool test_c2s = true;
    bool test_s2c = true;
    bool test_naive_extend = false;
    // bool 
    // C2S测试
    if (test_c2s)
    {
        // 构造XY-DFT辅助矩阵C
        ComplexMatrixGroup C = ComplexMatrixGroup(n, p);    // zero
        {
            
            for(int i=0, i5=1; i<n; i++, i5=i5*5%(4*n))
            {
                // i5 = 5^i
                
                for(int j=0; j<n; j++)
                {
                    // 计算zeta^{-5^i * j}
                    complex x = std::polar<double>(1, - M_PI_2 / n * (i5 * j % (4*n)));
                    for(int k=0; k<p-1; k++)
                    {
                        C.at(k, j, i) = x;
                    }
                }
            }
        }
        // 将C化为明文
        Plaintext ptC = Plaintext::from_cmat(C, mods, delta);
        // 执行XY-DFT
        Ciphertext m02 = mmkey.run_pc(
            ptC, ctkey.run(
                mmkey.run_cp(
                    ct1, ptC, ctx
                ), ctx
            ), ctx
        );// m02的缩放因子是delta^3
        // 执行W-DFT
        Ciphertext m03 = Ciphertext::zeros(n, p, mods); // m03的缩放因子是delta^4
        
        {

            for(int l=0, gamma=1; l<p-1; l++, gamma=gamma*3%p)
            {
                // 构造Rotate l的ksk
                RotateKey rkey = RotateKey::gen(sk, qo, ctx, mods, 0,0,l);
                Ciphertext rotated = rkey.run(m02, ctx);
                // 构造临时明文
                complex x = std::polar<double>(1, - (2 * M_PI / ((double)p)) * (double)gamma);
                x = (x-1.0)/((double)p);
                Plaintext tmppt = Plaintext::from_scalar(n, p, x, mods, delta);
                
                Ciphertext multed = rotated.mul_pt(tmppt, ctx);
                m03.add_(multed, ctx);
            }
        }
        Plaintext final_pt = m03.decrypt(sk, ctx);
        ComplexMatrixGroup result = final_pt.to_cmat(delta*delta*delta*delta);
        ComplexMatrixGroup encoded = mat1.encode();
        
        double max_diff = 0;
        for(int w=0;w<p-1;w++)
        {
            for(int x=0; x<n; x++)
            {
                for(int y=0; y<n; y++)
                {
                    double new_diff = std::abs(result.at(w, x, y)-encoded.at(w, x, y));
                    if(new_diff>max_diff)max_diff=new_diff;

                }
            }
        }
        std::cout << "exp_bootstrapping: BS(C2S): max_diff:" << max_diff << std::endl;
    }


    // S2C测试
    if (test_s2c)
    {
        // 构造XY-DFT辅助矩阵C2=C*
        ComplexMatrixGroup C2 = ComplexMatrixGroup(n, p);    // zero
        {
            
            for(int i=0, i5=1; i<n; i++, i5=i5*5%(4*n))
            {
                // i5 = 5^i
                
                for(int j=0; j<n; j++)
                {
                    // 注意在这里存在一个n倍关系
                    complex x = std::polar<double>(n, M_PI_2 / n * (i5 * j % (4*n)));
                    for(int k=0; k<p-1; k++)
                    {
                        C2.at(k, i, j) = x;
                    }
                }
            }
        }
        // 将C化为明文
        Plaintext ptC = Plaintext::from_cmat(C2, mods, delta);
        // 执行XY-DFT
        Ciphertext m02 = mmkey.run_pc(
            ptC, ctkey.run(
                mmkey.run_cp(
                    ct1, ptC, ctx
                ), ctx
            ), ctx
        );// m02的缩放因子是delta^3
        // 执行W-DFT
        Ciphertext m03 = Ciphertext::zeros(n, p, mods); // m03的缩放因子是delta^4
        
        {

            for(int l=0, gamma=1; l<p-1; l++, gamma=gamma*3%p)
            {
                // 构造Rotate l的ksk
                RotateKey rkey = RotateKey::gen(sk, qo, ctx, mods, 0,0,(p-1-l)%(p-1));
                Ciphertext rotated = rkey.run(m02, ctx);
                // 构造临时明文
                complex x = std::polar<double>(1, (2 * M_PI / ((double)p)) * (double)gamma);
                Plaintext tmppt = Plaintext::from_scalar(n, p, x, mods, delta);
                Ciphertext multed = rotated.mul_pt(tmppt, ctx);
                m03.add_(multed, ctx);
            }
        }
        Plaintext final_pt = m03.decrypt(sk, ctx);
        ComplexMatrixGroup result = final_pt.to_cmat(delta*delta*delta*delta);
        ComplexMatrixGroup decoded = mat1.decode();
        
        double max_diff = 0;
        for(int w=0;w<p-1;w++)
        {
            for(int x=0; x<n; x++)
            {
                for(int y=0; y<n; y++)
                {
                    double new_diff = std::abs(result.at(w, x, y)-decoded.at(w, x, y));
                    if(new_diff>max_diff)max_diff=new_diff;
                }
            }
        }
        std::cout << "exp_bootstrapping: BS(S2C): max_diff:" << max_diff << std::endl;
    }
    

    if (test_naive_extend)
    {
        // 构造一个明文
        ComplexMatrixGroup cmg_coeffs = ComplexMatrixGroup::random(100, n, p);
        Plaintext pt = Plaintext::_from_cmat_without_encoding(cmg_coeffs, mods, delta);
        Ciphertext ct = Ciphertext::encrypt(pt, sk, ctx);
        Ciphertext ct2 = ct.naive_moduli_extend({qo});
        Plaintext pt2 = ct2.decrypt(sk, ctx);
        ComplexMatrixGroup cmg_coeffs_2 = pt2._to_cmat_without_decoding(delta);
        // 检查差异
        printf("test_naive_extend diffs\n");
        double q0 = mods[0], q1 = mods[1], q2=mods[2];
        double qprod = q0*q1*q2; 
        for(int w=0; w<p-1; w++)
        {
            for(int x=0; x<n; x++)
            {
                for(int y=0; y<n; y++)
                {
                    complex diff = cmg_coeffs_2.at(w, x, y) -  cmg_coeffs.at(w, x, y);
                    // diff /= qprod;
                    diff = diff * delta / qprod;
                    printf("[%d %d %d] %.3lf %.3lf\n", w, x, y, diff.real(), diff.imag());
                    printf("[%d %d %d] %.3lf %.3lf\n", w, x, y, diff.real(), diff.imag());
                }
            }
        }
    }
    
}


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

static void check(Ciphertext& ct, const SecretKey& sk, const GentryPolyCtx& ctx, const char* name = nullptr)
{
    if(name==nullptr)name="NoName";
    double r = ct.check_abs(sk, ctx, 1.0);
    r = std::log2(r);
    printf("%s log(abs)=%.2lf\n",name, r);
}

void exp_sin()
{
    // 准备各个参数
    int n = 8;
    int p = 5;

    uint64_t qo = 576460752303421441;
    uint64_t qor = 19;
    std::vector<std::pair<uint64_t, uint64_t>> qrp = {
        {1099512345601, 7},
        {1099512693761, 3},
        {1099513216001, 3},
        {1099513999361, 6},
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
    double delta = double(1ULL<<40);    // 缩放因子和每一级模数都是2^40级

    GentryPolyCtx ctx(n, p, qrp);

    // 准备一个明文
    auto mods20 = get_mods_from_qrp(qrp, 0, 20);
    auto mods19 = get_mods_from_qrp(qrp, 0, 19);
    auto mods18 = get_mods_from_qrp(qrp, 0, 18);
    auto mods17 = get_mods_from_qrp(qrp, 0, 17);
    auto mods16 = get_mods_from_qrp(qrp, 0, 16);
    auto mods15 = get_mods_from_qrp(qrp, 0, 15);
    auto mods14 = get_mods_from_qrp(qrp, 0, 14);
    auto mods13 = get_mods_from_qrp(qrp, 0, 13);
    auto mods12 = get_mods_from_qrp(qrp, 0, 12);
    auto mods11 = get_mods_from_qrp(qrp, 0, 11);
    auto mods10 = get_mods_from_qrp(qrp, 0, 10);
    auto mods9  = get_mods_from_qrp(qrp, 0, 9);
    auto mods8  = get_mods_from_qrp(qrp, 0, 8);
    auto mods7  = get_mods_from_qrp(qrp, 0, 7);
    auto mods6  = get_mods_from_qrp(qrp, 0, 6);

    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random_imag(40*M_PI, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods20, delta);
    // 准备一个私钥
    SecretKey sk(n, p);
    // 加密成密文
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    check(ct1, sk, ctx, "ct1");
    // 现在ct1是待sin密文，缩放因子为delta
    // 除以2^r
    int r = 10;
    double exp2r = double(1ULL<<r);
    double exp2r_inv = 1/exp2r;
    std::cout << "exp2r_inv=" << exp2r_inv << std::endl;
    Plaintext pt_2r_inv = Plaintext::from_scalar(n, p, exp2r_inv, mods20, delta);
    Ciphertext x_0_0 = ct1.mul_pt(pt_2r_inv, ctx);
    check(x_0_0, sk, ctx, "x_0_0");
    Ciphertext x_0 = x_0_0.moduli_reduce(get_mods_from_qrp(qrp, 19, 20));    // 去除第19号模数
    check(x_0, sk, ctx, "x_0");
    Plaintext one_0 = Plaintext::from_scalar(n, p, 1, mods19, delta);
    Plaintext one_1 = Plaintext::from_scalar(n, p, 1, mods18, delta);
    MultKey mult_key_19 = MultKey::gen(sk, qo, ctx, mods19);
    MultKey mult_key_18 = MultKey::gen(sk, qo, ctx, mods18);
    Ciphertext x2_1 = mult_key_19.run(x_0,  x_0,  ctx).moduli_reduce(get_mods_from_qrp(qrp, 18, 19));
    check(x2_1, sk, ctx, "x2_1");
    // return;
    Ciphertext x_1 = x_0.mul_pt(one_0, ctx).moduli_reduce(get_mods_from_qrp(qrp, 18, 19));
    Ciphertext x4_2 = mult_key_18.run(x2_1, x2_1, ctx).moduli_reduce(get_mods_from_qrp(qrp, 17, 18));
    Ciphertext x3_2 = mult_key_18.run(x2_1, x_1, ctx).moduli_reduce(get_mods_from_qrp(qrp, 17, 18)).mul_int(4);
    Ciphertext x2_2 = x2_1.mul_pt(one_1, ctx).moduli_reduce(get_mods_from_qrp(qrp, 17, 18)).mul_int(12);
    Ciphertext x1_2 = x_1.mul_pt(one_1, ctx).moduli_reduce(get_mods_from_qrp(qrp, 17, 18)).mul_int(24);
    x4_2.add_(x3_2, ctx);
    x4_2.add_(x2_2, ctx);
    x4_2.add_(x1_2, ctx);
    Plaintext one_div_24 = Plaintext::from_scalar(n, p, 1.0/24.0, mods17, delta);
    Plaintext one_16 = Plaintext::from_scalar(n, p, 1.0, mods16, delta);
    Ciphertext y0 = x4_2.mul_pt(one_div_24, ctx).moduli_reduce(get_mods_from_qrp(qrp, 16, 17)).add_pt(one_16, ctx);
    // 对y0进行10次连续平方
    MultKey mult_key_16 = MultKey::gen(sk, qo, ctx, mods16);
    MultKey mult_key_15 = MultKey::gen(sk, qo, ctx, mods15);
    MultKey mult_key_14 = MultKey::gen(sk, qo, ctx, mods14);
    MultKey mult_key_13 = MultKey::gen(sk, qo, ctx, mods13);
    MultKey mult_key_12 = MultKey::gen(sk, qo, ctx, mods12);
    MultKey mult_key_11 = MultKey::gen(sk, qo, ctx, mods11);
    MultKey mult_key_10 = MultKey::gen(sk, qo, ctx, mods10);
    MultKey mult_key_9 = MultKey::gen(sk, qo, ctx, mods9);
    MultKey mult_key_8 = MultKey::gen(sk, qo, ctx, mods8);
    MultKey mult_key_7 = MultKey::gen(sk, qo, ctx, mods7);

    Ciphertext y1 = mult_key_16.run(y0, y0, ctx).moduli_reduce(get_mods_from_qrp(qrp, 15, 16));
    Ciphertext y2 = mult_key_15.run(y1, y1, ctx).moduli_reduce(get_mods_from_qrp(qrp, 14, 15));
    Ciphertext y3 = mult_key_14.run(y2, y2, ctx).moduli_reduce(get_mods_from_qrp(qrp, 13, 14));
    Ciphertext y4 = mult_key_13.run(y3, y3, ctx).moduli_reduce(get_mods_from_qrp(qrp, 12, 13));
    Ciphertext y5 = mult_key_12.run(y4, y4, ctx).moduli_reduce(get_mods_from_qrp(qrp, 11, 12));
    Ciphertext y6 = mult_key_11.run(y5, y5, ctx).moduli_reduce(get_mods_from_qrp(qrp, 10, 11));
    Ciphertext y7 = mult_key_10.run(y6, y6, ctx).moduli_reduce(get_mods_from_qrp(qrp, 9, 10));
    Ciphertext y8 = mult_key_9.run(y7, y7, ctx).moduli_reduce(get_mods_from_qrp(qrp, 8, 9));
    Ciphertext y9 = mult_key_8.run(y8, y8, ctx).moduli_reduce(get_mods_from_qrp(qrp, 7, 8));
    Ciphertext y10 =mult_key_7.run(y9, y9, ctx).moduli_reduce(get_mods_from_qrp(qrp, 6, 7));


    // 现在y10的虚部应该就是sin的取值了
    Plaintext decrypted = y10.decrypt(sk, ctx);
    ComplexMatrixGroup mat2 = decrypted.to_cmat(delta);
    // 比较结果
    for(int w=0; w<p-1; w++)
    {
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)
            {
                double got = mat2.at(w, x, y).imag();
                double expected = std::sin(mat1.at(w, x, y).imag());
                double diff = got - expected;
                printf("%d %d %d %.6lf %.6lf %.6lf\n", w, x, y, got, expected, diff);

            }
        }
    }

    
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
    static constexpr int low_ = 3;   // 我们假设输入密文的模数就是mods_all_[:low_] 实际上，low必须等于3
    static constexpr int top_ = 30;   // 我们假设最大密文的模数就是mods_all_[:top_]
    static constexpr int hig_ = 10;   // 我们假设输出密文的模数就是mods_all_[:hig_]

    // C2S
    CircledastKey mmkey_top_;
    ConjTransposeKey ctkey_top_;
    std::vector<RotateKey> rtkeys_top_;

    // S2C
    CircledastKey mmkey_s2c_;
    ConjTransposeKey ctkey_s2c_;
    std::vector<RotateKey> rtkeys_s2c_;


    // BS
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
        mmkey_top_(   CircledastKey::gen(sk, qo, ctx, slice(mods, 0, 30))),
        ctkey_top_(ConjTransposeKey::gen(sk, qo, ctx, slice(mods, 0, 30))),
        mmkey_s2c_(   CircledastKey::gen(sk, qo, ctx, slice(mods, 0, 13))),
        ctkey_s2c_(ConjTransposeKey::gen(sk, qo, ctx, slice(mods, 0, 13))),
        mtkey_c1_ (MultKey::gen(sk, qo, ctx, slice(mods, 0, 22))),
        mtkey_c2_ (MultKey::gen(sk, qo, ctx, slice(mods, 0, 21)))
    {
        n_ = n;
        p_ = p;
        ctx_ = &ctx;
        mods_all_ = mods;
        // 生成C2S的一系列RotateKey
        for(int l=0; l<p-1; l++)
        {
            rtkeys_top_.push_back(RotateKey::gen(sk, qo, ctx, slice(mods, 0, 30), 0,0,l));
            rtkeys_s2c_.push_back(RotateKey::gen(sk, qo, ctx, slice(mods, 0, 13), 0,0,l));
        }
        // 生成bs所需的一系列mtkey。共计10个，第一个的模数链长度是
        for(int i=0; i<10; i++)
        {
            int j = 19-i;
            mtkeys_bs_.push_back(MultKey::gen(sk, qo, ctx, slice(mods, 0, j)));
        }
    }

    /// _c2s的输入密文的模数链是top, 缩放因子是delta
    // 输出密文的模数链仍然是top，缩放因子是delta^4
    Ciphertext _c2s(const Ciphertext& input, double delta)
    {
        int n=n_, p=p_;
        std::vector<uint64_t> mods = input.get_moduli();
        assert(veccmp(mods, slice(mods_all_, 0, top_)));

        
        // 构造XY-DFT辅助矩阵C
        ComplexMatrixGroup C = ComplexMatrixGroup(n, p);    // zero
        {
            
            for(int i=0, i5=1; i<n; i++, i5=i5*5%(4*n))
            {
                // i5 = 5^i
                for(int j=0; j<n; j++)
                {
                    // 计算zeta^{-5^i * j}
                    complex x = std::polar<double>(1, - M_PI_2 / n * (i5 * j % (4*n)));
                    for(int k=0; k<p-1; k++)
                    {
                        C.at(k, j, i) = x;
                    }
                }
            }
        }
        // 将C化为明文
        Plaintext ptC = Plaintext::from_cmat(C, mods, delta);
        // 执行XY-DFT
        Ciphertext m02 = mmkey_top_.run_pc(
            ptC, ctkey_top_.run(
                mmkey_top_.run_cp(
                    input, ptC, *ctx_
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
            complex x = std::polar<double>(1, - (2 * M_PI / ((double)p)) * (double)gamma);
            x = (x-1.0)/((double)p);
            Plaintext tmppt = Plaintext::from_scalar(n, p, x, mods, delta);
            
            Ciphertext multed = rotated.mul_pt(tmppt, *ctx_);
            m03.add_(multed, *ctx_);
        }
        return m03;
    }

    Ciphertext _s2c(const Ciphertext& input, double delta)
    {
        int n=n_, p=p_;
        std::vector<uint64_t> mods = input.get_moduli();
        assert(veccmp(mods, slice(mods_all_, 0, hig_+3)));
        const GentryPolyCtx& ctx = *ctx_;
        // 构造XY-DFT辅助矩阵C2=C*
        ComplexMatrixGroup C2 = ComplexMatrixGroup(n, p);    // zero
        {
            
            for(int i=0, i5=1; i<n; i++, i5=i5*5%(4*n))
            {
                // i5 = 5^i
                
                for(int j=0; j<n; j++)
                {
                    // 注意在这里存在一个n倍关系
                    complex x = std::polar<double>(n, M_PI_2 / n * (i5 * j % (4*n)));
                    for(int k=0; k<p-1; k++)
                    {
                        C2.at(k, i, j) = x;
                    }
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
            complex x = std::polar<double>(1, (2 * M_PI / ((double)p)) * (double)gamma);
            Plaintext tmppt = Plaintext::from_scalar(n, p, x, mods, delta);
            Ciphertext multed = rotated.mul_pt(tmppt, ctx);
            m03.add_(multed, ctx);
        }
        return m03.moduli_reduce(slice(mods_all_, hig_, hig_+3));

    }


    Ciphertext _sin(const Ciphertext& input, double delta)
    {
        // 输入密文的模数链应该是30-3-4=23
        assert(veccmp(slice(mods_all_, 0, 23), input.get_moduli()));
        // 准备各个参数
        int n = n_;
        int p = p_;
        const GentryPolyCtx& ctx = *ctx_;
        

        // 准备一个明文
        std::vector<std::vector<uint64_t>> chains;
        // chains[i]是去掉倒数的i个元素后的模数链
        for(int i=23; i>0; i--)
            chains.push_back(slice(mods_all_, 0, i));

        std::vector<std::vector<uint64_t>> to_remove;
        // to_remove[i]只包括倒数第i个模数
        for(int i=23; i>0; i--)
            to_remove.push_back(slice(mods_all_, i-1, i));
        
        // 除以2^r，消耗一层
        int r = 10;
        double exp2r = double(1ULL<<r);
        double exp2r_inv = 1/exp2r;
        std::cout << "exp2r_inv=" << exp2r_inv << std::endl;
        Plaintext pt_2r_inv = Plaintext::from_scalar(n, p, exp2r_inv, chains[0], delta);
        Ciphertext x_0 = input.mul_pt(pt_2r_inv, ctx).moduli_reduce(to_remove[0]);
        Plaintext one_0 = Plaintext::from_scalar(n, p, 1, chains[1], delta);
        Plaintext one_1 = Plaintext::from_scalar(n, p, 1, chains[2], delta);
        Ciphertext x2_1 = mtkey_c1_.run(x_0,  x_0,  ctx).moduli_reduce(to_remove[1]);
        
        // return;
        Ciphertext x_1 = x_0.mul_pt(one_0, ctx).moduli_reduce(to_remove[1]);
        Ciphertext x4_2 = mtkey_c2_.run(x2_1, x2_1, ctx).moduli_reduce(to_remove[2]);
        Ciphertext x3_2 = mtkey_c2_.run(x2_1, x_1, ctx).moduli_reduce(to_remove[2]).mul_int(4);
        Ciphertext x2_2 = x2_1.mul_pt(one_1, ctx).moduli_reduce(to_remove[2]).mul_int(12);
        Ciphertext x1_2 = x_1.mul_pt(one_1, ctx).moduli_reduce(to_remove[2]).mul_int(24);
        x4_2.add_(x3_2, ctx);
        x4_2.add_(x2_2, ctx);
        x4_2.add_(x1_2, ctx);
        Plaintext one_div_24 = Plaintext::from_scalar(n, p, 1.0/24.0, chains[3], delta);
        Plaintext one_16 = Plaintext::from_scalar(n, p, 1.0, chains[4], delta);
        Ciphertext y0 = x4_2.mul_pt(one_div_24, ctx).moduli_reduce(to_remove[3]).add_pt(one_16, ctx);
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
        Ciphertext t1 = input.naive_moduli_extend(slice(mods_all_, low_+1, top_));
        return t1;
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
        {1099512345601, 7},
        {1099512693761, 3},
        {1099513216001, 3},
        {1099513999361, 6},
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
    double delta = double(1ULL<<40);    // 缩放因子和每一级模数都是2^40级

    GentryPolyCtx ctx(n, p, qrp);
    std::vector<uint64_t> mods = get_mods_from_qrp(qrp, 0, 50);

    // 准备一个明文
    printf("准备一个明文\n");

    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(100, n, p);
    // Plaintext pt1 = Plaintext::from_cmat(mat1, slice(mods, 0, 3), delta);
    Plaintext pt1 = Plaintext::_from_cmat_without_encoding(mat1, slice(mods, 0, 3), delta);
    // 准备一个私钥
    printf("准备一个私钥\n");
    SecretKey sk(n, p);
    printf("准备自举所需的公钥\n");
    BootstrapKey bsk(n, p, sk, qo, ctx, mods);
    // 加密成密文
    printf("加密成密文\n");
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);
    printf("开始自举\n");
    Ciphertext ct2 = bsk.main(ct1);
    printf("自举结束，检查取值\n");

    Plaintext decrypted = ct2.decrypt(sk, ctx);
    ComplexMatrixGroup mat2 = decrypted._to_cmat_without_decoding(delta);
    // 比较结果
    for(int w=0; w<p-1; w++)
    {
        for(int x=0; x<n; x++)
        {
            for(int y=0; y<n; y++)
            {
                complex got = mat2.at(w, x, y);
                complex expected = mat1.at(w, x, y);
                complex diff = got - expected;
                // printf("%d %d %d %.6lf %.6lf %.6lf\n", w, x, y, got, expected, diff);
                diff *= delta;
                diff /= mods[0];
                diff /= mods[1];
                diff /= mods[2];
                printf("%d %d %d %.6lf %.6lf\n", w, x, y, diff.real(), diff.imag());

            }
        }
    }
}