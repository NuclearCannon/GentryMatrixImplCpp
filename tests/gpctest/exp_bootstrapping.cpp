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
    double delta = 10000;

    GentryPolyCtx ctx(n, p, qrp);

    // 准备一个明文
    ComplexMatrixGroup mat1 = ComplexMatrixGroup::random(100, n, p);
    Plaintext pt1 = Plaintext::from_cmat(mat1, mods, delta);
    // 准备一个私钥
    SecretKey sk(n, p, mods);
    CircledastKey mmkey = CircledastKey::gen(sk, qo, ctx);
    ConjTransposeKey ctkey = ConjTransposeKey::gen(sk, qo, ctx);
    // 加密成密文
    Ciphertext ct1 = Ciphertext::encrypt(pt1, sk, ctx);

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
            RotateKey rkey = RotateKey::gen(sk, qo, ctx, 0,0,l);
            Ciphertext rotated = rkey.run(m02, ctx);
            // 构造临时明文
            ComplexMatrixGroup tmp(n, p);
            complex x = std::polar<double>(1, - (2 * M_PI / ((double)p)) * (double)gamma);
            x = (x-1.0)/((double)p);
            for(int k=0;k<p-1;k++)for(int i=0;i<n;i++)for(int j=0;j<n;j++)tmp.at(k,i,j)=x;
            Plaintext tmppt = Plaintext::from_cmat(tmp, mods, delta);
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