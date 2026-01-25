#include "FHE/key_switch.hpp"
#include "FHE/circledast_ct.hpp"
#include <memory>
#include <chrono>
#include <iostream>


std::pair<std::unique_ptr<KeySwitchingKey>, std::unique_ptr<KeySwitchingKey>> create_ksks_for_circledast(
    const ZiqArray& sk,
    const ZiqArrayContext* ctx_q,
    const ZiqArrayContext* ctx_Q,
    const fmpz_t B,
    int L
)
{
    ZiqArray s1 = sk.transpose().conj().w_inversion();
    ZiqArray s2 = sk.mul_poly(s1);
    return std::make_pair(
        std::make_unique<KeySwitchingKey>(s1, sk, ctx_q, ctx_Q, B, L),
        std::make_unique<KeySwitchingKey>(s2, sk, ctx_q, ctx_Q, B, L)
    );
}

/*
def circledast_ct(ct_u, ct_v, ksks):
    au, bu = ct_u
    av, bv = ct_v
    assert isinstance(au, ZqiXYW)
    assert isinstance(av, ZqiXYW)
    assert isinstance(bu, ZqiXYW)
    assert isinstance(bv, ZqiXYW)
    
    ksk1: KeySwitchingKey = ksks[0]
    ksk2: KeySwitchingKey = ksks[1]

    auav = au.circledast(av) # 此项需要以s(X,W)\overline{s}(Y, W^{-1})为KSK源
    aubv = au.circledast(bv) # 此项权重为s
    buav = bu.circledast(av) # 此项需要以\overline{s}(Y, W^{-1})为KSK源
    bubv = bu.circledast(bv) # 此项权重为1

    buav_ksed = ksk1.key_switch_big(buav)
    auav_ksed = ksk2.key_switch_big(auav)

    a = aubv + buav_ksed[0] + auav_ksed[0]
    b = bubv + buav_ksed[1] + auav_ksed[1]
    return a, b
*/

std::pair<ZiqArray, ZiqArray> circledast_ct(
    const ZiqArray& au,
    const ZiqArray& bu,
    const ZiqArray& av,
    const ZiqArray& bv,
    const KeySwitchingKey& ksk1,
    const KeySwitchingKey& ksk2
)
{
    // 需要先转化为half形式才能circledast
    ZiqArray auh = au.iw_ntt();
    ZiqArray avh = av.iw_ntt();
    ZiqArray buh = bu.iw_ntt();
    ZiqArray bvh = bv.iw_ntt();

    ZiqArray auav = auh.circledast(avh).iw_intt();
    ZiqArray aubv = auh.circledast(bvh).iw_intt();
    ZiqArray buav = buh.circledast(avh).iw_intt();
    ZiqArray bubv = buh.circledast(bvh).iw_intt();
    auto [buav0, buav1] = ksk1.key_switch_big_1(buav);
    auto [auav0, auav1] = ksk2.key_switch_big_1(auav);
    return std::make_pair(
        aubv.add(buav0).add(auav0),
        bubv.add(buav1).add(auav1)
    );
}


std::pair<ZiqArray, ZiqArray> circledast_ct_timer(
    const ZiqArray& au,
    const ZiqArray& bu,
    const ZiqArray& av,
    const ZiqArray& bv,
    const KeySwitchingKey& ksk1,
    const KeySwitchingKey& ksk2
)
{
    using namespace std::chrono;
    auto t0 = high_resolution_clock::now();
    // 第一段：iw_ntt
    ZiqArray auh = au.iw_ntt();
    ZiqArray avh = av.iw_ntt();
    ZiqArray buh = bu.iw_ntt();
    ZiqArray bvh = bv.iw_ntt();
    auto t1 = high_resolution_clock::now();
    std::cout << "[circledast_ct_timer] iw_ntt      耗时: " << duration_cast<microseconds>(t1 - t0).count() << " us" << std::endl;

    // 第二段：circledast
    ZiqArray auavh = auh.circledast(avh);
    ZiqArray aubvh = auh.circledast(bvh);
    ZiqArray buavh = buh.circledast(avh);
    ZiqArray bubvh = buh.circledast(bvh);
    auto t2 = high_resolution_clock::now();
    std::cout << "[circledast_ct_timer] circledast  耗时: " << duration_cast<microseconds>(t2 - t1).count() << " us" << std::endl;

    // 第三段：iw_intt
    ZiqArray auav = auavh.iw_intt();
    ZiqArray aubv = aubvh.iw_intt();
    ZiqArray buav = buavh.iw_intt();
    ZiqArray bubv = bubvh.iw_intt();
    auto t3 = high_resolution_clock::now();
    std::cout << "[circledast_ct_timer] iw_intt     耗时: " << duration_cast<microseconds>(t3 - t2).count() << " us" << std::endl;

    // 第四段：Key Switch
    auto [buav0, buav1] = ksk1.key_switch_big_1(buav);
    auto [auav0, auav1] = ksk2.key_switch_big_1(auav);
    auto t4 = high_resolution_clock::now();
    std::cout << "[circledast_ct_timer] KeySwitch   耗时: " << duration_cast<microseconds>(t4 - t3).count() << " us" << std::endl;

    return std::make_pair(
        aubv.add(buav0).add(auav0),
        bubv.add(buav1).add(auav1)
    );
}