#include "application.hpp"


void test_ntt();
void test_ntt_cuda();
void test_encrypt();
void test_encrypt_cuda();
void test_ks();
void test_ks_cuda();
void test_circledast();
void test_circledast_cuda();
int main()
{
    // test_ntt();
    // test_ntt_cuda();
    // test_encrypt();
    // test_encrypt_cuda();
    // test_ks();
    // test_ks_cuda();
    // test_circledast();
    // test_circledast_cuda();
    // Plaintext::test_pt_encode_and_decode();
    // Ciphertext::test_ct_encrypt_and_decrypt();
    // CircledastKey::test_pt_circledast_end2end();
    // CircledastKey::test_ct_circledast_end2end();
    // CircledastKey::test_ct_circledast_end2end(true);
    // ConjTransposeKey::test_pt_transpose();
    // ConjTransposeKey::test_ct_transpose();
    // MultKey::test_ct_mult();
    RotateKey::test_pt_rotate();
    return 0;
}