#ifndef ECC_BIGNUMOPERATION_H
#define ECC_BIGNUMOPERATION_H

#include <stdint.h>

typedef struct BigNum {
    uint32_t* num;
    int size;       // number of 32-bit words in num
    int sign;       // 1 for positive, 0 for negative
}BigNum;

void str_to_hex(const char* str, uint32_t* hex);
void BigNum_copy(BigNum* dest, BigNum* src);
int BigNum_cmp(BigNum* a, BigNum* b);
void BigNum_add(BigNum* a, BigNum* b, BigNum* result);
void mod_add(BigNum* a, BigNum*b, BigNum* p, BigNum* result);
void BigNum_sub_abs(BigNum* larger, BigNum* smaller, BigNum* result);
void BigNum_sub(BigNum* a, BigNum* b, BigNum* result);
void mod_sub(BigNum* a, BigNum* b, BigNum* p, BigNum* result);
void BigNum_OS(BigNum*a, BigNum* b, BigNum* result);
void BigNum_OS_16(BigNum*a, BigNum* b, BigNum* result);
void BigNum_PS(BigNum*a, BigNum* b, BigNum* result);
void BigNum_PS_32(BigNum*a, BigNum* b, BigNum* result);
void BigNum_squaring(BigNum* a, BigNum* result);
void BigNum_Fast_Reduction(BigNum* a, BigNum* p, BigNum* result);
void BigNum_mod_sub_inverse(BigNum* a, BigNum* b, BigNum* p, BigNum* result);
void BigNum_Inverse(BigNum* a, BigNum* p, BigNum* result);
void BigNum_right_shift1(BigNum* x);

#endif
