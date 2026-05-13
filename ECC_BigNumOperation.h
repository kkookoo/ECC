#ifndef ECC_BIGNUMOPERATION_H
#define ECC_BIGNUMOPERATION_H

#include <stdint.h>

typedef struct BigNum {
    uint32_t* num;
    int size;       // number of 32-bit words in num
    int sign;       // 1 for positive, 0 for negative
}BigNum;

void str_to_hex(const char* str, uint32_t* hex);
int BigNum_cmp(BigNum* a, BigNum* b);
void BigNum_add(BigNum* a, BigNum* b, BigNum* result);
void mod_add(BigNum* a, BigNum*b, BigNum* p, BigNum* result);
void BigNum_sub(BigNum* a, BigNum* b, BigNum* result);
void mod_sub(BigNum* a, BigNum* b, BigNum* p, BigNum* result);

#endif