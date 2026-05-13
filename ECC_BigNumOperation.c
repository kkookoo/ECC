#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ECC_BigNumOperation.h"

void str_to_hex(const char* str, uint32_t* hex)         // transform string to hex 
{
    uint32_t val[8];
    size_t len = strlen(str);

    for (int i = len; i > 0; i -= 8) 
    {
        for (int j = 0; j < 8; j++)
        {
            val[j] = (str[i - 1 - j] > '9') ? (str[i - 1 - j] - 'A' + 10) : (str[i - 1 - j] - '0');
        }

        hex[(len - i) / 8] = (val[7] << 28) | (val[6] << 24) | 
                             (val[5] << 20) | (val[4] << 16) | 
                             (val[3] << 12) | (val[2] << 8)  | 
                             (val[1] << 4)  | val[0];
    }
}

int BigNum_cmp(BigNum* a, BigNum* b)            // 1 : a > b, 0 : a == b, -1 : a < b
{
    if (a->size > b->size) return 1;
    if (a->size < b->size) return -1;

    for(int i = a->size - 1; i >= 0; i--)
    {
        if (a->num[i] > b->num[i]) return 1;
        if (a->num[i] < b->num[i]) return -1;
    }

    return 0;
}

void BigNum_add(BigNum* a, BigNum* b, BigNum* result)
{
    int carry = 0;
    result->size = (a->size > b->size) ? a->size : b->size; 

    for (int i = 0; i < result->size; i++)
    {
        uint64_t x = (i < a->size) ? a->num[i] : 0x00000000;        // overflow 방지 위해 uint64_t로 선언
        uint64_t y = (i < b->size) ? b->num[i] : 0x00000000;

        uint64_t sum = x + y + carry;
        
        result->num[i] = (uint32_t)sum;

        carry = (sum > 0xffffffff) ? 1 : 0;
    }

    if (carry) 
    {
        result->num[result->size] = (uint32_t)carry;
        result->size++;
    }
}

void mod_add(BigNum* a, BigNum*b, BigNum* p, BigNum* result)
{
    BigNum_add(a, b, result);

    if(result->size > p->size || BigNum_cmp(result, p) >= 0)
    {
        BigNum_sub(result, p, result);
    }
}

void BigNum_sub(BigNum* a, BigNum* b, BigNum* result)
{
    int borrow = 0;
    result->size = (a->size > b->size) ? a->size : b->size;

    if(BigNum_cmp(a, b) < 0) result-> sign = 0;
    if(BigNum_cmp(a, b) >= 0) result-> sign = 1;

    for (int i = 0; i < result->size; i++)
    {
        uint32_t x = (i < a->size) ? a->num[i] : 0x00000000;
        uint32_t y = (i < b->size) ? b->num[i] : 0x00000000;    

        uint32_t diff = x - y - borrow;

        result->num[i] = diff;

        borrow = (x < (uint64_t)y + borrow) ? 1 : 0;
    }
}

void mod_sub(BigNum* a, BigNum* b, BigNum* p, BigNum* result)
{
    BigNum_sub(a, b, result);

    if(result->sign == 0)
    {
        BigNum_add(result, p, result);
    }
}