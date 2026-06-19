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

void BigNum_copy(BigNum* dest, BigNum* src)
{
    dest->size = 8;
    dest->sign = src->sign;

    for (int i = 0; i < 8; i++)
        dest->num[i] = src->num[i];
}

void BigNum_size(BigNum* a)
{
    int n = a->size;

    while (n > 1 && a->num[n - 1] == 0)
        n--;

    a->size = n;
}

int BigNum_cmp(BigNum* a, BigNum* b)    // 1 : a > b, 0 : a == b, -1 : a < b
{                                       // 부채널 공격 방지를 위해 큰 수 비교 함수 구현 (비트 마스킹 연산)
    int result = 0;

    BigNum_size(a);
    BigNum_size(b);

    if (a->size > b->size) return 1;
    if (a->size < b->size) return -1;

    for (int i = a->size - 1; i >= 0; i--)
    {
        uint32_t gt = (a->num[i] > b->num[i]);
        uint32_t lt = (a->num[i] < b->num[i]);
        uint32_t keep = (uint32_t)(-(int32_t)(result == 0));
        result |= (int)(keep & (gt - lt));
    }

    return (result > 0) - (result < 0);
}


void BigNum_add(BigNum* a, BigNum* b, BigNum* result)
{
    int carry = 0;
    result->size = (a->size > b->size) ? a->size : b->size; 

    for (int i = 0; i < result->size; i++)
    {
        uint64_t x = (i < a->size) ? a->num[i] : 0x00000000;        // overflow 방지를 위해 uint64_t로 선언
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
    result->sign = 1;

    while(result->size > p->size || BigNum_cmp(result, p) >= 0)
    {
        BigNum_sub(result, p, result);
        result->sign = 1;
    }
}

void BigNum_sub_abs(BigNum* larger, BigNum* smaller, BigNum* result)
  {
      int borrow = 0;
      result->size = larger->size;

      for (int i = 0; i < result->size; i++)
      {
          uint32_t x = (i < larger->size) ? larger->num[i] : 0;
          uint32_t y = (i < smaller->size) ? smaller->num[i] : 0;

          uint32_t diff = x - y - borrow;
          result->num[i] = diff;

          borrow = (x < (uint64_t)y + borrow) ? 1 : 0;
      }

      while (result->size > 1 && result->num[result->size - 1] == 0)
          result->size--;
  }

  void BigNum_sub(BigNum* a, BigNum* b, BigNum* result)
  {
      int cmp = BigNum_cmp(a, b);

      if (cmp >= 0)
      {
          result->sign = 1;
          BigNum_sub_abs(a, b, result);
      }
      else
      {
          result->sign = 0;
          BigNum_sub_abs(b, a, result);
      }
  }

void mod_sub(BigNum* a, BigNum* b, BigNum* p, BigNum* result)
{
    uint32_t diff_num[9] = {0};
    BigNum diff;
    int cmp;

    diff.num = diff_num;
    diff.size = 1;
    diff.sign = 1;

    cmp = BigNum_cmp(a, b);

    if (cmp >= 0)
        BigNum_sub_abs(a, b, &diff);
    else
        BigNum_sub_abs(b, a, &diff);

    diff.sign = 1;

    while(diff.size > p->size || BigNum_cmp(&diff, p) >= 0)
    {
        BigNum_sub_abs(&diff, p, &diff);
        diff.sign = 1;
    }

    if (cmp >= 0 || (diff.size == 1 && diff.num[0] == 0))
    {
        for (int i = 0; i < p->size; i++)
            result->num[i] = (i < diff.size) ? diff.num[i] : 0;
    }
    else
    {
        BigNum_sub_abs(p, &diff, result);
    }

    result->sign = 1;
    result->size = p->size;
}

void BigNum_OS(BigNum* a, BigNum* b, BigNum* result)
{
    int n = a->size + b->size;
    result->sign = 1;

    for (int i = 0; i < n; i++)
        result->num[i] = 0;

    for (int i = 0; i < a->size; i++)
    {
        uint64_t carry = 0;

        for (int j = 0; j < b->size; j++)
        {
            uint64_t prod = (uint64_t)a->num[i] * b->num[j];
            uint64_t cur = (uint64_t)result->num[i + j] + (uint32_t)prod + carry;

            result->num[i + j] = (uint32_t)cur;
            carry = (prod >> 32) + (cur >> 32);
        }

        int k = i + b->size;

        while (carry)
        {
            uint64_t tail = (uint64_t)result->num[k] + (uint32_t)carry;
            result->num[k] = (uint32_t)tail;
            carry = (carry >> 32) + (tail >> 32);
            k++;
        }
    }

    result->size = n;
}

void BigNum_OS_16(BigNum* a, BigNum* b, BigNum* result)
{
    int n = a->size + b->size;
    result->sign = 1;

    for (int i = 0; i < n; i++)
        result->num[i] = 0;

    for (int i = 0; i < a->size; i++)
    {
        uint64_t carry = 0;

        for (int j = 0; j < b->size; j++)
        {
            uint32_t a1 = a->num[i] & 0x0000ffff;
            uint32_t a2 = (a->num[i] >> 16) & 0x0000ffff;
            uint32_t b1 = b->num[j] & 0x0000ffff;
            uint32_t b2 = (b->num[j] >> 16) & 0x0000ffff;

            uint32_t x0 = a1 * b1;
            uint32_t x1 = a1 * b2;
            uint32_t x2 = a2 * b1;
            uint32_t x3 = a2 * b2;

            uint32_t mid  = (x0 >> 16) + (x1 & 0x0000ffff) + (x2 & 0x0000ffff);
            uint32_t low  = (x0 & 0x0000ffff) | ((mid & 0x0000ffff) << 16);
            uint32_t high = x3 + (x1 >> 16) + (x2 >> 16) + (mid >> 16);

            uint64_t cur = (uint64_t)result->num[i + j] + low + (uint32_t)carry;
            result->num[i + j] = (uint32_t)cur;
            carry = (uint64_t)high + (carry >> 32) + (cur >> 32);
        }

        int k = i + b->size;

        while (carry)
        {
            uint64_t tail = (uint64_t)result->num[k] + (uint32_t)carry;
            result->num[k] = (uint32_t)tail;
            carry = (carry >> 32) + (tail >> 32);
            k++;
        }
    }

    result->size = n;
}

void BigNum_PS(BigNum* a, BigNum* b, BigNum* result)
{
    result->sign = 1;

    int n = a->size + b->size;
    for (int i = 0; i < result->size; i++)
        result->num[i] = 0;

    uint32_t c0 = 0, c1 = 0, c2 = 0;

    for (int k = 0; k < n - 1; k++)
    {
        int i_start = (k < b->size) ? 0 : (k - b->size + 1);
        int i_end   = (k < a->size - 1) ? k : (a->size - 1);

        for (int i = i_start; i <= i_end; i++)
        {
            int j = k - i;                  // i + j = k

            uint32_t a1 = a->num[i] & 0x0000ffff;
            uint32_t a2 = (a->num[i] >> 16) & 0x0000ffff;
            uint32_t b1 = b->num[j] & 0x0000ffff;
            uint32_t b2 = (b->num[j] >> 16) & 0x0000ffff;

            uint32_t x0 = a1 * b1;
            uint32_t x1 = a1 * b2;
            uint32_t x2 = a2 * b1;
            uint32_t x3 = a2 * b2;

            uint32_t mid  = (x0 >> 16) + (x1 & 0x0000ffff) + (x2 & 0x0000ffff);
            uint32_t low  = (x0 & 0x0000ffff) | ((mid & 0x0000ffff) << 16);
            uint32_t high = x3 + (x1 >> 16) + (x2 >> 16) + (mid >> 16);

            uint64_t s = (uint64_t)c0 + low;
            c0 = (uint32_t)s;

            uint64_t t = (uint64_t)c1 + high + (uint32_t)(s >> 32);
            c1 = (uint32_t)t;

            c2 += (uint32_t)(t >> 32);
        }

        result->num[k] = c0;
        c0 = c1;
        c1 = c2;
        c2 = 0;
    }

    result->num[n - 1] = c0;

    result->size = n;
}

void BigNum_PS_32(BigNum* a, BigNum* b, BigNum* result)
{
    result->sign = 1;

    int n = a->size + b->size;
    for (int i = 0; i < result->size; i++)
        result->num[i] = 0;

    uint32_t c0 = 0, c1 = 0, c2 = 0;

    for (int k = 0; k < n - 1; k++)
    {
        int i_start = (k < b->size) ? 0 : (k - b->size + 1);
        int i_end   = (k < a->size - 1) ? k : (a->size - 1);

        for (int i = i_start; i <= i_end; i++)
        {
            int j = k - i;
            uint64_t prod = (uint64_t)a->num[i] * b->num[j];

            uint64_t s = (uint64_t)c0 + (uint32_t)prod;
            c0 = (uint32_t)s;

            uint64_t t = (uint64_t)c1 + (prod >> 32) + (s >> 32);
            c1 = (uint32_t)t;
            c2 += (uint32_t)(t >> 32);
        }

        result->num[k] = c0;
        c0 = c1;
        c1 = c2;
        c2 = 0;
    }

    result->num[n - 1] = c0;
    result->size = n;
}

void BigNum_PS_32_unrolled(BigNum* a, BigNum* b, BigNum* result)
{
    if (a->size != 8 || b->size != 8)
    {
        BigNum_PS_32(a, b, result);
        return;
    }

    result->sign = 1;

    for (int i = 0; i < result->size; i++)
        result->num[i] = 0;

    uint32_t c0 = 0, c1 = 0, c2 = 0;

#define ADDMUL(i, j) do {                                                     \
        uint64_t prod = (uint64_t)a->num[(i)] * b->num[(j)];                  \
        uint64_t s = (uint64_t)c0 + (uint32_t)prod;                           \
        c0 = (uint32_t)s;                                                     \
        uint64_t t = (uint64_t)c1 + (prod >> 32) + (s >> 32);                 \
        c1 = (uint32_t)t;                                                     \
        c2 += (uint32_t)(t >> 32);                                            \
    } while (0)

#define STORE(k) do {                                                         \
        result->num[(k)] = c0;                                                \
        c0 = c1;                                                              \
        c1 = c2;                                                              \
        c2 = 0;                                                               \
    } while (0)

    ADDMUL(0, 0); STORE(0);
    ADDMUL(0, 1); ADDMUL(1, 0); STORE(1);
    ADDMUL(0, 2); ADDMUL(1, 1); ADDMUL(2, 0); STORE(2);
    ADDMUL(0, 3); ADDMUL(1, 2); ADDMUL(2, 1); ADDMUL(3, 0); STORE(3);
    ADDMUL(0, 4); ADDMUL(1, 3); ADDMUL(2, 2); ADDMUL(3, 1); ADDMUL(4, 0); STORE(4);
    ADDMUL(0, 5); ADDMUL(1, 4); ADDMUL(2, 3); ADDMUL(3, 2); ADDMUL(4, 1); ADDMUL(5, 0); STORE(5);
    ADDMUL(0, 6); ADDMUL(1, 5); ADDMUL(2, 4); ADDMUL(3, 3); ADDMUL(4, 2); ADDMUL(5, 1); ADDMUL(6, 0); STORE(6);
    ADDMUL(0, 7); ADDMUL(1, 6); ADDMUL(2, 5); ADDMUL(3, 4); ADDMUL(4, 3); ADDMUL(5, 2); ADDMUL(6, 1); ADDMUL(7, 0); STORE(7);
    ADDMUL(1, 7); ADDMUL(2, 6); ADDMUL(3, 5); ADDMUL(4, 4); ADDMUL(5, 3); ADDMUL(6, 2); ADDMUL(7, 1); STORE(8);
    ADDMUL(2, 7); ADDMUL(3, 6); ADDMUL(4, 5); ADDMUL(5, 4); ADDMUL(6, 3); ADDMUL(7, 2); STORE(9);
    ADDMUL(3, 7); ADDMUL(4, 6); ADDMUL(5, 5); ADDMUL(6, 4); ADDMUL(7, 3); STORE(10);
    ADDMUL(4, 7); ADDMUL(5, 6); ADDMUL(6, 5); ADDMUL(7, 4); STORE(11);
    ADDMUL(5, 7); ADDMUL(6, 6); ADDMUL(7, 5); STORE(12);
    ADDMUL(6, 7); ADDMUL(7, 6); STORE(13);
    ADDMUL(7, 7); STORE(14);

    result->num[15] = c0;
    result->size = 16;

#undef STORE
#undef ADDMUL
}

void BigNum_squaring(BigNum* a, BigNum* result)
{
    result->sign = 1;

    int n = 2 * a->size;
    for (int i = 0; i < result->size; i++)
        result->num[i] = 0;

    uint32_t c0 = 0, c1 = 0, c2 = 0;

    for (int k = 0; k < n - 1; k++)
    {
        int i_start = (k < a->size) ? 0 : (k - a->size + 1);
        int i_end   = (k < a->size - 1) ? k : (a->size - 1);

        for (int i = i_start; i <= i_end; i++)
        {
            int j = k - i;                  // i + j = k
            if(i > j) break;

            uint32_t a1 = a->num[i] & 0x0000ffff;
            uint32_t a2 = (a->num[i] >> 16) & 0x0000ffff;
            uint32_t b1 = a->num[j] & 0x0000ffff;
            uint32_t b2 = (a->num[j] >> 16) & 0x0000ffff;

            uint32_t x0 = a1 * b1;
            uint32_t x1 = a1 * b2;
            uint32_t x2 = a2 * b1;
            uint32_t x3 = a2 * b2;

            uint32_t mid  = (x0 >> 16) + (x1 & 0x0000ffff) + (x2 & 0x0000ffff);
            uint32_t low  = (x0 & 0x0000ffff) | ((mid & 0x0000ffff) << 16);
            uint32_t high = x3 + (x1 >> 16) + (x2 >> 16) + (mid >> 16);

            if(i < j)
            {
                uint64_t s = (uint64_t)c0 + 2 * (uint64_t)(low);
                c0 = (uint32_t)s;
                uint64_t t = (uint64_t)c1 + 2 * (uint64_t)(high) + (uint32_t)(s >> 32);
                c1 = (uint32_t)t;
                c2 += (uint32_t)(t >> 32);
            }
            else if(i == j)
            {
                uint64_t s = (uint64_t)c0 + low;
                c0 = (uint32_t)s;
                uint64_t t = (uint64_t)c1 + high + (uint32_t)(s >> 32);
                c1 = (uint32_t)t;
                c2 += (uint32_t)(t >> 32);
            }
        }
        
        result->num[k] = c0;
        c0 = c1;
        c1 = c2;
        c2 = 0;
    }

    result->num[n - 1] = c0;

    result->size = n;
}

void BigNum_Fast_Reduction(BigNum* a, BigNum* p, BigNum* result)
{
    BigNum s1, s2, s3, s4, s5, s6, s7, s8, s9;

    s1.size = 8;
    s2.size = 8;
    s3.size = 8;
    s4.size = 8;
    s5.size = 8;
    s6.size = 8;
    s7.size = 8;
    s8.size = 8;
    s9.size = 8;

    s1.num = malloc(sizeof(uint32_t) * 9);
    s2.num = malloc(sizeof(uint32_t) * 9);
    s3.num = malloc(sizeof(uint32_t) * 9);
    s4.num = malloc(sizeof(uint32_t) * 9);
    s5.num = malloc(sizeof(uint32_t) * 9);
    s6.num = malloc(sizeof(uint32_t) * 9);
    s7.num = malloc(sizeof(uint32_t) * 9);
    s8.num = malloc(sizeof(uint32_t) * 9);
    s9.num = malloc(sizeof(uint32_t) * 9);

    result->sign = 1;
    result->size = 8;

    for(int i = 0; i < 8; i++) s1.num[i] = a->num[i];

    s2.num[0] = 0x00000000; s2.num[1] = 0x00000000; s2.num[2] = 0x00000000; s2.num[3] = a->num[11];
    s2.num[4] = a->num[12]; s2.num[5] = a->num[13]; s2.num[6] = a->num[14]; s2.num[7] = a->num[15];

    s3.num[0] = 0x00000000; s3.num[1] = 0x00000000; s3.num[2] = 0x00000000; s3.num[3] = a->num[12];
    s3.num[4] = a->num[13]; s3.num[5] = a->num[14]; s3.num[6] = a->num[15]; s3.num[7] = 0x00000000;

    s4.num[0] = a->num[8]; s4.num[1] = a->num[9]; s4.num[2] = a->num[10]; s4.num[3] = 0x00000000;
    s4.num[4] = 0x00000000; s4.num[5] = 0x00000000; s4.num[6] = a->num[14]; s4.num[7] = a->num[15];

    s5.num[0] = a->num[9]; s5.num[1] = a->num[10]; s5.num[2] = a->num[11]; s5.num[3] = a->num[13];
    s5.num[4] = a->num[14]; s5.num[5] = a->num[15]; s5.num[6] = a->num[13]; s5.num[7] = a->num[8];

    s6.num[0] = a->num[11]; s6.num[1] = a->num[12]; s6.num[2] = a->num[13]; s6.num[3] = 0x00000000;
    s6.num[4] = 0x00000000; s6.num[5] = 0x00000000; s6.num[6] = a->num[8]; s6.num[7] = a->num[10];

    s7.num[0] = a->num[12]; s7.num[1] = a->num[13]; s7.num[2] = a->num[14]; s7.num[3] = a->num[15];
    s7.num[4] = 0x00000000; s7.num[5] = 0x00000000; s7.num[6] = a->num[9]; s7.num[7] = a->num[11];

    s8.num[0] = a->num[13]; s8.num[1] = a->num[14]; s8.num[2] = a->num[15]; s8.num[3] = a->num[8];
    s8.num[4] = a->num[9]; s8.num[5] = a->num[10]; s8.num[6] = 0x00000000; s8.num[7] = a->num[12];

    s9.num[0] = a->num[14]; s9.num[1] = a->num[15]; s9.num[2] = 0x00000000; s9.num[3] = a->num[9];
    s9.num[4] = a->num[10]; s9.num[5] = a->num[11]; s9.num[6] = 0x00000000; s9.num[7] = a->num[13];

    mod_add(&s2, &s2, p, &s2);
    mod_add(&s3, &s3, p, &s3);
    mod_add(&s1, &s2, p, &s1);
    mod_add(&s1, &s3, p, &s1);
    mod_add(&s1, &s4, p, &s1);
    mod_add(&s1, &s5, p, &s1);
    mod_sub(&s1, &s6, p, &s1);
    mod_sub(&s1, &s7, p, &s1);
    mod_sub(&s1, &s8, p, &s1);
    mod_sub(&s1, &s9, p, result);

    free(s1.num);
    free(s2.num);
    free(s3.num);
    free(s4.num);
    free(s5.num);
    free(s6.num);
    free(s7.num);
    free(s8.num);
    free(s9.num);
}

void BigNum_right_shift1(BigNum* x)
{
    uint32_t carry = 0;

    for (int i = x->size - 1; i >= 0; i--)
    {
        uint32_t new_carry = x->num[i] & 1;
        x->num[i] = (x->num[i] >> 1) | (carry << 31);
        carry = new_carry;
    }

    while (x->size > 1 && x->num[x->size - 1] == 0)
        x->size--;
}

void BigNum_mod_sub_inverse(BigNum* a, BigNum* b, BigNum* p, BigNum* result)        
{
    BigNum tmp;

    tmp.size = 8;
    tmp.sign = 1;
    tmp.num = calloc(9, sizeof(uint32_t));

    if(BigNum_cmp(a, b) >= 0)
    {
        BigNum_sub(a, b, result);
    }
    else
    {
        BigNum_sub(b, a, &tmp);
        BigNum_sub(p, &tmp, result);
    }

    result->sign = 1;
    result->size = 8;

    free(tmp.num);
}

void BigNum_Inverse(BigNum* a, BigNum* p, BigNum* result)
{
    BigNum u, v, x1, x2;

    u.size = 8;
    v.size = 8;
    x1.size = 8;
    x2.size = 8;

    u.num = calloc(9, sizeof(uint32_t));
    v.num = calloc(9, sizeof(uint32_t));
    x1.num = calloc(9, sizeof(uint32_t));
    x2.num = calloc(9, sizeof(uint32_t));

    u.sign = 1;
    v.sign = 1;
    x1.sign = 1;
    x2.sign = 1;

    result->sign = 1;

    for(int i = 0; i < 8; i++)
    {
        u.num[i] = a->num[i];
        v.num[i] = p->num[i];
        if(i == 0) x1.num[i] = 1;
        else x1.num[i] = 0;
        x2.num[i] = 0;
    }

    while((u.size != 1 || u.num[0] != 1) && (v.size != 1 || v.num[0] != 1))
    {
        while((u.num[0] & 1) == 0)
        {
            BigNum_right_shift1(&u);

            if((x1.num[0] & 1) == 0)
            {
                BigNum_right_shift1(&x1);
            }
            else
            {
                BigNum_add(&x1, p, &x1);
                BigNum_right_shift1(&x1);
            }
        }

        while((v.num[0] & 1) == 0)
        {
            BigNum_right_shift1(&v);

            if((x2.num[0] & 1) == 0)
            {
                BigNum_right_shift1(&x2);
            }
            else
            {
                BigNum_add(&x2, p, &x2);
                BigNum_right_shift1(&x2);
            }
        }

        if(BigNum_cmp(&u, &v) >= 0)
        {
            BigNum_sub(&u, &v, &u);
            BigNum_mod_sub_inverse(&x1, &x2, p, &x1);
        }
        
        else
        {
            BigNum_sub(&v, &u, &v);
            BigNum_mod_sub_inverse(&x2, &x1, p, &x2);
        }
    }

    if(u.size == 1 && u.num[0] == 1)
    {
        if(x1.sign == 0)
        {
            BigNum_add(&x1, p, &x1);
        }

        for(int i = 0; i < 8; i++)
            result->num[i] = x1.num[i];
    }
    else
    {
        if(x2.sign == 0)
        {
            BigNum_add(&x2, p, &x2);
        }
        
        for(int i = 0; i < 8; i++)
            result->num[i] = x2.num[i];
    }

    free(u.num);
    free(v.num);
    free(x1.num);
    free(x2.num);
}
