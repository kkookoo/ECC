#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "ECC_BigNumOperation.h"
#include "Elliptic_arithmetic.h"

#define MAX_LINE 256

BigNum P = 
{
    .num = (uint32_t[]){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 
                        0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF},
    .size = 8,
    .sign = 1
};

void trim(char* s)
{
    while (*s == ' ') memmove(s, s+1, strlen(s));
    s[strcspn(s, "\n")] = 0;
}

void print_words(BigNum* x)
{
    for (int i = 7; i >= 0; i--)
    {
        if (i < x->size)
            printf("%08x", x->num[i]);
        else
            printf("%08x", 0);
    }
}

int BigNum_equal_8(BigNum* a, BigNum* b)
{
    for (int i = 0; i < 8; i++)
    {
        uint32_t a_word = (i < a->size) ? a->num[i] : 0;
        uint32_t b_word = (i < b->size) ? b->num[i] : 0;

        if (a_word != b_word)
            return 0;
    }

    return 1;
}

void set_bignum_from_hex(BigNum* x, char* hex, uint32_t* buf)
{
    x->size = 8;
    x->sign = 1;
    x->num = buf;
    memset(buf, 0, sizeof(uint32_t) * 8);
    str_to_hex(hex, buf);
}

int main()
{
    FILE* fp = fopen("affine_affine_add_vectors.txt", "r");
    if (!fp)
    {
        printf("file open error\n");
        return 1;
    }

    char line[MAX_LINE];
    char p_x_hex[80], p_y_hex[80], q_x_hex[80], q_y_hex[80];
    char r_x_hex[80], r_y_hex[80];
    int id = -1;
    int p_infinity = 0;
    int q_infinity = 0;
    int r_infinity = 0;

    int test = 0, fail = 0;

    while (fgets(line, sizeof(line), fp))
    {
        if (strncmp(line, "ID:", 3) == 0)
        {
            id = atoi(line + 3);
        }
        else if (strncmp(line, "P_INFINITY:", 11) == 0)
        {
            p_infinity = atoi(line + 11);
        }
        else if (strncmp(line, "P_X:", 4) == 0)
        {
            strcpy(p_x_hex, line + 4);
            trim(p_x_hex);
        }
        else if (strncmp(line, "P_Y:", 4) == 0)
        {
            strcpy(p_y_hex, line + 4);
            trim(p_y_hex);
        }
        else if (strncmp(line, "Q_INFINITY:", 11) == 0)
        {
            q_infinity = atoi(line + 11);
        }
        else if (strncmp(line, "Q_X:", 4) == 0)
        {
            strcpy(q_x_hex, line + 4);
            trim(q_x_hex);
        }
        else if (strncmp(line, "Q_Y:", 4) == 0)
        {
            strcpy(q_y_hex, line + 4);
            trim(q_y_hex);
        }
        else if (strncmp(line, "R_INFINITY:", 11) == 0)
        {
            r_infinity = atoi(line + 11);
        }
        else if (strncmp(line, "R_X:", 4) == 0)
        {
            strcpy(r_x_hex, line + 4);
            trim(r_x_hex);
        }
        else if (strncmp(line, "R_Y:", 4) == 0)
        {
            strcpy(r_y_hex, line + 4);
            trim(r_y_hex);

            uint32_t p_x_buf[8], p_y_buf[8], q_x_buf[8], q_y_buf[8];
            uint32_t exp_x_buf[8], exp_y_buf[8];
            BigNum exp_x, exp_y;
            Point p_input, q_input, result;

            set_bignum_from_hex(&p_input.x, p_x_hex, p_x_buf);
            set_bignum_from_hex(&p_input.y, p_y_hex, p_y_buf);
            p_input.infinity = p_infinity;
            set_bignum_from_hex(&q_input.x, q_x_hex, q_x_buf);
            set_bignum_from_hex(&q_input.y, q_y_hex, q_y_buf);
            q_input.infinity = q_infinity;

            set_bignum_from_hex(&exp_x, r_x_hex, exp_x_buf);
            set_bignum_from_hex(&exp_y, r_y_hex, exp_y_buf);

            result.x.size = 8;
            result.y.size = 8;
            result.x.sign = 1;
            result.y.sign = 1;
            result.x.num = calloc(9, sizeof(uint32_t));
            result.y.num = calloc(9, sizeof(uint32_t));

            point_add(&p_input, &q_input, &result, &P);

            test++;

            int ok = ((int)result.infinity == r_infinity) &&
                     BigNum_equal_8(&result.x, &exp_x) &&
                     BigNum_equal_8(&result.y, &exp_y);

            printf("ID: %d\n", id);
            printf("P_INFINITY: %d\n", p_infinity);
            printf("P_X:        ");
            print_words(&p_input.x);
            printf("\nP_Y:        ");
            print_words(&p_input.y);
            printf("\nQ_INFINITY: %d", q_infinity);
            printf("\nQ_X:        ");
            print_words(&q_input.x);
            printf("\nQ_Y:        ");
            print_words(&q_input.y);
            printf("\nActual Infinity:   %u", result.infinity);
            printf("\nExpected Infinity: %d", r_infinity);
            printf("\nActual X:   ");
            print_words(&result.x);
            printf("\nExpected X: ");
            print_words(&exp_x);
            printf("\nActual Y:   ");
            print_words(&result.y);
            printf("\nExpected Y: ");
            print_words(&exp_y);
            printf("\nRESULT:     %s\n\n", ok ? "PASS" : "FAIL");

            if (!ok)
            {
                fail++;
            }

            free(result.x.num);
            free(result.y.num);
        }
    }

    fclose(fp);

    printf("Total: %d, Pass: %d, Fail: %d\n", test, test - fail, fail);
    printf("\n");

    return 0;
}
