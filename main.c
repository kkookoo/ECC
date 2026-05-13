#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "ECC_BigNumOperation.h"

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

int main()
{
    FILE* fp = fopen("mod_sub_test_vectors.txt", "r");
    if (!fp)
    {
        printf("file open error\n");
        return 1;
    }

    char line[MAX_LINE];
    char input1[300], input2[300], input3[300];

    int test = 0, fail = 0;

    while (fgets(line, sizeof(line), fp))
    {
        if (strncmp(line, "A:", 2) == 0)
        {
            strcpy(input1, line + 2);
            trim(input1);
        }
        else if (strncmp(line, "B:", 2) == 0)
        {
            strcpy(input2, line + 2);
            trim(input2);
        }
        else if (strncmp(line, "R:", 2) == 0)
        {
            strcpy(input3, line + 2);
            trim(input3);

            printf("A:        %s\n", input1);
            printf("B:        %s\n", input2);
            printf("R:        %s\n", input3);
            printf("Expected: %s\n", input3);
            printf("\n");

            BigNum num1, num2, expected, result;

            num1.size = strlen(input1) / 8;
            num1.num = malloc(sizeof(uint32_t) * num1.size);
            str_to_hex(input1, num1.num);

            num2.size = strlen(input2) / 8;
            num2.num = malloc(sizeof(uint32_t) * num2.size);
            str_to_hex(input2, num2.num);

            expected.size = strlen(input3) / 8;
            expected.num = malloc(sizeof(uint32_t) * expected.size);
            str_to_hex(input3, expected.num);

            result.size = (num1.size > num2.size) ? num1.size : num2.size;
            result.num = malloc(sizeof(uint32_t) * (result.size + 1));

            // mod_add(&num1, &num2, &P, &result);
            mod_sub(&num1, &num2, &P, &result);

            test++;

            // 비교
            int ok = 1;
            for (int i = 0; i < expected.size; i++)
            {
                if (result.num[i] != expected.num[i])
                {
                    ok = 0;
                    break;
                }
            }

            if (!ok)
            {
                fail++;
                printf("FAIL #%d\n", test);
                printf("Actual: ");
                for (int i = result.size - 1; i >= 0; i--)
                    printf("%08x", result.num[i]);
                printf("\n\n");
            }

            free(num1.num);
            free(num2.num);
            free(expected.num);
            free(result.num);
        }
    }

    fclose(fp);

    printf("Total: %d, Pass: %d, Fail: %d\n", test, test - fail, fail);

    return 0;
}