#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
#include "ECC_BigNumOperation.h"
#include "Elliptic_arithmetic.h"

#define N 10000
#define WARMUP 1000

BigNum P = {
    .num = (uint32_t[]){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000,
                        0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF},
    .size = 8,
    .sign = 1
};

static BigNum a[N], b[N], red_in[N];
static uint32_t a_buf[N][8];
static uint32_t b_buf[N][8];
static uint32_t red_in_buf[N][16];
static Point_Jacobian jacobian_points[N];
static uint32_t jacobian_x_buf[N][8];
static uint32_t jacobian_y_buf[N][8];
static uint32_t jacobian_z_buf[N][8];
static Point affine_points_p[N], affine_points_q[N];

typedef struct BenchResult {
    uint64_t total;
    double avg;
} BenchResult;

static inline uint64_t read_cycles(void)
{
    unsigned int aux;
    _mm_lfence();
    uint64_t t = __rdtscp(&aux);
    _mm_lfence();
    return t;
}

static uint32_t xorshift32(void)
{
    static uint32_t x = 0x12345678;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

static void random_bn_256(BigNum* x)
{
    x->size = 8;
    x->sign = 1;

    for (int i = 0; i < 8; i++)
        x->num[i] = xorshift32();

    x->num[7] &= 0x7fffffff;

    if (x->num[0] == 0)
        x->num[0] = 1;
}

static double bench_mod_add(void)
{
    BigNum r;
    uint32_t r_buf[9];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 8;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        mod_add(&a[i], &b[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        mod_add(&a[i], &b[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_mod_sub(void)
{
    BigNum r;
    uint32_t r_buf[9];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 8;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        mod_sub(&a[i], &b[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        mod_sub(&a[i], &b[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_os(void)
{
    BigNum r;
    uint32_t r_buf[16];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 16;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_OS(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_OS(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static BenchResult make_bench_result(uint64_t total)
{
    BenchResult result;

    result.total = total;
    result.avg = (double)total / N;

    return result;
}

static BenchResult bench_os_16(void)
{
    BigNum r;
    uint32_t r_buf[16];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 16;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_OS_16(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_OS_16(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return make_bench_result(end - start);
}

static BenchResult bench_os_16_reduction(void)
{
    BigNum product, reduced;
    uint32_t product_buf[16];
    uint32_t reduced_buf[9];
    volatile uint32_t sink = 0;

    product.num = product_buf;
    product.size = 16;
    product.sign = 1;
    reduced.num = reduced_buf;
    reduced.size = 8;
    reduced.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_OS_16(&a[i], &b[i], &product);
        BigNum_Fast_Reduction(&product, &P, &reduced);
        sink ^= reduced.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_OS_16(&a[i], &b[i], &product);
        BigNum_Fast_Reduction(&product, &P, &reduced);
        sink ^= reduced.num[0];
    }

    uint64_t end = read_cycles();
    return make_bench_result(end - start);
}

static BenchResult bench_fast_reduction_total(void)
{
    BigNum r;
    uint32_t r_buf[9];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 8;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_Fast_Reduction(&red_in[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_Fast_Reduction(&red_in[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return make_bench_result(end - start);
}

static void print_bench_result(const char* label, BenchResult result)
{
    printf("%-24s\ntotal: %llu cycles\navg: %.2f cycles/op\n",
           label, (unsigned long long)result.total, result.avg);
}

static double bench_ps_32(void)
{
    BigNum r;
    uint32_t r_buf[16];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 16;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_PS_32(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_PS_32(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_ps(void)
{
    BigNum r;
    uint32_t r_buf[16];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 16;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_PS(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_PS(&a[i], &b[i], &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_squaring(void)
{
    BigNum r;
    uint32_t r_buf[16];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 16;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_squaring(&a[i], &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_squaring(&a[i], &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_reduction(void)
{
    BigNum r;
    uint32_t r_buf[9];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 8;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_Fast_Reduction(&red_in[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_Fast_Reduction(&red_in[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_inverse(void)
{
    BigNum r;
    uint32_t r_buf[8];
    volatile uint32_t sink = 0;

    r.num = r_buf;
    r.size = 8;
    r.sign = 1;

    for (int i = 0; i < WARMUP; i++) {
        BigNum_Inverse(&a[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        BigNum_Inverse(&a[i], &P, &r);
        sink ^= r.num[0];
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static void free_jacobian_result(Point_Jacobian* point)
{
    free(point->x.num);
    free(point->y.num);
    free(point->z.num);
}

static double bench_point_double_jacobian(void)
{
    volatile uint32_t sink = 0;

    for (int i = 0; i < WARMUP; i++) {
        Point_Jacobian r = {0};
        point_double_jacobian(&jacobian_points[i], &r, &P);
        sink ^= r.x.num[0] ^ r.y.num[0] ^ r.z.num[0];
        free_jacobian_result(&r);
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        Point_Jacobian r = {0};
        point_double_jacobian(&jacobian_points[i], &r, &P);
        sink ^= r.x.num[0] ^ r.y.num[0] ^ r.z.num[0];
        free_jacobian_result(&r);
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_point_add_jacobian(void)
{
    volatile uint32_t sink = 0;

    for (int i = 0; i < WARMUP; i++) {
        Point_Jacobian r = {0};
        point_add_jacobian(&jacobian_points[i], &affine_points_q[i], &r, &P);
        sink ^= r.x.num[0] ^ r.y.num[0] ^ r.z.num[0];
        free_jacobian_result(&r);
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        Point_Jacobian r = {0};
        point_add_jacobian(&jacobian_points[i], &affine_points_q[i], &r, &P);
        sink ^= r.x.num[0] ^ r.y.num[0] ^ r.z.num[0];
        free_jacobian_result(&r);
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

static double bench_point_add(void)
{
    volatile uint32_t sink = 0;

    for (int i = 0; i < WARMUP; i++) {
        Point r = {0};
        uint32_t rx[9] = {0};
        uint32_t ry[9] = {0};
        r.x.num = rx;
        r.y.num = ry;
        r.x.size = 8;
        r.y.size = 8;
        r.x.sign = 1;
        r.y.sign = 1;

        point_add(&affine_points_p[i], &affine_points_q[i], &r, &P);
        sink ^= r.x.num[0] ^ r.y.num[0] ^ r.infinity;
    }

    uint64_t start = read_cycles();

    for (int i = 0; i < N; i++) {
        Point r = {0};
        uint32_t rx[9] = {0};
        uint32_t ry[9] = {0};
        r.x.num = rx;
        r.y.num = ry;
        r.x.size = 8;
        r.y.size = 8;
        r.x.sign = 1;
        r.y.sign = 1;

        point_add(&affine_points_p[i], &affine_points_q[i], &r, &P);
        sink ^= r.x.num[0] ^ r.y.num[0] ^ r.infinity;
    }

    uint64_t end = read_cycles();
    return (double)(end - start) / N;
}

int main(void)
{
    for (int i = 0; i < N; i++) {
        a[i].num = a_buf[i];
        b[i].num = b_buf[i];

        random_bn_256(&a[i]);
        random_bn_256(&b[i]);
    }

    for (int i = 0; i < N; i++) {
        jacobian_points[i].x.num = jacobian_x_buf[i];
        jacobian_points[i].y.num = jacobian_y_buf[i];
        jacobian_points[i].z.num = jacobian_z_buf[i];
        random_bn_256(&jacobian_points[i].x);
        random_bn_256(&jacobian_points[i].y);
        random_bn_256(&jacobian_points[i].z);
    }

    for (int i = 0; i < N; i++) {
        affine_points_p[i].x = a[i];
        affine_points_p[i].y = b[i];
        affine_points_p[i].infinity = 0;

        affine_points_q[i].x = jacobian_points[i].x;
        affine_points_q[i].y = jacobian_points[i].y;
        affine_points_q[i].infinity = 0;

        if (BigNum_cmp(&affine_points_p[i].x, &affine_points_q[i].x) == 0)
            affine_points_q[i].x.num[0] ^= 1;
    }

    for (int i = 0; i < N; i++) {
        red_in[i].num = red_in_buf[i];
        red_in[i].size = 16;
        red_in[i].sign = 1;
        BigNum_OS_16(&a[i], &b[i], &red_in[i]);
    }

    // printf("mod_add    : %.2f cycles/op\n", bench_mod_add());
    // printf("mod_sub    : %.2f cycles/op\n", bench_mod_sub());
    // printf("OS              : %.2f cycles/op\n", bench_os());
    // printf("PS_32           : %.2f cycles/op\n", bench_ps_32());
    print_bench_result("BigNum_OS", bench_os_16());
    printf("\n");
    print_bench_result("BigNum_OS + Fast_Reduction", bench_os_16_reduction());
    printf("\n");
    print_bench_result("Fast_Reduction", bench_fast_reduction_total());
    // printf("PS              : %.2f cycles/op\n", bench_ps());
    // printf("squaring   : %.2f cycles/op\n", bench_squaring());
    // printf("reduction  : %.2f cycles/op\n", bench_reduction());
    // printf("inverse    : %.2f cycles/op\n", bench_inverse());
    // printf("double_jacobian  : %.2f cycles/op\n", bench_point_double_jacobian());
    // printf("add_jacobian     : %.2f cycles/op\n", bench_point_add_jacobian());
    // printf("point_add        : %.2f cycles/op\n", bench_point_add());

    return 0;
}
