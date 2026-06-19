#ifndef ELLIPTIC_ARITHMETIC_H
#define ELLIPTIC_ARITHMETIC_H

#include "ECC_BigNumOperation.h"

typedef struct Point {
    BigNum x;
    BigNum y;
    uint32_t infinity;
} Point;

typedef struct Point_Jacobian {
    BigNum x;
    BigNum y;
    BigNum z;
} Point_Jacobian;

int BigNum_Is_Zero(BigNum* a);
void mod_half(BigNum* a, BigNum* p, BigNum* result);
void BigNum_mod(BigNum* a, BigNum* b, BigNum* result, BigNum* p);
void BigNum_sqr(BigNum* a, BigNum* result, BigNum* p);
void affine_to_jacobian(Point* affine, Point_Jacobian* jacobian);
void jacobian_to_affine(Point_Jacobian* jacobian, Point* affine, BigNum* p);
void point_double_jacobian(Point_Jacobian* P, Point_Jacobian* result, BigNum* p);
void point_add_jacobian(Point_Jacobian* P, Point* Q, Point_Jacobian* result, BigNum* p);
void scalar_mul(Point* P, BigNum* k, Point* result, BigNum* p);
void point_add(Point* P, Point* Q, Point* result, BigNum* p);

#endif