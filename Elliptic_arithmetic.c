#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ECC_BigNumOperation.h"
#include "Elliptic_arithmetic.h"

int BigNum_Is_Zero(BigNum* a)
{
    for(int i = 0; i < a->size; i++)
    {
        if(a->num[i] != 0)
            return 0;
    }
    return 1;
}

void mod_half(BigNum* a, BigNum* p, BigNum* result)
{
    BigNum tmp;
    uint32_t tmp_buf[9];

    tmp.num = tmp_buf;
    tmp.size = 8;
    tmp.sign = 1;

    if ((a->num[0] & 1) == 0) {
        BigNum_copy(result, a);
        BigNum_right_shift1(result);
    } else {
        BigNum_add(a, p, &tmp);
        BigNum_right_shift1(&tmp);
        BigNum_copy(result, &tmp);
    }

    result->size = 8;
    result->sign = 1;
}

void BigNum_mod(BigNum* a, BigNum* b, BigNum* result, BigNum* p)
{
    BigNum tmp16;

    tmp16.size = 16;
    tmp16.sign = 1;
    tmp16.num = calloc(16, sizeof(uint32_t));

    BigNum_PS(a, b, &tmp16);
    BigNum_Fast_Reduction(&tmp16, p, result);
    
    free(tmp16.num);
}

void BigNum_sqr(BigNum* a, BigNum* result, BigNum* p)
{
    BigNum tmp16;

    tmp16.size = 16;
    tmp16.sign = 1;
    tmp16.num = calloc(16, sizeof(uint32_t));

    BigNum_squaring(a, &tmp16);
    BigNum_Fast_Reduction(&tmp16, p, result);
    
    free(tmp16.num);
}

void affine_to_jacobian(Point* affine, Point_Jacobian* jacobian)
{
    jacobian->x.size = 8;
    jacobian->y.size = 8;
    jacobian->z.size = 8;

    jacobian->x.num = calloc(8, sizeof(uint32_t));
    jacobian->y.num = calloc(8, sizeof(uint32_t));
    jacobian->z.num = calloc(8, sizeof(uint32_t));

    jacobian->x.sign = 1;
    jacobian->y.sign = 1;
    jacobian->z.sign = 1;

    if (affine->infinity)
    {
        jacobian->x.num[0] = 0;
        jacobian->y.num[0] = 0;
        jacobian->z.num[0] = 0;
    }
    else
    {
        for (int i = 0; i < affine->x.size && i < 8; i++)
            jacobian->x.num[i] = affine->x.num[i];

        for (int i = 0; i < affine->y.size && i < 8; i++)
            jacobian->y.num[i] = affine->y.num[i];

        jacobian->z.num[0] = 1;

        jacobian->x.sign = affine->x.sign;
        jacobian->y.sign = affine->y.sign;
    }
}

void jacobian_to_affine(Point_Jacobian* jacobian, Point* affine, BigNum* p)
{
    if(BigNum_Is_Zero(&jacobian->z))
    {
        affine->x.size = 1;
        affine->y.size = 1;
        affine->x.num = calloc(1, sizeof(uint32_t));
        affine->y.num = calloc(1, sizeof(uint32_t));
        affine->x.num[0] = 0;
        affine->y.num[0] = 0;
        affine->x.sign = 1;
        affine->y.sign = 1;
        affine->infinity = 1;
    }
    else
    {
        BigNum z_inv, z2_inv, z3_inv, tmp16;

        z_inv.size = 8;
        z2_inv.size = 8;
        z3_inv.size = 8;
        tmp16.size = 16;
        z_inv.sign = 1;
        z2_inv.sign = 1;
        z3_inv.sign = 1;
        tmp16.sign = 1;
        z_inv.num = calloc(8, sizeof(uint32_t));
        z2_inv.num = calloc(8, sizeof(uint32_t));
        z3_inv.num = calloc(8, sizeof(uint32_t));
        tmp16.num = calloc(16, sizeof(uint32_t));

        affine->x.size = 8;
        affine->y.size = 8;
        affine->x.num = calloc(8, sizeof(uint32_t));
        affine->y.num = calloc(8, sizeof(uint32_t));
        affine->x.sign = 1;
        affine->y.sign = 1;

        BigNum_Inverse(&jacobian->z, p, &z_inv);

        BigNum_squaring(&z_inv, &tmp16);
        BigNum_Fast_Reduction(&tmp16, p, &z2_inv);
        
        BigNum_PS(&z2_inv, &z_inv, &tmp16);
        BigNum_Fast_Reduction(&tmp16, p, &z3_inv);
        
        BigNum_PS(&jacobian->x, &z2_inv, &tmp16);
        BigNum_Fast_Reduction(&tmp16, p, &affine->x);

        BigNum_PS(&jacobian->y, &z3_inv, &tmp16);
        BigNum_Fast_Reduction(&tmp16, p, &affine->y);

        affine->infinity = 0;

        free(z_inv.num);
        free(z2_inv.num);
        free(z3_inv.num);
        free(tmp16.num);
    }
}

void point_double_jacobian(Point_Jacobian* P, Point_Jacobian* result, BigNum* p)    
{
    BigNum T1, T2, T3;

    T1.size = 8; T2.size = 8; T3.size = 8;
    T1.sign = 1; T2.sign = 1; T3.sign = 1;
    T1.num = calloc(9, sizeof(uint32_t));
    T2.num = calloc(9, sizeof(uint32_t));
    T3.num = calloc(9, sizeof(uint32_t));

    result->x.size = 8;
    result->y.size = 8;
    result->z.size = 8;

    result->x.sign = 1;
    result->y.sign = 1;
    result->z.sign = 1;

    result->x.num = calloc(9, sizeof(uint32_t));
    result->y.num = calloc(9, sizeof(uint32_t));
    result->z.num = calloc(9, sizeof(uint32_t));

    if(BigNum_Is_Zero(&P->z) == 0)
    {
        BigNum_sqr(&P->z, &T1, p);     
        mod_sub(&P->x, &T1, p, &T2);       
        mod_add(&P->x, &T1, p, &T1);        
        BigNum_mod(&T1, &T2, &T2, p);     
        BigNum_copy(&T1, &T2);
        mod_add(&T2, &T2, p, &T2);         
        mod_add(&T2, &T1, p, &T2);
        mod_add(&P->y, &P->y, p, &result->y);       
        BigNum_mod(&result->y, &P->z, &result->z, p);    
        BigNum_sqr(&result->y, &result->y, p);       
        BigNum_mod(&result->y, &P->x, &T3, p);        
        BigNum_sqr(&result->y, &result->y, p);       
        mod_half(&result->y, p, &result->y);       
        BigNum_sqr(&T2, &result->x, p);       
        mod_add(&T3, &T3, p, &T1);
        mod_sub(&result->x, &T1, p, &result->x);       
        mod_sub(&T3, &result->x, p, &T1);        
        BigNum_mod(&T1, &T2, &T1, p);        
        mod_sub(&T1, &result->y, p, &result->y);       
    }
    
    free(T1.num);
    free(T2.num);
    free(T3.num);
}

void point_add_jacobian(Point_Jacobian* P, Point* Q, Point_Jacobian* result, BigNum* p)
{
    BigNum T1, T2, T3, T4;

    T1.size = 8; T2.size = 8; T3.size = 8; T4.size = 8;
    T1.sign = 1; T2.sign = 1; T3.sign = 1; T4.sign = 1;
    T1.num = calloc(9, sizeof(uint32_t));
    T2.num = calloc(9, sizeof(uint32_t));
    T3.num = calloc(9, sizeof(uint32_t));
    T4.num = calloc(9, sizeof(uint32_t));

    result->x.size = 8;
    result->y.size = 8;
    result->z.size = 8;

    result->x.sign = 1;
    result->y.sign = 1;
    result->z.sign = 1;

    result->x.num = calloc(9, sizeof(uint32_t));
    result->y.num = calloc(9, sizeof(uint32_t));
    result->z.num = calloc(9, sizeof(uint32_t));

    if(Q->infinity)
    {
        BigNum_copy(&result->x, &P->x);
        BigNum_copy(&result->y, &P->y);
        BigNum_copy(&result->z, &P->z);
    }
    else if(BigNum_Is_Zero(&P->z))
    {
        BigNum_copy(&result->x, &Q->x);
        BigNum_copy(&result->y, &Q->y);
        result->z.num[0] = 1;
    }
    else
    {
        BigNum_sqr(&P->z, &T1, p);     
        BigNum_mod(&T1, &P->z, &T2, p);     
        BigNum_mod(&T1, &Q->x, &T1, p);     
        BigNum_mod(&T2, &Q->y, &T2, p);     
        mod_sub(&T1, &P->x, p, &T1);      
        mod_sub(&T2, &P->y, p, &T2);        
        if(BigNum_Is_Zero(&T1))
        {
            if(BigNum_Is_Zero(&T2))
            {
                point_double_jacobian(P, result, p);        // Q를 사용하지 않고, P를 doubling
                goto cleanup;
            }
            else
            {
                goto cleanup;
            }
        }
        BigNum_mod(&P->z, &T1, &result->z, p);    
        BigNum_sqr(&T1, &T3, p);    
        BigNum_mod(&T3, &T1, &T4, p);     
        BigNum_mod(&T3, &P->x, &T3, p);     
        mod_add(&T3, &T3, p, &T1);         
        BigNum_sqr(&T2, &result->x, p);     
        mod_sub(&result->x, &T1, p, &result->x);       
        mod_sub(&result->x, &T4, p, &result->x);       
        mod_sub(&T3, &result->x, p, &T3);        
        BigNum_mod(&T3, &T2, &T3, p);        
        BigNum_mod(&T4, &P->y, &T4, p);        
        mod_sub(&T3, &T4, p, &result->y);      
    }
 cleanup:
    free(T1.num);
    free(T2.num);
    free(T3.num);
    free(T4.num);
}

void scalar_mul(Point* P, BigNum* k, Point* result, BigNum* p)          // Ltr
{
    Point_Jacobian result_jacobian;

    result_jacobian.x.size = 8;
    result_jacobian.y.size = 8;
    result_jacobian.z.size = 8;

    result_jacobian.x.sign = 1;
    result_jacobian.y.sign = 1;
    result_jacobian.z.sign = 1;

    result_jacobian.x.num = calloc(9, sizeof(uint32_t));
    result_jacobian.y.num = calloc(9, sizeof(uint32_t));
    result_jacobian.z.num = calloc(9, sizeof(uint32_t));

    for (int i = k->size * 32 - 1; i >= 0; i--)
    {
        Point_Jacobian tmp = {0};

        point_double_jacobian(&result_jacobian, &tmp, p);
        BigNum_copy(&result_jacobian.x, &tmp.x);
        BigNum_copy(&result_jacobian.y, &tmp.y);
        BigNum_copy(&result_jacobian.z, &tmp.z);
        free(tmp.x.num);
        free(tmp.y.num);
        free(tmp.z.num);

        if ((k->num[i / 32] >> (i % 32)) & 1)
        {
            Point_Jacobian added = {0};

            point_add_jacobian(&result_jacobian, P, &added, p);
            BigNum_copy(&result_jacobian.x, &added.x);
            BigNum_copy(&result_jacobian.y, &added.y);
            BigNum_copy(&result_jacobian.z, &added.z);
            free(added.x.num);
            free(added.y.num);
            free(added.z.num);
        }
    }

    jacobian_to_affine(&result_jacobian, result, p);

    free(result_jacobian.x.num);
    free(result_jacobian.y.num);
    free(result_jacobian.z.num);
}

void point_add(Point* P, Point* Q, Point* result, BigNum* p)
{
    if (P->infinity)
    {
        BigNum_copy(&result->x, &Q->x);
        BigNum_copy(&result->y, &Q->y);
        result->infinity = Q->infinity;
    }
    else if (Q->infinity)
    {
        BigNum_copy(&result->x, &P->x);
        BigNum_copy(&result->y, &P->y);
        result->infinity = P->infinity;
    }
    else if (BigNum_cmp(&P->x, &Q->x) == 0)
    {
        if (BigNum_cmp(&P->y, &Q->y) == 0)
        {
            Point_Jacobian P_jacobian, R_jacobian;

            affine_to_jacobian(P, &P_jacobian);
            point_double_jacobian(&P_jacobian, &R_jacobian, p);
            jacobian_to_affine(&R_jacobian, result, p);

            free(P_jacobian.x.num);
            free(P_jacobian.y.num);
            free(P_jacobian.z.num);
            free(R_jacobian.x.num);
            free(R_jacobian.y.num);
            free(R_jacobian.z.num);
        }
        else
        {
            result->x.size = 1;
            result->y.size = 1;
            result->x.num = calloc(1, sizeof(uint32_t));
            result->y.num = calloc(1, sizeof(uint32_t));
            result->x.num[0] = 0;
            result->y.num[0] = 0;
            result->x.sign = 1;
            result->y.sign = 1;
            result->infinity = 1;
        }
    }
    else
    {
        BigNum lambda, tmp1, tmp2;

        lambda.size = 8; tmp1.size = 8; tmp2.size = 8;
        lambda.sign = 1; tmp1.sign = 1; tmp2.sign = 1;

        lambda.num = calloc(9, sizeof(uint32_t));
        tmp1.num = calloc(9, sizeof(uint32_t));
        tmp2.num = calloc(9, sizeof(uint32_t));

        BigNum inv;

        inv.size = 8; inv.sign = 1;
        inv.num = calloc(9, sizeof(uint32_t));

        mod_sub(&Q->y, &P->y, p, &tmp1);
        mod_sub(&Q->x, &P->x, p, &tmp2);
        BigNum_Inverse(&tmp2, p, &inv);
        BigNum_mod(&tmp1, &inv, &lambda, p);
        BigNum_sqr(&lambda, &tmp1, p);
        mod_sub(&tmp1, &P->x, p, &tmp1);
        mod_sub(&tmp1, &Q->x, p, &result->x);
        mod_sub(&P->x, &result->x, p, &tmp1);
        BigNum_mod(&tmp1, &lambda, &tmp2, p);
        mod_sub(&tmp2, &P->y, p, &result->y);

        result->infinity = 0;

        free(lambda.num);
        free(tmp1.num);
        free(tmp2.num);
        free(inv.num);
    }
}