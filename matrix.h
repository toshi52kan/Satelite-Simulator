#pragma once

#ifndef _MATRIX
#define _MATRIX

typedef struct {
	int rows;
	int columns;
	double *values;

}matrix_t;

typedef matrix_t* (*FcnEq)(matrix_t*, matrix_t*, matrix_t*);

matrix_t *matrix_init(matrix_t *target, int rows, int columns);
matrix_t *matrix_zero(matrix_t *target);
matrix_t *matrix_mul(matrix_t *a, matrix_t *b, matrix_t *result);
matrix_t *matrix_add(matrix_t *a, matrix_t *b, matrix_t *result);
matrix_t *matrix_alloc(matrix_t *target);
matrix_t *matrix_sub(matrix_t *a, matrix_t *b, matrix_t *result);
matrix_t *matrix_unit(matrix_t *c);
matrix_t *matrix_inverse(matrix_t *target);
matrix_t *matrix_transpose(matrix_t *target,matrix_t *result);
matrix_t *matrix_scalar(matrix_t *a, double b, matrix_t *c);
matrix_t *matrix_values_copy(matrix_t *target, matrix_t *result,int rows,int columns);




#endif