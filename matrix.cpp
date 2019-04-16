#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <conio.h>
#include <math.h>
#include "matrix.h"

/**
* matrix_t型を初期化します。
*
* @param target 対象とする行列
*/
matrix_t *matrix_init(matrix_t *target,int rows,int columns ) {
	target = (matrix_t *)malloc(sizeof(matrix_t));
	target->values = NULL;
	target->rows = rows;
	target->columns = columns;

	return target;
}

/**
* 行列の要素を0で初期化します。
*
* @param target 対象とする行列
* @return エラーがおきた際には SUCCESS 以外の値が返ります
*/

matrix_t *matrix_zero(matrix_t *target) {
	int i;
	for (i = 0; i < target->rows*target->columns; i++) {
		target->values[i] = 0;
	}
	return target;
}

matrix_t *matrix_mul(matrix_t *a, matrix_t *b, matrix_t *result) {
	int i, j, k;
	matrix_t cpy;
	cpy.rows = result->rows;
	cpy.columns = result->columns;
	cpy.values = (double *)malloc(sizeof(double)*result->rows*result->columns);
	//	cpy = *result;
	for (int i = 0; i < result->rows*result->columns; i++) {
		cpy.values[i] = result->values[i];
	}

	for (i = 0; i<result->rows; i++) {
		for (j = 0; j<result->columns; j++) {
			cpy.values[i * result->columns + j] = 0.0;
			for (k = 0; k<a->columns; k++) cpy.values[i*result->columns + j] += a->values[i*a->columns + k] * b->values[k*b->columns + j];
		//				printf("res:%f\n", cpy.values[i*result->columns + j]);
			
		}
	}

	for (int i = 0; i < result->rows*result->columns; i++) {
		result->values[i] = cpy.values[i];
	}

	free(cpy.values);
	return result;
}

matrix_t *matrix_add(matrix_t *a, matrix_t *b, matrix_t *result) {
	int i, j, k;
	matrix_t cpy;
	cpy.rows = result->rows;
	cpy.columns = result->columns;
	cpy.values = (double *)malloc(sizeof(double)*result->rows*result->columns);

	for (int i = 0; i < result->rows*result->columns; i++) {
		cpy.values[i] = result->values[i];
//				printf("k2:%f\n", result->values[i]);
	}
	for (i = 0; i<result->rows; i++) {
		for (j = 0; j<result->columns; j++) {
			cpy.values[i * result->columns + j] = 0.0;
			cpy.values[i*result->columns + j] = a->values[i*result->columns + j] + b->values[i*result->columns + j];
		}
	}
	for (int i = 0; i < result->rows*result->columns; i++) {
		result->values[i] = cpy.values[i];
//		printf("result:%f\n", cpy.values[i]);
	}
	free(cpy.values);
	return result;
}

matrix_t *matrix_alloc(matrix_t *target) {

	target->values = (double *)malloc(sizeof(double)*target->rows*target->columns);

	return target;
}

/*matrix_t *matrix_add(matrix_t *a, matrix_t *b) {
	int i, j, k;
	matrix_t *result = 0;
	result = matrix_init(result, a->columns, b->rows);
	result = matrix_alloc(result);

	for (i = 0; i<result->rows; i++) {
		for (j = 0; j<result->columns; j++) {
			result->values[i * result->columns + j] = 0.0;
			result->values[i*result->columns + j] = a->values[i*result->columns + j] + b->values[i*result->columns + j];
			printf("result:%f\n", result->values[i]);
		}
	}

	
	return result;
}*/

matrix_t *matrix_sub(matrix_t *a, matrix_t *b, matrix_t *result) {
	int i, j, k;
	matrix_t cpy;
	cpy.rows = result->rows;
	cpy.columns = result->columns;
	cpy.values = (double *)malloc(sizeof(double)*result->rows*result->columns);

	for (int i = 0; i < result->rows*result->columns; i++) {
		cpy.values[i] = result->values[i];
	}

	for (i = 0; i<result->rows; i++) {
		for (j = 0; j<result->columns; j++) {
			cpy.values[i * result->columns + j] = 0.0;
			cpy.values[i*result->columns + j] = a->values[i*result->columns + j] - b->values[i*result->columns + j];
		}
	}

	for (int i = 0; i < result->rows*result->columns; i++) {
		result->values[i] = cpy.values[i];
	}
	free(cpy.values);
	return result;
}

matrix_t *matrix_unit(matrix_t *c) {
	int i;
//	c->values = 0;
	for (i = 0; i<c->rows; i++) {
		c->values[i*c->columns + i] = 1;
	}
	return c;
}



matrix_t *matrix_scalar(matrix_t *a, double b, matrix_t *result) {
	matrix_t cpy;
	cpy.rows = result->rows;
	cpy.columns = result->columns;
	cpy.values = (double *)malloc(sizeof(double)*result->rows*result->columns);

	for (int i = 0; i < result->rows*result->columns; i++) {
		cpy.values[i] = result->values[i];
	}

	for (int i = 0; i<result->rows; i++) {
		for (int j = 0; j < result->columns; j++) {
			cpy.values[i * result->columns + j] = 0.0;
			cpy.values[i*result->columns + j] = a->values[i*result->columns + j] * b;
		}
	}

	for (int i = 0; i < result->rows*result->columns; i++) {
		result->values[i] = cpy.values[i];
	}
	free(cpy.values);

	return result;
}

matrix_t *matrix_values_copy(matrix_t *target, matrix_t *result,int rows,int columns) {
	for (int i = 0; i < rows*columns; i++) {
		result->values[i] = target->values[i];
	}
	return result;
}


matrix_t *matrix_inverse(matrix_t *target) {
	int pivotrow = 0;
	int init_pivotrow = 0;
	double max;
	double tmp;
	matrix_t *L = 0;
	L = matrix_init(L, 3, 3);
	L = matrix_alloc(L);
	L = matrix_zero(L);

	matrix_t *U = 0;
	U = matrix_init(U, 3, 3);
	U = matrix_alloc(U);
	U = matrix_zero(U);

/*	matrix_t cpy;
	cpy.rows = target->rows;
	cpy.columns = target->columns;
	cpy.values = (double *)malloc(sizeof(double)*target->rows*target->columns);

	for (int i = 0; i < target->rows*target->columns; i++) {
		cpy.values[i] = target->values[i];
	}*/

	max = fabs(target->values[0]);
	for (int i = 0; i < target->rows; i++) {
		pivotrow = i;
		if (i == 0) {
			for (int j = i + 1; j < target->rows; j++) {
				if (max < fabs(target->values[j*target->columns])) {
					max = fabs(target->values[j*target->columns]);
					pivotrow = j;
				}
			}

			if (init_pivotrow != pivotrow) {
				for (int j = 0; j < target->rows; j++) {
					tmp = target->values[init_pivotrow*target->columns + j];
					target->values[init_pivotrow*target->columns + j] = target->values[pivotrow*target->columns + j];
					target->values[pivotrow*target->columns + j] = tmp;
				}

			}
		}

		for (int j = i+1; j<target->rows; j++) {
			target->values[j * target->rows + i] /= target->values[i * target->rows + i];
			for (int k = i+1; k < target->columns; k++) {
				target->values[j*target->columns + k] -= target->values[j*target->columns + i] * target->values[i*target->columns + k];
			}
				
		}


	}
	for (int j = 0; j < 9; j++) {
		printf("LU:%f\n", target->values[j]);
	}
	L->values[0] = 1.0;
	L->values[1] = 0.0;
	L->values[2] = 0.0;
	L->values[3] = -target->values[3];
	L->values[4] = 1.0;
	L->values[5] = 0.0;
	L->values[6] = target->values[3] * target->values[7] - target->values[6];
	L->values[7] = -target->values[7];
	L->values[8] = 1.0;

	U->values[0] = target->values[4] * target->values[8];
	U->values[1] = -target->values[1] * target->values[8];
	U->values[2] = target->values[1] * target->values[5] - target->values[2]*target->values[4];
	U->values[3] = 0.0;
	U->values[4] = target->values[0] * target->values[8];
	U->values[5] = -target->values[0] * target->values[5];
	U->values[6] = 0.0;
	U->values[7] = 0.0;
	U->values[8] = target->values[0] * target->values[4];

/*	U->values[0] = 1.0/ target->values[0];
	U->values[1] = -target->values[1] / (target->values[0] * target->values[4]);
	U->values[2] = (target->values[1] * target->values[5]) - (target->values[2] * target->values[4]) / (target->values[0] * target->values[4] * target->values[8]);
	U->values[3] = 0.0;
	U->values[4] = 1.0 / target->values[4];
	U->values[5] = target->values[5] / (target->values[0] * target->values[4]);
	U->values[6] = 0.0;
	U->values[7] = 0.0;
	U->values[8] = 1.0 / target->values[8];*/
	for (int i = 0; i < target->rows*target->columns; i++) {
		U->values[i] = U->values[i] / (target->values[0] * target->values[4] * target->values[8]);
	}

	matrix_mul(U, L, target);
	for (int j = 0; j < 9; j++) {
		printf("swap:%f\n", target->values[j]);
	}
	for (int j = 0; j < 9; j++) {
		printf("L:%f\n", L->values[j]);
		printf("U:%f\n", U->values[j]);
	}

	free(L);
	free(U);
	return target;

}


matrix_t *matrix_transpose(matrix_t *target,matrix_t *result) {
	matrix_t cpy;
	cpy.rows = target->rows;
	cpy.columns = target->columns;
	cpy.values = (double *)malloc(sizeof(double)*target->rows*target->columns);

	for (int i = 0; i < target->rows*target->columns; i++) {
		cpy.values[i] = target->values[i];
	}

	int buf = target->rows;
	result->rows = target->columns;
	result->columns = buf;

	for (int i = 0; i<result->rows; i++) {
		for (int j = 0; j<result->columns; j++) {
			result->values[i*result->columns + j] = cpy.values[j*result->rows + i];
//			printf("ret:%f\n", target->values[j*target->columns + i]);
		}
	}

	free(cpy.values);
	return result;
}
