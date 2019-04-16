#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <conio.h>
#include <math.h>
#include "matrix.h"
#include "variable.h"

#define n 3
#define ROW 3
#define COLUMN 3

#define TIMESTEP 0.1
#define VAR_STT_X sin(7.0 / 3.0 / 3600.0 / 2.0 / 180.0 * M_PI)
#define VAR_STT_Y sin(7.0 / 3.0 / 3600.0 / 2.0 / 180.0 * M_PI)
#define VAR_STT_Z sin(77.0 / 3.0 / 3600.0 / 2.0 / 180.0 * M_PI)
#define RATE_STT 1
#define SYSTEM_NOISE 0.0023*M_PI/180
#define OBSERVE_NOISE 0.0000005*M_PI/180
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */

						/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti<N; mti++) {
		mt[mti] =
			(1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		/* 2002/01/09 modified by Makoto Matsumoto             */
		mt[mti] &= 0xffffffffUL;
		/* for >32 bit machines */
	}
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
	int i, j, k;
	init_genrand(19650218UL);
	i = 1; j = 0;
	k = (N>key_length ? N : key_length);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
			+ init_key[j] + j; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++; j++;
		if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
		if (j >= key_length) j = 0;
	}
	for (k = N - 1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
			- i; /* non linear */
		mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
		i++;
		if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
	}

	mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N) { /* generate N words at one time */
		int kk;

		if (mti == N + 1)   /* if init_genrand() has not been called, */
			init_genrand(5489UL); /* a default initial seed is used */

		for (kk = 0; kk<N - M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk<N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
	return (long)(genrand_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
	return genrand_int32()*(1.0 / 4294967295.0);
	/* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
	return genrand_int32()*(1.0 / 4294967296.0);
	/* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
	return (((double)genrand_int32()) + 0.5)*(1.0 / 4294967296.0);
	/* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
	unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
	return(a*67108864.0 + b)*(1.0 / 9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */


double Uniform(void) {
	return genrand_real3();
}
double rand_normal(double mu, double sigma) {
	double z = sqrt(-2.0*log(Uniform())) * sin(2.0*M_PI*Uniform());
	return mu + sigma*z;
}

matrix_t *SkewSymMat(matrix_t *a, matrix_t *b) {

	b->values[0] = 0.0;
	b->values[1] = -a->values[2];
	b->values[2] = a->values[1];
	b->values[3] = a->values[2];
	b->values[4] = 0.0;
	b->values[5] = -a->values[0];
	b->values[6] = -a->values[1];
	b->values[7] = a->values[0];
	b->values[8] = 0.0;

	return b;
}

matrix_t *q_matrix(matrix_t *x, int sw, matrix_t *result) {

	if (sw == 1) {
		result->values[0] = x->values[0];
		result->values[1] = -x->values[3];
		result->values[2] = x->values[2];
		result->values[3] = x->values[3];
		result->values[4] = x->values[0];
		result->values[5] = -x->values[1];
		result->values[6] = -x->values[2];
		result->values[7] = x->values[1];
		result->values[8] = x->values[0];
		result->values[9] = -x->values[1];
		result->values[10] = -x->values[2];
		result->values[11] = -x->values[3];
	}

	if (sw == 2) {
		result->values[0] = x->values[0];
		result->values[1] = -x->values[3];
		result->values[2] = x->values[2];
		result->values[3] = x->values[1];
		result->values[4] = x->values[3];
		result->values[5] = x->values[0];
		result->values[6] = -x->values[1];
		result->values[7] = x->values[2];
		result->values[8] = -x->values[2];
		result->values[9] = x->values[1];
		result->values[10] = x->values[0];
		result->values[11] = x->values[3];
		result->values[12] = -x->values[1];
		result->values[13] = -x->values[2];
		result->values[14] = -x->values[3];
		result->values[15] = -x->values[0];
	}
	return result;
}

matrix_t *Kinematics(matrix_t *quat, matrix_t *omega, matrix_t *result) {

	matrix_t *ret = 0;
	double form = 0;
	ret = matrix_init(ret, 4, 3);
	ret = matrix_alloc(ret);
	ret = matrix_zero(ret);

	q_matrix(quat, 1, ret);

	matrix_scalar(ret, 0.5, ret);
	/*	for (int i = 0; i < 12; i++) {
	printf("ret:%f\n", ret->values[i]);
	}*/
	matrix_mul(ret, omega, result);
	//	matrix_add(ret, quat, result);
	form = 1 - ((result->values[1] * result->values[1]) + (result->values[2] * result->values[2]) + (result->values[3] * result->values[3]));

	result->values[0] = sqrt(form);
	//	printf("form:%f\n", result->values[0]);
	free(ret);

	return result;
}


matrix_t *DynEq(matrix_t *omega, matrix_t *torque, matrix_t *result) {

	matrix_t *hoge = 0;
	hoge = (matrix_t *)malloc(sizeof(matrix_t));
	hoge->rows = 3;
	hoge->columns = 3;
	hoge->values = (double *)malloc(sizeof(double)*hoge->rows*hoge->columns);

	matrix_t *ada = 0;
	ada = (matrix_t *)malloc(sizeof(matrix_t));
	ada->rows = 3;
	ada->columns = 1;
	ada->values = (double *)malloc(sizeof(double)*ada->rows*ada->columns);

	matrix_t *inv_I = 0;
	inv_I = (matrix_t *)malloc(sizeof(matrix_t));
	inv_I->rows = 3;
	inv_I->columns = 3;
	inv_I->values = (double *)malloc(sizeof(double)*inv_I->rows*inv_I->columns);

	/*	inv_I->values[0] = 1;
	inv_I->values[1] = -1;
	inv_I->values[2] = 0;
	inv_I->values[3] = 0;
	inv_I->values[4] = 1;
	inv_I->values[5] = -1;
	inv_I->values[6] = 0;
	inv_I->values[7] = 0;
	inv_I->values[8] = 1;*/

	inv_I->values[0] = 1;
	inv_I->values[1] = 0;
	inv_I->values[2] = 0;
	inv_I->values[3] = 0;
	inv_I->values[4] = 1;
	inv_I->values[5] = 0;
	inv_I->values[6] = 0;
	inv_I->values[7] = 0;
	inv_I->values[8] = 1;
	matrix_zero(hoge);
	SkewSymMat(omega, hoge);

	matrix_zero(ada);
	matrix_mul(hoge, I, hoge);
	matrix_mul(hoge, omega, ada);

	matrix_sub(torque, ada, ada);
	for (int i = 0; i < 3; i++) {
		//		printf("ada:%f\n", ada->values[i]);
	}
	//	matrix_zero(ada);
	matrix_mul(inv_I, ada, result);


	free(inv_I);
	free(hoge);
	free(ada);

	return result;
}

matrix_t *RK4(FcnEq func, matrix_t *omega, matrix_t *torque, matrix_t *result) {

	matrix_t *arg = 0;
	arg = matrix_init(arg, result->rows, result->columns);
	arg = matrix_alloc(arg);

	matrix_t *k1 = 0;
	k1 = matrix_init(k1, result->rows, result->columns);
	k1 = matrix_alloc(k1);

	matrix_t *k2 = 0;
	k2 = matrix_init(k2, result->rows, result->columns);
	k2 = matrix_alloc(k2);

	matrix_t *k3 = 0;
	k3 = matrix_init(k3, result->rows, result->columns);
	k3 = matrix_alloc(k3);

	matrix_t *k4 = 0;
	k4 = matrix_init(k4, result->rows, result->columns);
	k4 = matrix_alloc(k4);

	func(omega, torque, k1);
	k1 = matrix_scalar(k1, TIMESTEP, k1);
	matrix_scalar(k1, 0.5, arg);
	for (int i = 0; i < 3; i++) {
		//		printf("k1:%f\n", k1->values[i]);
		//		printf("omega:%f\n",omega->values[i]);
	}
	matrix_add(omega, arg, arg);


	func(arg, torque, k2);
	matrix_zero(arg);
	k2 = matrix_scalar(k2, TIMESTEP, k2);

	matrix_scalar(k2, 0.5, arg);
	matrix_add(omega, arg, arg);

	func(arg, torque, k3);
	matrix_zero(arg);
	k3 = matrix_scalar(k3, TIMESTEP, k3);
	for (int i = 0; i < 3; i++) {
		//		printf("k3[%d]:%f\n",i, k3->values[i]);
		//		printf("omega:%f\n", omega->values[i]);
	}
	matrix_scalar(k3, 1.0, arg);
	matrix_add(omega, arg, arg);

	func(arg, torque, k4);
	matrix_zero(arg);
	k4 = matrix_scalar(k4, TIMESTEP, k4);

	matrix_zero(arg);
	matrix_add(k1, matrix_scalar(k2, 2.0, arg), arg);
	matrix_add(arg, k3, arg);
	matrix_add(arg, k3, arg);
	matrix_add(arg, k4, arg);
	matrix_scalar(arg, 1.0 / 6.0, arg);
	matrix_add(arg, omega, result);
	/*	for (int i = 0; i < 3; i++) {
	printf("result:%f\n", result->values[i]);
	}*/

	free(arg);
	free(k1);
	free(k2);
	free(k3);
	free(k4);

	return omega;

}

matrix_t *PIDcontrol(matrix_t *input, matrix_t *result) {
	double kp, kd;
	kp = 2;
	matrix_t *target = 0;
	target = matrix_init(target, 3, 1);
	target = matrix_alloc(target);
	target = matrix_zero(target);

	matrix_t *ret = 0;
	ret = matrix_init(ret, 3, 1);
	ret = matrix_alloc(ret);
	ret = matrix_zero(ret);

	matrix_sub(target, input, ret);

	matrix_scalar(ret, kp, result);
	for (int i = 0; i < 3; i++) {
		//		printf("ret:%f\n", result->values[i]);
	}
	free(target);
	free(ret);

	return result;
}

matrix_t *ExtendKalmanFilter(FcnEq func, matrix_t *quat, matrix_t *omega, matrix_t *result) {

	double form = 0.0;
	func(quat, omega, quat); //状態量の時間更新


	matrix_t *ret = 0;
	ret = matrix_init(ret, 6, 6);
	ret = matrix_alloc(ret);
	ret = matrix_zero(ret);

	matrix_t *ret2 = 0;
	ret2 = matrix_init(ret2, 3, 6);
	ret2 = matrix_alloc(ret2);
	ret2 = matrix_zero(ret2);

	matrix_t *ret3 = 0;
	ret3 = matrix_init(ret3, 6, 3);
	ret3 = matrix_alloc(ret3);
	ret3 = matrix_zero(ret3);

	matrix_t *tmp = 0;
	tmp = matrix_init(tmp, 6, 6);
	tmp = matrix_alloc(tmp);
	tmp = matrix_zero(tmp);

	matrix_t *hoge = 0;
	hoge = matrix_init(hoge, 6, 6);
	hoge = matrix_alloc(hoge);
	hoge = matrix_zero(hoge);

	matrix_t *ada = 0;
	ada = matrix_init(ada, 3, 3);
	ada = matrix_alloc(ada);
	ada = matrix_zero(ada);

	matrix_t *scala = 0;
	scala = matrix_init(scala, 4, 4);
	scala = matrix_alloc(scala);
	scala = matrix_zero(scala);

	matrix_t *X_delta = 0;
	X_delta = matrix_init(X_delta, 6, 1);
	X_delta = matrix_alloc(X_delta);
	X_delta = matrix_zero(X_delta);

	matrix_t *C = 0;
	C = matrix_init(C, 3, 6);
	C = matrix_alloc(C);
	C = matrix_zero(C);

	C->values[0] = 1.0;
	C->values[7] = 1.0;
	C->values[14] = 1.0;

	matrix_t *CT = 0;
	CT = matrix_init(CT, 3, 6);
	CT = matrix_alloc(CT);
	CT = matrix_zero(CT);

	CT->values[0] = 1.0;
	CT->values[7] = 1.0;
	CT->values[14] = 1.0;

	matrix_t *AdT = 0;
	AdT = matrix_init(AdT, 6, 6);
	AdT = matrix_alloc(AdT);
	AdT = matrix_zero(AdT);

	matrix_t *BdT = 0;
	BdT = matrix_init(BdT, 6, 6);
	BdT = matrix_alloc(BdT);
	BdT = matrix_zero(BdT);


	Ad->values[1] = omega->values[2] * TIMESTEP;
	Ad->values[2] = -omega->values[1] * TIMESTEP;
	Ad->values[3] = -0.5*TIMESTEP;
	Ad->values[6] = -omega->values[2] * TIMESTEP;
	Ad->values[8] = omega->values[0] * TIMESTEP;
	Ad->values[10] = -0.5*TIMESTEP;
	Ad->values[12] = omega->values[1] * TIMESTEP;
	Ad->values[13] = -omega->values[0] * TIMESTEP;
	Ad->values[17] = -0.5*TIMESTEP;

	STT->values[1] = quat->values[0] * rand_normal(0.0, 0.0001);
	STT->values[2] = quat->values[0] * rand_normal(0.0, 0.0001);
	STT->values[3] = quat->values[0] * rand_normal(0.0, 0.0005);
	STT->values[0] = STT->values[0] = 1 - ((STT->values[1] * STT->values[1]) + (STT->values[2] * STT->values[2]) + (STT->values[3] * STT->values[3]));

/*	for (int i = 0; i < 36; i++) {
		AdT->values[i] = Ad->values[i];
	}
	
	for (int i = 0; i < 36; i++) {
		BdT->values[i] = Bd->values[i];
	}*/

	matrix_mul(matrix_mul(matrix_transpose(Ad,AdT), P, ret), Ad, ret);
	matrix_mul(matrix_mul(matrix_transpose(Bd,BdT), Q, hoge), Bd, hoge);
	matrix_add(ret, hoge, P); //事前共分散行列の計算

/*	for (int j = 0; j < 3; j++) {
		printf("omega:%f\n", omega->values[j]);
	}*/
	matrix_mul(C, P, ret2);

	matrix_mul(ret2, matrix_transpose(C, CT), ada);

	matrix_add(ada, R, ada);
	for (int j = 0; j < 9; j++) {
		printf("ada:%f\n", ada->values[j]);
	}
	matrix_inverse(ada);

	matrix_mul(matrix_mul(P, matrix_transpose(C,CT), ret3), ada, kalmangain); //カルマンゲインの計算
	for (int j = 0; j < 18; j++) {
		printf("kalman:%f\n", ret3->values[j]);
	}
	q_matrix(quat, 2, scala);
//	matrix_mul(scala, STT, quat_delta); // 適当にquatにしてる
	matrix_mul(scala, STT, quat_delta); // 適当にquatにしてる
	X_delta->values[0] = quat_delta->values[1];
	X_delta->values[1] = quat_delta->values[2];
	X_delta->values[2] = quat_delta->values[3];

//	matrix_values_copy(quat_delta, q_delta_dummy,3,1);

	matrix_mul(kalmangain, X_delta, X_delta); //観測値との差の重みづけ
/*	for (int j = 0; j < 3; j++) {
		printf("P:%f\n", q_delta_dummy->values[j]);
		//		printf("P:%f\n", quat_delta->values[j]);
	}*/
	quat_delta->values[1] = X_delta->values[0];
	quat_delta->values[2] = X_delta->values[1];
	quat_delta->values[3] = X_delta->values[2];
//	matrix_values_copy(q_delta_dummy,quat_delta,3,1);
	gyro_bias->values[0] = X_delta->values[3];
	gyro_bias->values[1] = X_delta->values[4];
	gyro_bias->values[2] = X_delta->values[5];

	quat_delta->values[0] = 1 - ((quat_delta->values[1] * quat_delta->values[1]) + (quat_delta->values[2] * quat_delta->values[2]) + (quat_delta->values[3] * quat_delta->values[3]));

//	omega->values[0] = omega->values[0]-gyro_bias->values[0];

//	quat_delta->values[0] = sqrt(quat_delta->values[0]);

//	matrix_mul(matrix_transpose(C), quat_delta, quat_delta);

//	quat_delta->values[0] = (1 - quat_delta->values[1] * quat_delta->values[1] + quat_delta->values[2] * quat_delta->values[2] + quat_delta->values[3] * quat_delta->values[3]);

	matrix_mul(scala, quat_delta, result); //状態量の観測更新


									   //	matrix_sub(I, matrix_mul(kalmangain, C, ret),ret);

	matrix_mul(kalmangain, C, tmp);
/*	for (int j = 0; j < 36; j++) {
		printf("tmp:%f\n", tmp->values[j]);
	}*/
	matrix_sub(I_kalman, tmp, tmp);


	matrix_mul(tmp, P, P); //事後共分散行列の更新
/*	for (int j = 0; j < 36; j++) {
		printf("P:%f\n", P->values[j]);
	}*/
	free(ret);
	free(hoge);
	free(ada);
	free(scala);
	free(X_delta);
	free(tmp);
	free(ret2);
	free(ret3);
	free(C);
	free(AdT);
	free(BdT);
	free(CT);


	//結果格納用の構造体の種類がたくさん必要になる
	//構造体の値の要素を取り出してコピーする機能がない

	return result;
}

int main() {

	FILE *fd;
	matrix_t *A = 0;
	matrix_t *timestep = 0;
	matrix_t *T = 0;
	matrix_t *omega = 0;
	matrix_t *quat = 0;
	matrix_t *out_trq = 0;
	matrix_t *in_trq = 0;

	P = matrix_init(P, 6, 6);
	P = matrix_alloc(P);
	P = matrix_zero(P);
	P = matrix_unit(P);
/*	P->values[0] = 0.1;
	P->values[7] = 10;
	P->values[14] = 10;
	P->values[21] = TIMESTEP;
	P->values[28] = TIMESTEP;
	P->values[35] = TIMESTEP;*/
	matrix_scalar(P, 100.0, P);

	for (int j = 0; j < 36; j++) {
		printf("P:%f\n", P->values[j]);
	}
	I_kalman = matrix_init(I_kalman, 6, 6);
	I_kalman = matrix_alloc(I_kalman);
	I_kalman = matrix_zero(I_kalman);
	matrix_unit(I_kalman);


	Ad = matrix_init(Ad, 6, 6);
	Ad = matrix_alloc(Ad);
	Ad = matrix_zero(Ad);
	Ad = matrix_unit(Ad);


	Bd = matrix_init(Bd, 6, 6);
	Bd = matrix_alloc(Bd);
	Bd = matrix_zero(Bd);
	Bd->values[0] = -0.5*TIMESTEP;
	Bd->values[7] = -0.5*TIMESTEP;
	Bd->values[14] = -0.5*TIMESTEP;
	Bd->values[21] = TIMESTEP;
	Bd->values[28] = TIMESTEP;
	Bd->values[35] = TIMESTEP;

	Q = matrix_init(Q, 6, 6);
	Q = matrix_alloc(Q);
	Q = matrix_zero(Q);
	/*	Q->values[0] = SYSTEM_NOISE*SYSTEM_NOISE*TIMESTEP;
	Q->values[7] = SYSTEM_NOISE*SYSTEM_NOISE *TIMESTEP;
	Q->values[14] = SYSTEM_NOISE*SYSTEM_NOISE*TIMESTEP;
	Q->values[21] = OBSERVE_NOISE*OBSERVE_NOISE*TIMESTEP;
	Q->values[28] = OBSERVE_NOISE*OBSERVE_NOISE*TIMESTEP;
	Q->values[35] = OBSERVE_NOISE*OBSERVE_NOISE*TIMESTEP;*/
/*	Q->values[0] = 1.1;
	Q->values[7] =1.1;
	Q->values[14] = 1.1;
	Q->values[21] = 1.1;
	Q->values[28] = 1.1;
	Q->values[35] = 1.1;*/

	Q->values[0] = 0.0001;
	Q->values[7] = 0.0001;
	Q->values[14] = 0.0001;
	Q->values[21] = 0.000001;
	Q->values[28] = 0.000001;
	Q->values[35] = 0.000001;



	R = matrix_init(R, 3, 3);
	R = matrix_alloc(R);
	R = matrix_zero(R);
	/*	R->values[0] = VAR_STT_X*VAR_STT_X*RATE_STT;
	R->values[4] = VAR_STT_Y*VAR_STT_Y*RATE_STT;
	R->values[8] = VAR_STT_Z*VAR_STT_Z*RATE_STT;*/

/*	R->values[0] = 1.1;
	R->values[4] = 1.1;
	R->values[8] = 1.1;*/

	R->values[0] = 0.000001;
	R->values[4] = 0.000001;
	R->values[8] = 0.000005;

	kalmangain = matrix_init(kalmangain, 6, 3);
	kalmangain = matrix_alloc(kalmangain);
	kalmangain = matrix_zero(kalmangain);
//	kalmangain = matrix_unit(kalmangain);

	quat_delta = matrix_init(quat_delta, 4, 1);
	quat_delta = matrix_alloc(quat_delta);
	quat_delta = matrix_zero(quat_delta);

	STT = matrix_init(STT, 4, 1);
	STT = matrix_alloc(STT);
	STT = matrix_zero(STT);
	

	out_trq = (matrix_t *)malloc(sizeof(matrix_t));
	out_trq = matrix_init(out_trq, 3, 1);
	out_trq = matrix_alloc(out_trq);
	out_trq = matrix_zero(out_trq);

	in_trq = (matrix_t *)malloc(sizeof(matrix_t));
	in_trq = matrix_init(in_trq, 3, 1);
	in_trq = matrix_alloc(in_trq);
	in_trq = matrix_zero(in_trq);

	quat = (matrix_t *)malloc(sizeof(matrix_t));
	quat = matrix_init(quat, 4, 1);
	quat = matrix_alloc(quat);
	quat->values[0] = 1.0;
	quat->values[1] = 0.0;
	quat->values[2] = 0.0;
	quat->values[3] = 0.0;

	gyro_bias = (matrix_t *)malloc(sizeof(matrix_t));
	gyro_bias = matrix_init(gyro_bias, 3, 1);
	gyro_bias = matrix_alloc(gyro_bias);
	gyro_bias = matrix_zero(gyro_bias);


	I = (matrix_t *)malloc(sizeof(matrix_t));
	I = matrix_init(I, 3, 3);
	I = matrix_alloc(I);

	I->values[0] = 1.0;
	I->values[1] = 0;
	I->values[2] = 0;
	I->values[3] = 1.1;
	I->values[4] = 0;
	I->values[5] = 0;
	I->values[6] = 1.2;
	I->values[7] = 0;
	I->values[8] = 0;




	A = (matrix_t *)malloc(sizeof(matrix_t));
	A->rows = 4;
	A->columns = 1;
	A->values = (double *)malloc(sizeof(double)*A->rows*A->columns);
	matrix_zero(A);
	A->values[0] = 1.0;

	timestep = (matrix_t *)malloc(sizeof(matrix_t));
	timestep->rows = 1;
	timestep->columns = 1;
	timestep->values = (double *)malloc(sizeof(double)*timestep->rows*timestep->columns);
	timestep->values[0] = 0.1;

	T = (matrix_t *)malloc(sizeof(matrix_t));
	T->rows = 3;
	T->columns = 1;
	T->values = (double *)malloc(sizeof(double)*T->rows*T->columns);

	omega = (matrix_t *)malloc(sizeof(matrix_t));
	omega->rows = 3;
	omega->columns = 1;
	omega->values = (double *)malloc(sizeof(double)*omega->rows*omega->columns);

	omega->values[0] = 0.0001;
	omega->values[1] = 0.0002;
	omega->values[2] = 0.003;

	in_trq->values[0] = 0.00003;
	in_trq->values[1] = 0.00003;
	in_trq->values[2] = 0.00003;



		fd = fopen("sim.csv", "w");
	for (int i = 0; i < 100; i++) {
		T->values[0] = 0.001*sin(i);
		T->values[1] = 0.002*cos(i);
		T->values[2] = 0.003*sin(i);

		RK4(DynEq, omega, in_trq, omega);
//		Kinematics(A, omega, A);
//		RK4(Kinematics, A, omega, A);
		ExtendKalmanFilter(Kinematics, quat, omega,quat);
		PIDcontrol(omega, out_trq);
		matrix_add(T, out_trq, in_trq);

		for (int i = 0; i < 4; i++) {
			printf("quat[%d]:%f\n", i, quat->values[i]);
		}
/*		for (int i = 0; i < 3; i++) {
			printf("omega[%d]:%f\n", i, N->values[i]);
		}*/
		/*		for (int i = 0; i < 3; i++) {
		printf("torque[%d]:%f\n",i, out_trq->values[i]);
		//			printf("quat[%d]:%f\n",i, quat->values[i]);

		}*/
				fprintf(fd, "%f,%f,%f,%f,%f,%f,%f,%f,%f\n", STT->values[0], STT->values[1], STT->values[2], STT->values[3], quat->values[0],
					quat->values[1], quat->values[2], quat->values[3], P->values[35]);

	}
		fclose(fd);

	/*	for (int i = 0; i < 9; i++) {
	printf("ret:%f\n", I->values[i]);
	}*/

	free(Ad);
	free(Bd);
	free(Q);
	free(C);
	free(I_kalman);
	free(R);
	free(kalmangain);
	free(quat_delta);
	free(STT);
	free(A); //メモリ解放
	free(P);
	free(quat);
	free(out_trq);
	free(in_trq);
	free(T);
	free(omega);
	free(I);
	free(timestep);
	free(gyro_bias);
	_getch();
	return 0;
}