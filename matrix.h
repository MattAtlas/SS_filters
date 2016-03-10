#ifndef MATRIX_H_
#define MATRIX_H_


typedef struct SqrMatrix{
	int order;
	float** mat;
}SqrMatrix;

typedef struct Matrix{
	int rows;
	int cols;
	float** mat;
}Matrix;

typedef struct SS_filter{
	float dt;
	int states;
	int inputs;
	int outputs;
	SqrMatrix F;
	Matrix G;
	Matrix H;
	float* X0;
	float* X1;
	float* Y;
}SS_filter;

SqrMatrix CreateSqrMatrix(int n);
Matrix CreateMatrix(int rows, int cols);
float* CreateVector(int n);
float* marchFilter(SS_filter* sys, float input[]);
//int factorial(int num);

SqrMatrix discrete_F(SqrMatrix A, float h);

Matrix discrete_G(SqrMatrix A, Matrix B, float h);
SS_filter CreateSSfilter(SqrMatrix A, Matrix B, Matrix C,float dt);
#endif
