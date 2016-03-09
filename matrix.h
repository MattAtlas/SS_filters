#ifndef MATRIX_H_
#define MATRIX_H_


typedef struct SqrMatrix{
	int order;
	float** mat;
}SqrMatrix;


SqrMatrix CreateSqrMatrix(int n);
float * CreateVector(int n);
int factorial(int num);

SqrMatrix discrete_F(SqrMatrix A, float h);

float * discrete_G(SqrMatrix A, float B[], float h);

#endif
