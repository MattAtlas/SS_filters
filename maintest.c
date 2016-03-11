#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"


int main()
{
	float h = 0.1;

	Matrix A = CreateSqrMatrix(3);
	Matrix B = CreateColumnVector(3);
	Matrix C = CreateRowVector(3);
	Matrix u = CreateColumnVector(1);

	A.mat[0][0] = -5;
	A.mat[0][1] = -5;
	A.mat[0][2] = -5;
	A.mat[1][0] = 1;
	A.mat[2][1] = 1;
	
	B.mat[0][0] = 1;

	C.mat[0][1] = 1;
	C.mat[0][2] = 5;	
	
	u.mat[0][0] = 1;
	
	SS_filter sys = CreateSSfilter(A,B,C,h);
	
	sys.X0.mat[0][0] = 0;
	sys.X0.mat[1][0] = 0;
	sys.X0.mat[2][0] = 0;	
	
	for (int i=0;i<250;i++){
		marchFilter(&sys,u);
	}
	
	printMatrix(&sys.F);
	
	printMatrix(&sys.G);
	
	getDetMatrix(&A);
	
	return 0;
}
