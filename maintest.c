#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"


int main()
{
	int i,j;
	float h = 0.01;
	float in = 1;
	float* u;
	u = &in;
	//float* y;
	SqrMatrix A = CreateSqrMatrix(3);
	Matrix B = CreateMatrix(3,1);
	Matrix C = CreateMatrix(1,3);
	C.mat[0][1] = 1;
	C.mat[0][2] = 5;
	
	B.mat[0][0] = 1;
	
	A.mat[0][0] = -5;
	A.mat[0][1] = -5;
	A.mat[0][2] = -5;
	
	A.mat[1][0] = 1;
	A.mat[2][1] = 1;
	
	SS_filter sys;
	sys = CreateSSfilter(A,B,C,h);
	
	
	sys.X0[0] = 0;
	sys.X0[1] = 0;
	sys.X0[2] = 0;
	
	float* y;
	
	for (int i=0;i<100;i++){
		y = marchFilter(sys,u);
	}
	printf("\n");
	for (int i=0;i<sys.states;i++){
		for (int j=0;j<sys.states;j++){
			printf("%0.10f\t",sys.F.mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for (int i=0;i<sys.states;i++){
		for (int j=0;j<sys.inputs;j++){
			printf("%0.10f\t",sys.G.mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("y = %f\n",*y);
	printf("states = %d\ninputs = %d\noutputs = %d\n",sys.states,sys.inputs,sys.outputs);
	return 0;
}
