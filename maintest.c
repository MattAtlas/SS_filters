#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"


int main()
{
	int i,j;
	float h = 0.01;
	float B[3] = {0,0,1};

	float * G = CreateVector(3);
	SqrMatrix A = CreateSqrMatrix(3);
	
	A.mat[0][0] = 3.14;
	A.mat[0][1] = 5.54;
	
	A.mat[1][0] = 5.;
	A.mat[1][2] = 10.;
	
	A.mat[2][0] = -12;
	A.mat[2][1] = -2;

	SqrMatrix F = CreateSqrMatrix(3);
	F = discrete_F(A,h);
	G = discrete_G(A,B,h);
	printf("\n");
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			printf("%0.10f\t",F.mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for (int i=0;i<3;i++){
		printf("%0.10f\n",G[i]);
	}
	printf("\n");
	return 0;
}
