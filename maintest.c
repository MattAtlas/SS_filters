#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
	
	float buffer[250] = {0};
	

	FILE *fp = fopen("file.txt", "w+");
	if (fp == NULL){
		printf("Error opening file!\n");
		exit(-1);
	}

	for (int i=0;i<250;i++){
		marchFilter(&sys,u);
		buffer[i] = sys.X0.mat[2][0]*5;
		
	}
	for (int i=0;i<250;i++){
	fprintf(fp,"%f\n",buffer[i]);
	}
	
	fclose(fp);
	Matrix Test = CreateSqrMatrix(3);
	
	invertMatrix(&sys.F,&Test);
	
	printMatrix(Test);
	
	
	
	return 0;
}
