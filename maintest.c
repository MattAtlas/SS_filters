#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"


int main()
{
	float h = 0.1;

	Matrix num = CreateRowVector(2);
	Matrix den = CreateRowVector(4);
	num.mat[0][0] = 1;
	num.mat[0][1] = 5;
	den.mat[0][0] = 1;
	den.mat[0][1] = 5;
	den.mat[0][2] = 5;
	den.mat[0][3] = 5;
	
	
	CT_SS_filter CT_sys = CreateCTSSfilter(3,1,1);
	tf2ss(num, den, &CT_sys);

	Matrix u = CreateColumnVector(1);

	u.mat[0][0] = 1;
	
	DT_SS_filter DT_sys = CreateDTSSfilter(&CT_sys,h);

	
	float buffer[250] = {0};
	

	FILE *fp = fopen("file.txt", "w+");
	if (fp == NULL){
		printf("Error opening file!\n");
		exit(-1);
	}

	for (int i=0;i<250;i++){
		marchFilter(&DT_sys,u);
		buffer[i] = DT_sys.X0.mat[2][0]*5;
		
	}
	for (int i=0;i<250;i++){
	fprintf(fp,"%f\n",buffer[i]);
	}
	
	fclose(fp);

	Matrix Test = CreateSqrMatrix(3);
	
	invertMatrix(&DT_sys.F,&Test);
	
	printMatrix(CT_sys.A);
	printMatrix(CT_sys.B);
	printMatrix(CT_sys.C);
	printMatrix(DT_sys.F);
	printMatrix(DT_sys.G);
	printMatrix(DT_sys.H);

	
	
	
	
	return 0;
}
