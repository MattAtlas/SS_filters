//matrix.c
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>


SqrMatrix CreateSqrMatrix(int n){
	SqrMatrix mat;
	mat.order = n;
	mat.mat = (float**)malloc(mat.order * sizeof(float*));
		for (int i=0; i<mat.order; i++)
			mat.mat[i] = (float*)calloc(mat.order, sizeof(float));
	return mat;
}


Matrix CreateMatrix(int rows, int cols){
	Matrix mat;
	mat.rows = rows;
	mat.cols = cols;
	mat.mat = (float**)malloc(mat.rows * sizeof(float*));
		for (int i=0; i<mat.rows; i++)
			mat.mat[i] = (float*)calloc(mat.cols, sizeof(float));
	return mat;
}


SS_filter CreateSSfilter(SqrMatrix A, Matrix B, Matrix C, float dt){
	SS_filter filter;
	
	filter.states  = A.order;
	filter.inputs  = B.cols;
	filter.outputs = C.rows;
	
	filter.X0 = CreateVector(filter.states);
	filter.X1 = CreateVector(filter.states);
	filter.Y  = CreateVector(filter.outputs);
	
	filter.F = CreateSqrMatrix(filter.states);
	filter.G = CreateMatrix(filter.states,filter.inputs);
	filter.H = CreateMatrix(filter.outputs,filter.states);
	
	filter.F = discrete_F(A, dt);
	filter.G = discrete_G(A, B, dt);
	filter.H = C;

	return filter;
}

float* marchFilter(SS_filter* sys, float input[]){
	float* output = CreateVector(&sys->outputs);
	
	float sum1 = 0;
	float sum2 = 0;
	
	for (int i=0;i<sys->states;i++){
		for (int j=0;j<sys->states;j++){
			sum1 = sum1 + sys->F->mat[i][j]*sys->X0[j];
		}
		for (int k=0;k<sys->inputs;k++){
			sum2 = sum2 + sys->G->mat[i][k]*input[k];
		}
		sys->X1[i] = sum1 + sum2;
		printf("X1 = %f\n",sys->X1[i]);
		sum1 = 0;
		sum2 = 0;
	}
	for (int i=0;i<sys->outputs;i++){
		for (int j=0;j<sys->states;j++){
			sum1 = sum1 + sys->H->mat[i][j]*sys->X0[j];
		}
		output[i] = sum1;
		sum1 = 0;
	}
	
	sys->X0 = sys->X1;
	return output;
}





float* CreateVector(int n){
	float* vector;
	vector = malloc(n*sizeof(float));
	return vector;
}


//simple factorial function for input of >= 0
// int factorial(int num){
	// int i;
	// int out = 1;
	// if (num < 0)
		// printf("factorial can't be negative");
	// else {
		// for (i=1;i<num;i++){
			// out = out*i;
		// }
	// }
	// return out;
// }


// Using CT A matrix and time step, get DT F matrix
SqrMatrix discrete_F(SqrMatrix A, float h){
	int m = A.order;
	SqrMatrix sumold = CreateSqrMatrix(m);
	SqrMatrix sumnew = CreateSqrMatrix(m);
	SqrMatrix result = CreateSqrMatrix(m);
	float sum = 0;

	// initialize identity matrix
	for (int i=0;i<m;i++){		
		sumold.mat[i][i] = 1;	// A^0 first element of sum
		result.mat[i][i] = 1;
	}
	// N = order of exponential expansion sum
	for (int N=1;N<5;N++){
		for (int i=0;i<m;i++){
			for (int j=0;j<m;j++){	
				for (int k=0;k<m;k++){
					// do the matrix multiplication
					sum = sum + sumold.mat[i][k]*A.mat[k][j];
				}
				// save mult sum to new location
				sumnew.mat[i][j] = h*sum/N;
				sum = 0; 	// re-initialize sum for next loop
			}
		}	
		for (int i=0;i<m;i++){
			for (int j=0;j<m;j++){
				//matrix exponential expansion sum over N
				result.mat[i][j] = result.mat[i][j] + sumnew.mat[i][j];				
				sumold.mat[i][j] = sumnew.mat[i][j];	
			}
		}
	}
	return result;
}

	

// Using CT A and B matrices and time step, get DT G matrix
Matrix discrete_G(SqrMatrix A, Matrix B, float h){
	
	int m = A.order;
	SqrMatrix sumold = CreateSqrMatrix(m);
	SqrMatrix sumnew = CreateSqrMatrix(m);
	SqrMatrix result = CreateSqrMatrix(m);
	float sum = 0;
	Matrix G = CreateMatrix(B.rows,B.cols);

	// initialize identity matrix
	for (int i=0;i<m;i++){		
		sumold.mat[i][i] = 1;	// A^0 first element of sum
		result.mat[i][i] = 1;
	}
	// N = order of exponential expansion sum
	for (int N=1;N<5;N++){				
		for (int i=0;i<m;i++){
			for (int j=0;j<m;j++){	
				for (int k=0;k<m;k++){
					// do the matrix multiplication
					sum = sum + sumold.mat[i][k]*A.mat[k][j];
				}
				// save mult sum to new location
				sumnew.mat[i][j] = h*sum/(N+1);// factorial starts 1 higher than in F
				sum = 0; 	// re-initialize sum for next loop
			}
		}	
		for (int i=0;i<m;i++){
			for (int j=0;j<m;j++){
				//matrix exponential expansion sum over N
				result.mat[i][j] = result.mat[i][j] + sumnew.mat[i][j];				
				sumold.mat[i][j] = sumnew.mat[i][j];	
			}
		}
	}
	for (int i=0;i<m;i++){
		for (int j=0;j<B.cols;j++){
			for (int k=0;k<B.rows;k++){
			sum = sum + result.mat[i][k]*B.mat[k][j];
			}
			G.mat[i][j] = h*sum;
			sum = 0;
		}
	}
	return G;
}
















