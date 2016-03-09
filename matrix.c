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


float * CreateVector(int n){
	float * vector;
	vector = malloc(n*sizeof(float));
	return vector;
}


//simple factorial function for input of >= 0
int factorial(int num){
	int i;
	int out = 1;
	if (num < 0)
		printf("factorial can't be negative");
	else {
		for (i=1;i<num;i++){
			out = out*i;
		}
	}
	return out;
}


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
float * discrete_G(SqrMatrix A, float B[], float h){
	int m = A.order;
	SqrMatrix sumold = CreateSqrMatrix(m);
	SqrMatrix sumnew = CreateSqrMatrix(m);
	SqrMatrix result = CreateSqrMatrix(m);
	float sum = 0;
	float * G = malloc(m*sizeof *G);
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
		for (int j=0;j<m;j++){
			sum = sum + result.mat[i][j]*B[j];
		}
		G[i] = h*sum;
		sum = 0;
	}
	return G;
}
