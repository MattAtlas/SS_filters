//matrix.c
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>

#define DEBUG

// #ifdef DEBUG
// printf("hello world\n");
// #endif

	

// matrix allocation
////////////////////////////////////////////////////////////
Matrix CreateMatrix(int rows, int cols){
	Matrix A;
	if(rows<1 || cols<1){
		printf("error creating matrix, row or col must be >=1");
		return A;
	}
	A.rows = rows;
	A.cols = cols;
	A.mat = (float**)malloc(rows * sizeof(float*));
		for (int i=0; i<rows; i++)
			A.mat[i] = (float*)calloc(cols, sizeof(float));
	return A;
}

Matrix CreateSqrMatrix(int n){
	Matrix A = CreateMatrix(n,n);
	return A;
}

Matrix CreateRowVector(int n){
	Matrix A = CreateMatrix(1,n);
	return A;
}

Matrix CreateColumnVector(int n){
	Matrix A = CreateMatrix(n,1);
	return A;	
} 

int setMatrixEntry(Matrix* A, int row, int col, float val){
	if(A==NULL){
		printf("error, matrix is null pointer\n");
		return -1;
	}
	if(row<0||row>(A->rows-1)){
		printf("Error, row out of bounds\n");
		return -1;
	}
	if(col<0||col>(A->cols-1)){
		printf("Error, col out of bounds\n");
		return -1;
	}
	A->mat[row-1][col-1] = val;
	return 0;
}


int getMatrixEntry(Matrix* A, int row, int col, float* val){
	if(A==NULL){
		printf("error, A is null pointer\n");
		return -1;
	}
	if(row<0||row>(A->rows-1)){
		printf("Error, row out of bounds\n");
		return -1;
	}
	if(col<0||col>(A->cols-1)){
		printf("Error, col out of bounds\n");
		return -1;
	}
	*val = A->mat[row-1][col-1];
	return 0;
}


int getVectorEntry(Matrix* A, int pos, float* val){
	if(A==NULL){
		printf("error, A is null pointer\n");
		return -1;
	}
	if(A->rows!=1 && A->cols!=1){
		printf("error: Matrix A is not a vector\n");
		return -1;
	}
	
	int row, col;
	if(A->cols==1){
		row = pos-1;
		col = 0;
	}
	else {
		col = pos-1;
	}	row = 0;
	
	*val = A->mat[row][col];
	return 0;
}
////////////////////////////////////////////////////////////






SS_filter CreateSSfilter(Matrix A, Matrix B, Matrix C, float dt){
	SS_filter filter;
	
	// error handling ////////////

	
	
	
	filter.states  = A.rows;
	filter.inputs  = B.cols;
	filter.outputs = C.rows;
	
	filter.X0 = CreateColumnVector(filter.states);
	filter.X1 = CreateColumnVector(filter.states);
	filter.Y  = CreateRowVector(filter.outputs);
	
	filter.F = CreateSqrMatrix(filter.states);
	filter.G = CreateMatrix(filter.states,filter.inputs);
	filter.H = CreateMatrix(filter.outputs,filter.states);
	
	filter.F = C2D_A2F(A, dt);
	filter.G = C2D_B2G(A, B, dt);
	filter.H = C;
	
	filter.saturation_en = 0;
	filter.saturation_high = CreateColumnVector(filter.inputs);
	filter.saturation_low  = CreateColumnVector(filter.inputs);
	
	filter.K = CreateMatrix(filter.inputs,filter.states);

	return filter;
}

int marchFilter(SS_filter* sys, Matrix input){
	
	Matrix FX;	// temporary matrix structs
	Matrix Gu;
	
	if (input.rows != sys->G.cols){
		printf("Error: input vector size mismatch");
		return -1;
	}
	
	if (sys->saturation_en == 1){
		for (int i=0;i<input.rows;i++){
			if (input.mat[i][0] < sys->saturation_low.mat[i][0]){
			input.mat[i][0] = sys->saturation_low.mat[i][0];}
			if (input.mat[i][0] > sys->saturation_high.mat[i][0]){
			input.mat[i][0] = sys->saturation_high.mat[i][0];}
		}
	}
	
	multiplyMatrices(&sys->F, &sys->X0, &FX);
	multiplyMatrices(&sys->G, &input, &Gu);
	
	addMatrices(&FX, &Gu, &sys->X1);

	printf("y = %f\n",sys->X0.mat[2][0]*sys->H.mat[0][2]);

	sys->X0 = sys->X1;
	return 0;
}

int multiplyMatrices(Matrix* A, Matrix* B, Matrix* out){
	if (A->cols != B->rows){
		printf("Invalid matrix sizes");
		return -1;
	}
	Matrix result = CreateMatrix(A->rows, B->cols);
	*out = result;	
	
	float sum = 0;
	
	for (int i=0;i<(A->rows);i++){
		for (int j=0;j<(B->cols);j++){	
			for (int k=0;k<(A->cols);k++){
				// do the matrix multiplication
				sum = sum + A->mat[i][k]*B->mat[k][j];
			}
			// save mult sum to new location
			out->mat[i][j] = sum;
			sum = 0; 	// re-initialize sum for next loop
		}
	}
	return 0;
}

int scalarMultiply(Matrix A, float s){
	
	for (int i=0;i<(A.rows);i++){
		for (int j=0;j<(A.cols);j++){	
			A.mat[i][j] = s*A.mat[i][j];
		}
	}
	return 0;
}

int addMatrices(Matrix* A, Matrix* B, Matrix* out){
	if ((A->rows != B->rows)||(A->cols != B->cols)){
		printf("Invalid matrix sizes");
		return -1;
	}
	
	Matrix result = CreateMatrix(A->rows, A->cols);
	*out = result;
	
	for (int i=0;i<(A->rows);i++){
		for (int j=0;j<(A->cols);j++){	
			out->mat[i][j] = A->mat[i][j] + B->mat[i][j];
		}
	}
	return 0;
}

void printMatrix(Matrix A){
	printf("\n");
	for (int i=0;i<A.rows;i++){
		for (int j=0;j<A.cols;j++){
			printf("%f\t",A.mat[i][j]);
		}	
		printf("\n");
	}
}	


// copy information of one matrix to a new memory location 
Matrix duplicateMatrix(Matrix A){
	Matrix temp = CreateMatrix(A.rows,A.cols);
	for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
			temp.mat[i][j] = A.mat[i][j];
		}
	}
	return temp;
}


float getDetMatrix(Matrix A){

	if (A.rows != A.cols){
		printf("Error: Matrix is not square");
		return -1;
	}
	float ratio,det;

	Matrix temp = duplicateMatrix(A);
	
	for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.rows;j++){
            if(j>i){
				ratio = temp.mat[j][i]/temp.mat[i][i];
                for(int k=0;k<A.rows;k++){
                    temp.mat[j][k] = temp.mat[j][k] - ratio * temp.mat[i][k];
                }
            }
        }
    }
	det = 1; //storage for determinant
    for(int i=0;i<A.rows;i++) det = det*temp.mat[i][i];

    return det;  
}



int invertMatrix(Matrix* A, Matrix* out){
	
	float det,coDet;
	
	det = getDetMatrix(*A);
	
	printf("det = %f\n",det);
	if (det == 0){
		printf("Error: Matrix is not invertable");
		return -1;
	}
	
	Matrix result = CreateSqrMatrix(A->rows);
	*out = result;
	
	Matrix cofactors = CreateSqrMatrix(A->rows - 1);
	int i,j,ii,jj,i1,j1;

	for (i=0;i<A->rows;i++){					// current row of A to test
		for (j=0;j<A->rows;j++){				// current col of A to test

			i1 = 0;								// index for cofactor row
			for (ii=0;ii<A->rows;ii++){			// count up thru # of rows of A
				if (ii == i) continue;			// if = to current row of A.. skip
										
				j1 = 0;							// index for cofactor col
				for (jj=0;jj<A->rows;jj++){		// count up thru # of cols of A
					if (jj == j) continue;		// if = to current col of A.. skip
						
					cofactors.mat[i1][j1] = A->mat[ii][jj]; 	// place proper element in new matrix
					j1++;
				}
				i1++;
			}
			coDet = getDetMatrix(cofactors);
			out->mat[j][i] = (pow(-1.0,i+j+2.0) * coDet)/det; 	// saves as transpose
		}
	}
	return 0;
}


// Using CT A matrix and time step, get DT F matrix
Matrix C2D_A2F(Matrix A, float h){
	int m = A.rows;
	Matrix sumold = CreateSqrMatrix(m);
	Matrix sumnew = CreateSqrMatrix(m);
	Matrix result = CreateSqrMatrix(m);
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
Matrix C2D_B2G(Matrix A, Matrix B, float h){
	
	int m = A.rows;
	Matrix sumold = CreateSqrMatrix(m);
	Matrix sumnew = CreateSqrMatrix(m);
	Matrix result = CreateSqrMatrix(m);
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














