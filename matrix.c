//matrix.c
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>

#define DEBUG



	

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

// Sets up all matrices for system
CT_SS_filter CreateCTSSfilter(int states, int inputs, int outputs){
	CT_SS_filter filter;
	
	// error handling ////////////
	filter.states  = states;
	filter.inputs  = inputs;
	filter.outputs = outputs;
	
	filter.A = CreateSqrMatrix(filter.states);
	filter.B = CreateMatrix(filter.states,filter.inputs);
	filter.C = CreateMatrix(filter.outputs,filter.states);
	
	filter.X0 = CreateColumnVector(filter.states);
	filter.X1 = CreateColumnVector(filter.states);
	filter.Y  = CreateRowVector(filter.outputs);

	
	filter.saturation_en = 0;
	filter.saturation_flag = 0;
	filter.saturation_high = CreateColumnVector(filter.inputs);
	filter.saturation_low  = CreateColumnVector(filter.inputs);
	
	filter.K = CreateMatrix(filter.inputs,filter.states);
	filter.L = CreateMatrix(filter.states,filter.outputs);
	
	return filter;
}



// Sets up all matrices for system
DT_SS_filter CreateDTSSfilter(CT_SS_filter* CT_sys, float dt){
	DT_SS_filter filter;
	
	// error handling ////////////
	filter.dt = dt;

	filter.states  = CT_sys->A.rows;
	filter.inputs  = CT_sys->B.cols;
	filter.outputs = CT_sys->C.rows;
	
	filter.X0 = CreateColumnVector(filter.states);
	filter.X1 = CreateColumnVector(filter.states);
	filter.Y  = CreateRowVector(filter.outputs);
	
	filter.F = CreateSqrMatrix(filter.states);
	filter.G = CreateMatrix(filter.states,filter.inputs);
	filter.H = CreateMatrix(filter.outputs,filter.states);

	// DT model of system	
	filter.F = C2D_A2F(CT_sys->A, filter.dt);
	filter.G = C2D_B2G(CT_sys->A, CT_sys->B, filter.dt);
	filter.H = CT_sys->C;
	
	filter.saturation_en = 0;
	filter.saturation_flag = 0;
	filter.saturation_high = CreateColumnVector(filter.inputs);
	filter.saturation_low  = CreateColumnVector(filter.inputs);
	
	filter.K = CreateMatrix(filter.inputs,filter.states);
	filter.L = CreateMatrix(filter.states,filter.outputs);
	
	return filter;
}

int tf2ss(Matrix b, Matrix a, CT_SS_filter* CT_sys){
	
	// make sure it's proper
	if (b.cols >= a.cols){
		printf("Error: Improper transfer function\n");
		return -1;
	}
	
	Matrix b_temp = CreateRowVector(a.cols);
	// allocate memory for zero coeffs
	if (b.cols < a.cols){
		
		for(int i=1;i<=b.cols;i++){
			b_temp.mat[0][a.cols - i] = b.mat[0][b.cols - i];
		}
		
		// Matrix b = CreateRowVector(a.cols);
		// for(int i=0;i<a.cols;i++){
			// b.mat[0][i] = temp.mat[0][i];
		// }
	}
	
	
	if (a.mat[0][0] != 1){
		float coeff = a.mat[0][0];
		for (int i=0;i<=a.cols;i++){	
			a.mat[0][i] = a.mat[0][i] / coeff;	
			b_temp.mat[0][i] = b_temp.mat[0][i] / coeff;
		}
	}

	// Fill in A matrix
	for (int i=0;i<a.cols-2;i++){
		CT_sys->A.mat[i+1][i] = 1;			// fill in lower identity
	}
	for (int i=0;i<a.cols-1;i++){
		CT_sys->A.mat[0][i] = -a.mat[0][i+1];	// fill in top row
	}
	printMatrix(CT_sys->A);
	
	CT_sys->B.mat[0][0] = 1;
	
	for (int i=0;i<a.cols-1;i++){
		CT_sys->C.mat[0][i] = b_temp.mat[0][i+1];
	}
	
	return 0;
}






/*
int saturate(DT_SS_filter* sys, Matrix* input){
	
	for (int i=0;i<sys->inputs;i++){
		if (&input[i] > sys->saturation_high[i]){
			&input[i] = sys->saturation_high[i];
			
			sys->saturation_flag = 1;
		}
		else if(&input[i] < sys->saturation_low[i]){
				&input[i] = sys->saturation_low[i];
			
			sys->saturation_flag = 1;
		}
		else{
			sys->saturation_flag = 0;
		}
	}
	
	return 0;
}
*/

int marchFilter(DT_SS_filter* DT_sys, Matrix input){
	
	Matrix FX;	// temporary matrix structs
	Matrix Gu;
	DT_sys->X0 = DT_sys->X1;
	
	if (input.rows != DT_sys->G.cols){
		printf("Error: input vector size mismatch");
		return -1;
	}
	
	if (DT_sys->saturation_en == 1){
		for (int i=0;i<input.rows;i++){
			if (input.mat[i][0] < DT_sys->saturation_low.mat[i][0]){
			input.mat[i][0] = DT_sys->saturation_low.mat[i][0];}
			if (input.mat[i][0] > DT_sys->saturation_high.mat[i][0]){
			input.mat[i][0] = DT_sys->saturation_high.mat[i][0];}
		}
	}
	
	multiplyMatrices(&DT_sys->F, &DT_sys->X1, &FX);
	multiplyMatrices(&DT_sys->G, &input, &Gu);
	
	addMatrices(&FX, &Gu, &DT_sys->X0);

	#ifdef DEBUG
	printf("y = %f\n",DT_sys->X1.mat[2][0]*DT_sys->H.mat[0][2]);
	#endif
	


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












