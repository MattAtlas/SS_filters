#ifndef MATRIX_H_
#define MATRIX_H_


typedef struct Matrix{
	int rows;
	int cols;
	float** mat;
}Matrix;

typedef struct SS_filter{
	float dt;
	int states;
	int inputs;
	int outputs;
	int saturation_en;
	Matrix saturation_high;
	Matrix saturation_low;
	Matrix F;
	Matrix G;
	Matrix H; 
	Matrix X0;
	Matrix X1;
	Matrix Y;
	Matrix K;
}SS_filter;

// all arrays, 1D & 2D are stored in a matrix struct
// like matlab. Use handy functions
Matrix CreateMatrix(int rows, int cols);
int setMatrixEntry(Matrix* A, int row, int col, float val);
int getMatrixEntry(Matrix* A, int row, int col, float* val);
int getVectorEntry(Matrix* A, int pos, float* val);
Matrix CreateSqrMatrix(int n);
Matrix CreateRowVector(int n);
Matrix CreateColumnVector(int n);
Matrix duplicateMatrix(Matrix A);
// Matrix operations
int multiplyMatrices(Matrix* A, Matrix* B, Matrix* out);
int addMatrices(Matrix* A, Matrix* B, Matrix* out);
int invertMatrix(Matrix* A, Matrix* out);
int scalarMultiply(Matrix A, float s);
float getDetMatrix(Matrix A);

void printMatrix(Matrix A);

// convertions from continuous to discrete time
Matrix C2D_A2F(Matrix A, float h);
Matrix C2D_B2G(Matrix A, Matrix B, float h);

// construct and allocate memory for a discrete time filter
SS_filter CreateSSfilter(Matrix A, Matrix B, Matrix C,float dt);

// March discrete filter forward using a new input
// returns 0 or -1 for success or failure
// user can grab new output from sys
int marchFilter(SS_filter* sys, Matrix input);

#endif
