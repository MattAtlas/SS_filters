/*******************************************************************************
* filter_math.c
*
* Matt Atlas 2016
*******************************************************************************/

#include "../robotics_cape.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG


/*******************************************************************************
* 
*
* Sets up all matrices for system
*******************************************************************************/
CT_SS_filter_t createCTSSfilter(int states, int inputs, int outputs){
	CT_SS_filter_t filter;
	
	// error handling
	filter.states  = states;
	filter.inputs  = inputs;
	filter.outputs = outputs;
	
	filter.A = createSquareMatrix(filter.states);
	filter.B = createMatrix(filter.states,filter.inputs);
	filter.C = createMatrix(filter.outputs,filter.states);
	
	filter.Xold = createVector(filter.states);
	filter.Xnew = createVector(filter.states);
	filter.Y    = createVector(filter.outputs);

	
	filter.saturation_en   = 0;
	filter.saturation_flag = 0;
	filter.saturation_high = createVector(filter.inputs);
	filter.saturation_low  = createVector(filter.inputs);
	
	filter.K = createMatrix(filter.inputs,filter.states);
	filter.L = createMatrix(filter.states,filter.outputs);
	
	return filter;
}


/*******************************************************************************
* 
*
* Sets up all matrices for system
*******************************************************************************/
DT_SS_filter_t createDTSSfilter(CT_SS_filter_t* CT_sys, float dt){
	DT_SS_filter_t filter;
	
	// error handling
	filter.dt = dt;

	filter.states  = CT_sys->A.rows;
	filter.inputs  = CT_sys->B.cols;
	filter.outputs = CT_sys->C.rows;
	
	filter.Xold = createVector(filter.states);
	filter.Xnew = createVector(filter.states);
	filter.Y    = createVector(filter.outputs);
	
	filter.F = createSqrareMatrix(filter.states);
	filter.G = createMatrix(filter.states,filter.inputs);
	filter.H = createMatrix(filter.outputs,filter.states);

	// DT model of system	
	filter.F = C2D_A2F(CT_sys->A, filter.dt);
	filter.G = C2D_B2G(CT_sys->A, CT_sys->B, filter.dt);
	filter.H = CT_sys->C;
	
	filter.saturation_en   = 0;
	filter.saturation_flag = 0;
	filter.saturation_high = createVector(filter.inputs);
	filter.saturation_low  = createVector(filter.inputs);
	
	filter.K = createMatrix(filter.inputs,filter.states);
	filter.L = createMatrix(filter.states,filter.outputs);
	
	return filter;
}



/*******************************************************************************
* 
*
* 
*******************************************************************************/
int tf2ss(vector_t b, vector_t a, CT_SS_filter_t* CT_sys){
	int i;
	// make sure it's proper
	if (b.len >= a.len){
		printf("Error: Improper transfer function\n");
		return -1;
	}
	
	// make new b vector the same length as a
	vector_t b_temp = createVector(a.len);
	// allocate memory for zero coeffs

	// put the values of b in new b starting at the end	
	for(i=1;i<=b.len;i++){
		b_temp.data[a.len - i] = b.data[b.len - i];
	}

	if(a.data[0] != 1){
		float coeff = a.data[0];
		for(i=0;i<=a.len;i++){	
			a.data[i] = a.data[i] / coeff;	
			b_temp.data[i] = b_temp.data[i] / coeff;
		}
	}

	// Fill in A matrix
	for(i=0;i<a.len-2;i++){
		CT_sys->A.data[i+1][i] = 1;			// fill in lower identity
	}
	for(i=0;i<a.len-1;i++){
		CT_sys->A.data[0][i] = -a.data[0][i+1];	// fill in top row
	}
	printMatrix(CT_sys->A);
	
	CT_sys->B.data[0][0] = 1;
	
	for(i=0;i<a.len-1;i++){
		CT_sys->C.data[0][i] = b_temp.data[i+1];
	}
	
	return 0;
}





/*******************************************************************************
* 
*
* 
*******************************************************************************/
// int saturate(DT_SS_filter_t* sys, matrix_t* input){
	
	// for (int i=0;i<sys->inputs;i++){
		// if (&input[i] > sys->saturation_high[i]){
			// &input[i] = sys->saturation_high[i];
			
			// sys->saturation_flag = 1;
		// }
		// else if(&input[i] < sys->saturation_low[i]){
				// &input[i] = sys->saturation_low[i];
			
			// sys->saturation_flag = 1;
		// }
		// else{
			// sys->saturation_flag = 0;
		// }
	// }
	
	// return 0;
// }

/*******************************************************************************
* 
*
* 
*******************************************************************************/
int marchFilter(DT_SS_filter_t* DT_sys, vector_t input){
	
	matrix_t FX;	// temporary matrix structs
	matrix_t Gu;
	int i;
	
	DT_sys->Xold = DT_sys->Xnew;
	
	if(input.len != DT_sys->G.cols){
		printf("Error: input vector size mismatch");
		return -1;
	}
	
	if(DT_sys->saturation_en == 1){
		for(i=0;i<input.len;i++){
			if(input.data[i] < DT_sys->saturation_low.data[i]){
				input.data[i] = DT_sys->saturation_low.data[i];
				DT_sys->saturation_flag = 1;
			}
			if(input.data[i] > DT_sys->saturation_high.data[i]){
				input.data[i] = DT_sys->saturation_high.data[i];
				DT_sys->saturation_flag = 1;
			}
		}
	}
	
	FX = matrixMultiply(DT_sys->F, DT_sys->Xold);
	Gu = matrixMultiply(DT_sys->G, input);
	
	DT_sys->Xnew = matrixAdd(FX, Gu);

	#ifdef DEBUG
	// output y is estimate of sensor data
	printf("y = %f\n",DT_sys->Xold.data[2][0]*DT_sys->H.data[0][2]);
	#endif
	
	return 0;
}



/*******************************************************************************
* matrix_t C2D_A2F(matrix_t A, float h)
*
* Using CT A matrix and time step, get DT F matrix
*******************************************************************************/
matrix_t C2D_A2F(matrix_t A, float h){
	int m = A.rows;
	int i,j,k,N;
	matrix_t sumold = createSqrareMatrix(m);
	matrix_t sumnew = createSqrareMatrix(m);
	matrix_t result = createSqrareMatrix(m);
	float sum;

	// initialize identity matrix
	for(i=0;i<m;i++){		
		sumold.data[i][i] = 1;	// A^0 first element of sum
		result.data[i][i] = 1;
	}
	// N = order of exponential expansion sum
	for(N=1;N<5;N++){
		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				sum = 0; 	// initialize sum for next loop
				for(k=0;k<m;k++){
					// do the matrix multiplication
					sum += sumold.data[i][k]*A.data[k][j];
				}
				// save mult sum to new location
				sumnew.data[i][j] = h*sum/N;
			}
		}	
		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				//matrix exponential expansion sum over N
				result.data[i][j] += sumnew.data[i][j];				
				sumold.data[i][j] = sumnew.data[i][j];	
			}
		}
	}
	return result;
}

	
/*******************************************************************************
* 
*
* Using CT A and B matrices and time step, get DT G matrix
*******************************************************************************/
matrix_t C2D_B2G(matrix_t A, matrix_t B, float h){
	int i,j,k,N;
	int m = A.rows;
	matrix_t sumold = CreateSqrMatrix(m);
	matrix_t sumnew = CreateSqrMatrix(m);
	matrix_t result = CreateSqrMatrix(m);
	float sum;
	matrix_t G = CreateMatrix(B.rows,B.cols);

	// initialize identity matrix
	for(i=0;i<m;i++){		
		sumold.data[i][i] = 1;	// A^0 first element of sum
		result.data[i][i] = 1;
	}
	// N = order of exponential expansion sum
	for(N=1;N<5;N++){				
		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				sum = 0; 			// initialize sum for next loop
				for(k=0;k<m;k++){
					// do the matrix multiplication
					sum += sumold.data[i][k]*A.data[k][j];
				}
				// save mult sum to new location
				sumnew.data[i][j] = h*sum/(N+1);// factorial starts 1 higher than in F
			}
		}	
		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				//matrix exponential expansion sum over N
				result.data[i][j] += sumnew.data[i][j];				
				sumold.data[i][j] = sumnew.data[i][j];	
			}
		}
	}
	for(i=0;i<m;i++){
		sum = 0;
		for(j=0;j<B.cols;j++){
			for(k=0;k<B.rows;k++){
			sum += result.data[i][k]*B.data[k][j];
			}
		G.data[i][j] = h*sum;
		}
	}
	return G;
}



vector_t convolve(vector_t v1, vector_t v2){
	
	int m,n,i,j,k;
	vector_t out;
	if(!v1.initialized || !v2.initialized){
		printf("ERROR: vector not initialized yet\n");
		return -1;
	}
	m = v1.len;
	n = v2.len;
	k = m+n-1;
	out = createVector(k);
	
	for(i=0;i<m;i++){
		for(j=0;j<n;i++){
			out(i+j) += u(i) * v(j);
		}
	}
	return out;	
}
