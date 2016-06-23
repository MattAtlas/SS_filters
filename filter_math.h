#ifndef MATRIX_H_
#define MATRIX_H_

// add this into Robotics_Cape.h

typedef struct CT_SS_filter_t{

	int states;
	int inputs;
	int outputs;
	int saturation_en;
	int saturation_flag;
	vector_t saturation_high;
	vector_t saturation_low;
	matrix_t A;
	matrix_t B;
	matrix_t C; 
	vector_t Xnew;
	vector_t Xold;
	matrix_t Y;
	matrix_t K;
	matrix_t L;
}CT_SS_filter_t;



typedef struct DT_SS_filter_t{
	float dt;
	int states;
	int inputs;
	int outputs;
	int saturation_en;
	int saturation_flag;
	vector_t saturation_high;
	vector_t saturation_low;
	matrix_t F;
	matrix_t G;
	matrix_t H; 
	vector_t Xnew;
	vector_t Xold;
	matrix_t Y;
	matrix_t K;
	matrix_t L;
}DT_SS_filter_t;


int tf2ss(matrix_t num, matrix_t den, CT_SS_filter_t* CT_sys);

matrix_t C2D_A2F(matrix_t A, float h);
matrix_t C2D_B2G(matrix_t A, matrix_t B, float h);

CT_SS_filter_t createCTSSfilter(int states,int inputs,int outputs);

DT_SS_filter_t createDTSSfilter(CT_SS_filter_t* CT_sys,float dt);

int marchFilter(DT_SS_filter_t* sys, matrix_t input);

vector_t convolve(vector_t v1, vector_t v2);
//int saturate(SS_filter* sys, matrix_t* input);

#endif
