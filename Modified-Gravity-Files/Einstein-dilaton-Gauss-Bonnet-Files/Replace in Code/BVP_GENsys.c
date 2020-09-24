#include "BVP_header.h"

int BVP_GENsys(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[])
{
	int i,j,k,l, status;
	const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;
	const int nnzJac = (*params).nnzparam;
	const double r_H = (*params).r_Hparam;
	const double M = (*params).Mparam;
	const double alpha = (*params).alphaparam;
	const double beta = (*params).betaparam;
	const size_t slen = r+4;
	int pn;
	const double sgn = -1.0;
	int a, b; //a for x, b for y(theta)
	int I, onexside, twoxside;
	int J, oneyside, twoyside;
	double xk, yk;
	const int dimx = 1;
	const int dimy = n;
	int nnztemp;
	
	int *sx;
	sx = (int*)calloc(slen, sizeof(int));
	int *sy;
	sy = (int*)calloc(slen, sizeof(int));
	
	//----- Finite difference stencil -----
	double *a0d0, *a0d1, *a0d2;
	a0d0 = (double*)calloc(slen, sizeof(double));
	a0d1 = (double*)calloc(slen, sizeof(double));
	a0d2 = (double*)calloc(slen, sizeof(double));
	double *b0d0, *b0d1, *b0d2;
	b0d0 = (double*)calloc(slen, sizeof(double));
	b0d1 = (double*)calloc(slen, sizeof(double));
	b0d2 = (double*)calloc(slen, sizeof(double));
	
	double *a1d0, *a1d1, *a1d2;
	a1d0 = (double*)calloc(slen, sizeof(double));
	a1d1 = (double*)calloc(slen, sizeof(double));
	a1d2 = (double*)calloc(slen, sizeof(double));
	double *b1d0, *b1d1, *b1d2;
	b1d0 = (double*)calloc(slen, sizeof(double));
	b1d1 = (double*)calloc(slen, sizeof(double));
	b1d2 = (double*)calloc(slen, sizeof(double));
	
	double *a2d0, *a2d1, *a2d2;
	a2d0 = (double*)calloc(slen, sizeof(double));
	a2d1 = (double*)calloc(slen, sizeof(double));
	a2d2 = (double*)calloc(slen, sizeof(double));
	double *b2d0, *b2d1, *b2d2;
	b2d0 = (double*)calloc(slen, sizeof(double));
	b2d1 = (double*)calloc(slen, sizeof(double));
	b2d2 = (double*)calloc(slen, sizeof(double));
	
	double *a3d0, *a3d1, *a3d2;
	a3d0 = (double*)calloc(slen, sizeof(double));
	a3d1 = (double*)calloc(slen, sizeof(double));
	a3d2 = (double*)calloc(slen, sizeof(double));
	double *b3d0, *b3d1, *b3d2;
	b3d0 = (double*)calloc(slen, sizeof(double));
	b3d1 = (double*)calloc(slen, sizeof(double));
	b3d2 = (double*)calloc(slen, sizeof(double));
	
	double *a4d0, *a4d1, *a4d2;
	a4d0 = (double*)calloc(slen, sizeof(double));
	a4d1 = (double*)calloc(slen, sizeof(double));
	a4d2 = (double*)calloc(slen, sizeof(double));
	double *b4d0, *b4d1, *b4d2;
	b4d0 = (double*)calloc(slen, sizeof(double));
	b4d1 = (double*)calloc(slen, sizeof(double));
	b4d2 = (double*)calloc(slen, sizeof(double));
	
	double *a5d0, *a5d1, *a5d2;
	a5d0 = (double*)calloc(slen, sizeof(double));
	a5d1 = (double*)calloc(slen, sizeof(double));
	a5d2 = (double*)calloc(slen, sizeof(double));
	double *b5d0, *b5d1, *b5d2;
	b5d0 = (double*)calloc(slen, sizeof(double));
	b5d1 = (double*)calloc(slen, sizeof(double));
	b5d2 = (double*)calloc(slen, sizeof(double));
	
	double *a6d0, *a6d1, *a6d2;
	a6d0 = (double*)calloc(slen, sizeof(double));
	a6d1 = (double*)calloc(slen, sizeof(double));
	a6d2 = (double*)calloc(slen, sizeof(double));
	double *b6d0, *b6d1, *b6d2;
	b6d0 = (double*)calloc(slen, sizeof(double));
	b6d1 = (double*)calloc(slen, sizeof(double));
	b6d2 = (double*)calloc(slen, sizeof(double));
	
	double *a7d0, *a7d1, *a7d2;
	a7d0 = (double*)calloc(slen, sizeof(double));
	a7d1 = (double*)calloc(slen, sizeof(double));
	a7d2 = (double*)calloc(slen, sizeof(double));
	double *b7d0, *b7d1, *b7d2;
	b7d0 = (double*)calloc(slen, sizeof(double));
	b7d1 = (double*)calloc(slen, sizeof(double));
	b7d2 = (double*)calloc(slen, sizeof(double));
	
	double *a8d0, *a8d1, *a8d2;
	a8d0 = (double*)calloc(slen, sizeof(double));
	a8d1 = (double*)calloc(slen, sizeof(double));
	a8d2 = (double*)calloc(slen, sizeof(double));
	double *b8d0, *b8d1, *b8d2;
	b8d0 = (double*)calloc(slen, sizeof(double));
	b8d1 = (double*)calloc(slen, sizeof(double));
	b8d2 = (double*)calloc(slen, sizeof(double));
	
	double *a9d0, *a9d1, *a9d2;
	a9d0 = (double*)calloc(slen, sizeof(double));
	a9d1 = (double*)calloc(slen, sizeof(double));
	a9d2 = (double*)calloc(slen, sizeof(double));
	double *b9d0, *b9d1, *b9d2;
	b9d0 = (double*)calloc(slen, sizeof(double));
	b9d1 = (double*)calloc(slen, sizeof(double));
	b9d2 = (double*)calloc(slen, sizeof(double));
	
	//----- Field derivatives initialization -----
	double u0d00k,  u0d10k,  u0d20k,  u0d01k,  u0d02k;
	double u0d00dk, u0d10dk, u0d20dk, u0d01dk, u0d02dk;
	
	double u1d00k,  u1d10k,  u1d20k,  u1d01k,  u1d02k;
	double u1d00dk, u1d10dk, u1d20dk, u1d01dk, u1d02dk;
	
	double u2d00k,  u2d10k,  u2d20k,  u2d01k,  u2d02k;
	double u2d00dk, u2d10dk, u2d20dk, u2d01dk, u2d02dk;
	
	double u3d00k,  u3d10k,  u3d20k,  u3d01k,  u3d02k;
	double u3d00dk, u3d10dk, u3d20dk, u3d01dk, u3d02dk;
	
	double u4d00k,  u4d10k,  u4d20k,  u4d01k,  u4d02k;
	double u4d00dk, u4d10dk, u4d20dk, u4d01dk, u4d02dk;
	
	double u5d00k,  u5d10k,  u5d20k,  u5d01k,  u5d02k;
	double u5d00dk, u5d10dk, u5d20dk, u5d01dk, u5d02dk;
	
	double u6d00k,  u6d10k,  u6d20k,  u6d01k,  u6d02k;
	double u6d00dk, u6d10dk, u6d20dk, u6d01dk, u6d02dk;
	
	double u7d00k,  u7d10k,  u7d20k,  u7d01k,  u7d02k;
	double u7d00dk, u7d10dk, u7d20dk, u7d01dk, u7d02dk;
	
	double u8d00k,  u8d10k,  u8d20k,  u8d01k,  u8d02k;
	double u8d00dk, u8d10dk, u8d20dk, u8d01dk, u8d02dk;
	
	double u9d00k,  u9d10k,  u9d20k,  u9d01k,  u9d02k;
	double u9d00dk, u9d10dk, u9d20dk, u9d01dk, u9d02dk;
	
	//----- Jacobian partials -----
	double FE0;
	double dFE0d0d00, dFE0d0d10, dFE0d0d01, dFE0d0d20, dFE0d0d02;
	double dFE0d1d00, dFE0d1d10, dFE0d1d01, dFE0d1d20, dFE0d1d02;
	double dFE0d2d00, dFE0d2d10, dFE0d2d01, dFE0d2d20, dFE0d2d02;
	double dFE0d3d00, dFE0d3d10, dFE0d3d01, dFE0d3d20, dFE0d3d02;
	double dFE0d4d00, dFE0d4d10, dFE0d4d01, dFE0d4d20, dFE0d4d02;
	double dFE0d5d00, dFE0d5d10, dFE0d5d01, dFE0d5d20, dFE0d5d02;
	double dFE0d6d00, dFE0d6d10, dFE0d6d01, dFE0d6d20, dFE0d6d02;
	double dFE0d7d00, dFE0d7d10, dFE0d7d01, dFE0d7d20, dFE0d7d02;
	double dFE0d8d00, dFE0d8d10, dFE0d8d01, dFE0d8d20, dFE0d8d02;
	double dFE0d9d00, dFE0d9d10, dFE0d9d01, dFE0d9d20, dFE0d9d02;
	
	double FE1;
	double dFE1d0d00, dFE1d0d10, dFE1d0d01, dFE1d0d20, dFE1d0d02;
	double dFE1d1d00, dFE1d1d10, dFE1d1d01, dFE1d1d20, dFE1d1d02;
	double dFE1d2d00, dFE1d2d10, dFE1d2d01, dFE1d2d20, dFE1d2d02;
	double dFE1d3d00, dFE1d3d10, dFE1d3d01, dFE1d3d20, dFE1d3d02;
	double dFE1d4d00, dFE1d4d10, dFE1d4d01, dFE1d4d20, dFE1d4d02;
	double dFE1d5d00, dFE1d5d10, dFE1d5d01, dFE1d5d20, dFE1d5d02;
	double dFE1d6d00, dFE1d6d10, dFE1d6d01, dFE1d6d20, dFE1d6d02;
	double dFE1d7d00, dFE1d7d10, dFE1d7d01, dFE1d7d20, dFE1d7d02;
	double dFE1d8d00, dFE1d8d10, dFE1d8d01, dFE1d8d20, dFE1d8d02;
	double dFE1d9d00, dFE1d9d10, dFE1d9d01, dFE1d9d20, dFE1d9d02;
	
	double FE2;
	double dFE2d0d00, dFE2d0d10, dFE2d0d01, dFE2d0d20, dFE2d0d02;
	double dFE2d1d00, dFE2d1d10, dFE2d1d01, dFE2d1d20, dFE2d1d02;
	double dFE2d2d00, dFE2d2d10, dFE2d2d01, dFE2d2d20, dFE2d2d02;
	double dFE2d3d00, dFE2d3d10, dFE2d3d01, dFE2d3d20, dFE2d3d02;
	double dFE2d4d00, dFE2d4d10, dFE2d4d01, dFE2d4d20, dFE2d4d02;
	double dFE2d5d00, dFE2d5d10, dFE2d5d01, dFE2d5d20, dFE2d5d02;
	double dFE2d6d00, dFE2d6d10, dFE2d6d01, dFE2d6d20, dFE2d6d02;
	double dFE2d7d00, dFE2d7d10, dFE2d7d01, dFE2d7d20, dFE2d7d02;
	double dFE2d8d00, dFE2d8d10, dFE2d8d01, dFE2d8d20, dFE2d8d02;
	double dFE2d9d00, dFE2d9d10, dFE2d9d01, dFE2d9d20, dFE2d9d02;
	
	double FE3;
	double dFE3d0d00, dFE3d0d10, dFE3d0d01, dFE3d0d20, dFE3d0d02;
	double dFE3d1d00, dFE3d1d10, dFE3d1d01, dFE3d1d20, dFE3d1d02;
	double dFE3d2d00, dFE3d2d10, dFE3d2d01, dFE3d2d20, dFE3d2d02;
	double dFE3d3d00, dFE3d3d10, dFE3d3d01, dFE3d3d20, dFE3d3d02;
	double dFE3d4d00, dFE3d4d10, dFE3d4d01, dFE3d4d20, dFE3d4d02;
	double dFE3d5d00, dFE3d5d10, dFE3d5d01, dFE3d5d20, dFE3d5d02;
	double dFE3d6d00, dFE3d6d10, dFE3d6d01, dFE3d6d20, dFE3d6d02;
	double dFE3d7d00, dFE3d7d10, dFE3d7d01, dFE3d7d20, dFE3d7d02;
	double dFE3d8d00, dFE3d8d10, dFE3d8d01, dFE3d8d20, dFE3d8d02;
	double dFE3d9d00, dFE3d9d10, dFE3d9d01, dFE3d9d20, dFE3d9d02;
	
	double FE4;
	double dFE4d0d00, dFE4d0d10, dFE4d0d01, dFE4d0d20, dFE4d0d02;
	double dFE4d1d00, dFE4d1d10, dFE4d1d01, dFE4d1d20, dFE4d1d02;
	double dFE4d2d00, dFE4d2d10, dFE4d2d01, dFE4d2d20, dFE4d2d02;
	double dFE4d3d00, dFE4d3d10, dFE4d3d01, dFE4d3d20, dFE4d3d02;
	double dFE4d4d00, dFE4d4d10, dFE4d4d01, dFE4d4d20, dFE4d4d02;
	double dFE4d5d00, dFE4d5d10, dFE4d5d01, dFE4d5d20, dFE4d5d02;
	double dFE4d6d00, dFE4d6d10, dFE4d6d01, dFE4d6d20, dFE4d6d02;
	double dFE4d7d00, dFE4d7d10, dFE4d7d01, dFE4d7d20, dFE4d7d02;
	double dFE4d8d00, dFE4d8d10, dFE4d8d01, dFE4d8d20, dFE4d8d02;
	double dFE4d9d00, dFE4d9d10, dFE4d9d01, dFE4d9d20, dFE4d9d02;
	
	double FE5;
	double dFE5d0d00, dFE5d0d10, dFE5d0d01, dFE5d0d20, dFE5d0d02;
	double dFE5d1d00, dFE5d1d10, dFE5d1d01, dFE5d1d20, dFE5d1d02;
	double dFE5d2d00, dFE5d2d10, dFE5d2d01, dFE5d2d20, dFE5d2d02;
	double dFE5d3d00, dFE5d3d10, dFE5d3d01, dFE5d3d20, dFE5d3d02;
	double dFE5d4d00, dFE5d4d10, dFE5d4d01, dFE5d4d20, dFE5d4d02;
	double dFE5d5d00, dFE5d5d10, dFE5d5d01, dFE5d5d20, dFE5d5d02;
	double dFE5d6d00, dFE5d6d10, dFE5d6d01, dFE5d6d20, dFE5d6d02;
	double dFE5d7d00, dFE5d7d10, dFE5d7d01, dFE5d7d20, dFE5d7d02;
	double dFE5d8d00, dFE5d8d10, dFE5d8d01, dFE5d8d20, dFE5d8d02;
	double dFE5d9d00, dFE5d9d10, dFE5d9d01, dFE5d9d20, dFE5d9d02;
	
	double FE6;
	double dFE6d0d00, dFE6d0d10, dFE6d0d01, dFE6d0d20, dFE6d0d02;
	double dFE6d1d00, dFE6d1d10, dFE6d1d01, dFE6d1d20, dFE6d1d02;
	double dFE6d2d00, dFE6d2d10, dFE6d2d01, dFE6d2d20, dFE6d2d02;
	double dFE6d3d00, dFE6d3d10, dFE6d3d01, dFE6d3d20, dFE6d3d02;
	double dFE6d4d00, dFE6d4d10, dFE6d4d01, dFE6d4d20, dFE6d4d02;
	double dFE6d5d00, dFE6d5d10, dFE6d5d01, dFE6d5d20, dFE6d5d02;
	double dFE6d6d00, dFE6d6d10, dFE6d6d01, dFE6d6d20, dFE6d6d02;
	double dFE6d7d00, dFE6d7d10, dFE6d7d01, dFE6d7d20, dFE6d7d02;
	double dFE6d8d00, dFE6d8d10, dFE6d8d01, dFE6d8d20, dFE6d8d02;
	double dFE6d9d00, dFE6d9d10, dFE6d9d01, dFE6d9d20, dFE6d9d02;
	
	double FE7;
	double dFE7d0d00, dFE7d0d10, dFE7d0d01, dFE7d0d20, dFE7d0d02;
	double dFE7d1d00, dFE7d1d10, dFE7d1d01, dFE7d1d20, dFE7d1d02;
	double dFE7d2d00, dFE7d2d10, dFE7d2d01, dFE7d2d20, dFE7d2d02;
	double dFE7d3d00, dFE7d3d10, dFE7d3d01, dFE7d3d20, dFE7d3d02;
	double dFE7d4d00, dFE7d4d10, dFE7d4d01, dFE7d4d20, dFE7d4d02;
	double dFE7d5d00, dFE7d5d10, dFE7d5d01, dFE7d5d20, dFE7d5d02;
	double dFE7d6d00, dFE7d6d10, dFE7d6d01, dFE7d6d20, dFE7d6d02;
	double dFE7d7d00, dFE7d7d10, dFE7d7d01, dFE7d7d20, dFE7d7d02;
	double dFE7d8d00, dFE7d8d10, dFE7d8d01, dFE7d8d20, dFE7d8d02;
	double dFE7d9d00, dFE7d9d10, dFE7d9d01, dFE7d9d20, dFE7d9d02;
	
	double FE8;
	double dFE8d0d00, dFE8d0d10, dFE8d0d01, dFE8d0d20, dFE8d0d02;
	double dFE8d1d00, dFE8d1d10, dFE8d1d01, dFE8d1d20, dFE8d1d02;
	double dFE8d2d00, dFE8d2d10, dFE8d2d01, dFE8d2d20, dFE8d2d02;
	double dFE8d3d00, dFE8d3d10, dFE8d3d01, dFE8d3d20, dFE8d3d02;
	double dFE8d4d00, dFE8d4d10, dFE8d4d01, dFE8d4d20, dFE8d4d02;
	double dFE8d5d00, dFE8d5d10, dFE8d5d01, dFE8d5d20, dFE8d5d02;
	double dFE8d6d00, dFE8d6d10, dFE8d6d01, dFE8d6d20, dFE8d6d02;
	double dFE8d7d00, dFE8d7d10, dFE8d7d01, dFE8d7d20, dFE8d7d02;
	double dFE8d8d00, dFE8d8d10, dFE8d8d01, dFE8d8d20, dFE8d8d02;
	double dFE8d9d00, dFE8d9d10, dFE8d9d01, dFE8d9d20, dFE8d9d02;
	
	double FE9;
	double dFE9d0d00, dFE9d0d10, dFE9d0d01, dFE9d0d20, dFE9d0d02;
	double dFE9d1d00, dFE9d1d10, dFE9d1d01, dFE9d1d20, dFE9d1d02;
	double dFE9d2d00, dFE9d2d10, dFE9d2d01, dFE9d2d20, dFE9d2d02;
	double dFE9d3d00, dFE9d3d10, dFE9d3d01, dFE9d3d20, dFE9d3d02;
	double dFE9d4d00, dFE9d4d10, dFE9d4d01, dFE9d4d20, dFE9d4d02;
	double dFE9d5d00, dFE9d5d10, dFE9d5d01, dFE9d5d20, dFE9d5d02;
	double dFE9d6d00, dFE9d6d10, dFE9d6d01, dFE9d6d20, dFE9d6d02;
	double dFE9d7d00, dFE9d7d10, dFE9d7d01, dFE9d7d20, dFE9d7d02;
	double dFE9d8d00, dFE9d8d10, dFE9d8d01, dFE9d8d20, dFE9d8d02;
	double dFE9d9d00, dFE9d9d10, dFE9d9d01, dFE9d9d20, dFE9d9d02;
	
	//----- Rezero system -----
	nnztemp = 0;
	for (i = 0; i < N; ++i)
	{
		bvec[i]  = 0.0;
		Dxvec[i] = 0.0;
		Dyvec[i] = 0.0;
	}
	for (k = 0; k < nnzJac; ++k)
	{
		JacVal[k] = 0.0;
		JacRow[k] = 0;
		JacCol[k] = 0;
	}
	
	nnztemp = BVP_GENBC(params, x, y, Psi, JacVal, JacRow, JacCol, bvec, Dxvec, Dyvec);
	
	
	//:::::::::::::::::::: Main Loop ::::::::::::::::::::
	for (i = 1; i < m-1; ++i)
	{
		for (j = 1; j < n-1; ++j)
		{
			k = (i*n+j)*p;
			I = i, J = j;
			xk = x[j];
			yk = y[i];
			
			//Zero derivatives of each field
			onexside = 0; twoxside = 0;
			oneyside = 0; twoyside = 0;
			u0d00k  = 0.0, u0d10k  = 0.0, u0d20k  = 0.0, u0d01k  = 0.0, u0d02k  = 0.0;
			u0d00dk = 0.0, u0d10dk = 0.0, u0d20dk = 0.0, u0d01dk = 0.0, u0d02dk = 0.0;
			
			u1d00k  = 0.0, u1d10k  = 0.0, u1d20k  = 0.0, u1d01k  = 0.0, u1d02k  = 0.0;
			u1d00dk = 0.0, u1d10dk = 0.0, u1d20dk = 0.0, u1d01dk = 0.0, u1d02dk = 0.0;
			
			u2d00k  = 0.0, u2d10k  = 0.0, u2d20k  = 0.0, u2d01k  = 0.0, u2d02k  = 0.0;
			u2d00dk = 0.0, u2d10dk = 0.0, u2d20dk = 0.0, u2d01dk = 0.0, u2d02dk = 0.0;
			
			u3d00k  = 0.0, u3d10k  = 0.0, u3d20k  = 0.0, u3d01k  = 0.0, u3d02k  = 0.0;
			u3d00dk = 0.0, u3d10dk = 0.0, u3d20dk = 0.0, u3d01dk = 0.0, u3d02dk = 0.0;
			
			u4d00k  = 0.0, u4d10k  = 0.0, u4d20k  = 0.0, u4d01k  = 0.0, u4d02k  = 0.0;
			u4d00dk = 0.0, u4d10dk = 0.0, u4d20dk = 0.0, u4d01dk = 0.0, u4d02dk = 0.0;
			
			u5d00k  = 0.0, u5d10k  = 0.0, u5d20k  = 0.0, u5d01k  = 0.0, u5d02k  = 0.0;
			u5d00dk = 0.0, u5d10dk = 0.0, u5d20dk = 0.0, u5d01dk = 0.0, u5d02dk = 0.0;
			
			u6d00k  = 0.0, u6d10k  = 0.0, u6d20k  = 0.0, u6d01k  = 0.0, u6d02k  = 0.0;
			u6d00dk = 0.0, u6d10dk = 0.0, u6d20dk = 0.0, u6d01dk = 0.0, u6d02dk = 0.0;
			
			u7d00k  = 0.0, u7d10k  = 0.0, u7d20k  = 0.0, u7d01k  = 0.0, u7d02k  = 0.0;
			u7d00dk = 0.0, u7d10dk = 0.0, u7d20dk = 0.0, u7d01dk = 0.0, u7d02dk = 0.0;
			
			u8d00k  = 0.0, u8d10k  = 0.0, u8d20k  = 0.0, u8d01k  = 0.0, u8d02k  = 0.0;
			u8d00dk = 0.0, u8d10dk = 0.0, u8d20dk = 0.0, u8d01dk = 0.0, u8d02dk = 0.0;
			
			u9d00k  = 0.0, u9d10k  = 0.0, u9d20k  = 0.0, u9d01k  = 0.0, u9d02k  = 0.0;
			u9d00dk = 0.0, u9d10dk = 0.0, u9d20dk = 0.0, u9d01dk = 0.0, u9d02dk = 0.0;
			
			//Compute stencil
			status = steninit(sy, &oneyside, &twoyside, I, r, m);
			status = steninit(sx, &onexside, &twoxside, J, r, n);
			
			//Compute Newton polynomial representation for each field
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+0, Psi, 0, b0d0, b0d1, b0d2, &u0d00dk, &u0d01dk, &u0d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+0, Psi, 0, a0d0, a0d1, a0d2, &u0d00dk, &u0d10dk, &u0d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+1, Psi, 1, b1d0, b1d1, b1d2, &u1d00dk, &u1d01dk, &u1d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+1, Psi, 1, a1d0, a1d1, a1d2, &u1d00dk, &u1d10dk, &u1d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+2, Psi, 2, b2d0, b2d1, b2d2, &u2d00dk, &u2d01dk, &u2d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+2, Psi, 2, a2d0, a2d1, a2d2, &u2d00dk, &u2d10dk, &u2d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+3, Psi, 3, b3d0, b3d1, b3d2, &u3d00dk, &u3d01dk, &u3d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+3, Psi, 3, a3d0, a3d1, a3d2, &u3d00dk, &u3d10dk, &u3d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+4, Psi, 4, b4d0, b4d1, b4d2, &u4d00dk, &u4d01dk, &u4d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+4, Psi, 4, a4d0, a4d1, a4d2, &u4d00dk, &u4d10dk, &u4d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+5, Psi, 5, b5d0, b5d1, b5d2, &u5d00dk, &u5d01dk, &u5d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+5, Psi, 5, a5d0, a5d1, a5d2, &u5d00dk, &u5d10dk, &u5d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+6, Psi, 6, b6d0, b6d1, b6d2, &u6d00dk, &u6d01dk, &u6d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+6, Psi, 6, a6d0, a6d1, a6d2, &u6d00dk, &u6d10dk, &u6d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+7, Psi, 7, b7d0, b7d1, b7d2, &u7d00dk, &u7d01dk, &u7d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+7, Psi, 7, a7d0, a7d1, a7d2, &u7d00dk, &u7d10dk, &u7d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+8, Psi, 8, b8d0, b8d1, b8d2, &u8d00dk, &u8d01dk, &u8d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+8, Psi, 8, a8d0, a8d1, a8d2, &u8d00dk, &u8d10dk, &u8d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+9, Psi, 9, b9d0, b9d1, b9d2, &u9d00dk, &u9d01dk, &u9d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+9, Psi, 9, a9d0, a9d1, a9d2, &u9d00dk, &u9d10dk, &u9d20dk, xk);
			
			//Calculate derivative of each field
			for (b = 0; b <= r+oneyside; ++b)
			{
				u0d01k += Psi[((I+sy[b])*n +J)*p +0]*b0d1[b];
				u0d02k += Psi[((I+sy[b])*n +J)*p +0]*b0d2[b];
				
				u1d01k += Psi[((I+sy[b])*n +J)*p +1]*b1d1[b];
				u1d02k += Psi[((I+sy[b])*n +J)*p +1]*b1d2[b];
				
				u2d01k += Psi[((I+sy[b])*n +J)*p +2]*b2d1[b];
				u2d02k += Psi[((I+sy[b])*n +J)*p +2]*b2d2[b];
				
				u3d01k += Psi[((I+sy[b])*n +J)*p +3]*b3d1[b];
				u3d02k += Psi[((I+sy[b])*n +J)*p +3]*b3d2[b];
				
				u4d01k += Psi[((I+sy[b])*n +J)*p +4]*b4d1[b];
				u4d02k += Psi[((I+sy[b])*n +J)*p +4]*b4d2[b];
				
				u5d01k += Psi[((I+sy[b])*n +J)*p +5]*b5d1[b];
				u5d02k += Psi[((I+sy[b])*n +J)*p +5]*b5d2[b];
				
				u6d01k += Psi[((I+sy[b])*n +J)*p +6]*b6d1[b];
				u6d02k += Psi[((I+sy[b])*n +J)*p +6]*b6d2[b];
				
				u7d01k += Psi[((I+sy[b])*n +J)*p +7]*b7d1[b];
				u7d02k += Psi[((I+sy[b])*n +J)*p +7]*b7d2[b];
				
				u8d01k += Psi[((I+sy[b])*n +J)*p +8]*b8d1[b];
				u8d02k += Psi[((I+sy[b])*n +J)*p +8]*b8d2[b];
				
				u9d01k += Psi[((I+sy[b])*n +J)*p +9]*b9d1[b];
				u9d02k += Psi[((I+sy[b])*n +J)*p +9]*b9d2[b];
				
			}
			for (a = 0; a <= r+onexside; ++a)
			{
				u0d10k += Psi[(I*n +J+sx[a])*p +0]*a0d1[a];
				u0d20k += Psi[(I*n +J+sx[a])*p +0]*a0d2[a];
				
				u1d10k += Psi[(I*n +J+sx[a])*p +1]*a1d1[a];
				u1d20k += Psi[(I*n +J+sx[a])*p +1]*a1d2[a];
				
				u2d10k += Psi[(I*n +J+sx[a])*p +2]*a2d1[a];
				u2d20k += Psi[(I*n +J+sx[a])*p +2]*a2d2[a];
				
				u3d10k += Psi[(I*n +J+sx[a])*p +3]*a3d1[a];
				u3d20k += Psi[(I*n +J+sx[a])*p +3]*a3d2[a];
				
				u4d10k += Psi[(I*n +J+sx[a])*p +4]*a4d1[a];
				u4d20k += Psi[(I*n +J+sx[a])*p +4]*a4d2[a];
				
				u5d10k += Psi[(I*n +J+sx[a])*p +5]*a5d1[a];
				u5d20k += Psi[(I*n +J+sx[a])*p +5]*a5d2[a];
				
				u6d10k += Psi[(I*n +J+sx[a])*p +6]*a6d1[a];
				u6d20k += Psi[(I*n +J+sx[a])*p +6]*a6d2[a];
				
				u7d10k += Psi[(I*n +J+sx[a])*p +7]*a7d1[a];
				u7d20k += Psi[(I*n +J+sx[a])*p +7]*a7d2[a];
				
				u8d10k += Psi[(I*n +J+sx[a])*p +8]*a8d1[a];
				u8d20k += Psi[(I*n +J+sx[a])*p +8]*a8d2[a];
				
				u9d10k += Psi[(I*n +J+sx[a])*p +9]*a9d1[a];
				u9d20k += Psi[(I*n +J+sx[a])*p +9]*a9d2[a];
				
			}
			
			for (a = 0; a <= r+onexside; ++a)
			{
				u0d00k += Psi[(I*n +J+sx[a])*p +0]*a0d0[a];
				u1d00k += Psi[(I*n +J+sx[a])*p +1]*a1d0[a];
				u2d00k += Psi[(I*n +J+sx[a])*p +2]*a2d0[a];
				u3d00k += Psi[(I*n +J+sx[a])*p +3]*a3d0[a];
				u4d00k += Psi[(I*n +J+sx[a])*p +4]*a4d0[a];
				u5d00k += Psi[(I*n +J+sx[a])*p +5]*a5d0[a];
				u6d00k += Psi[(I*n +J+sx[a])*p +6]*a6d0[a];
				u7d00k += Psi[(I*n +J+sx[a])*p +7]*a7d0[a];
				u8d00k += Psi[(I*n +J+sx[a])*p +8]*a8d0[a];
				u9d00k += Psi[(I*n +J+sx[a])*p +9]*a9d0[a];
			}
			
			//---------- Evaluate partial derivatives and structure for Jacobian ----------
			
			//FE0
			FE0       =       FE0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d0d00 = dFE0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d0d10 = dFE0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d0d01 = dFE0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d0d20 = dFE0d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d0d02 = dFE0d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d1d00 = dFE0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d1d10 = dFE0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d1d01 = dFE0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d1d20 = dFE0d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d1d02 = dFE0d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d2d00 = dFE0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d2d10 = dFE0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d2d01 = dFE0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d2d20 = dFE0d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d2d02 = dFE0d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d3d00 = dFE0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d3d10 = dFE0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d3d01 = dFE0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d3d20 = dFE0d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d3d02 = dFE0d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d4d00 = dFE0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d4d10 = dFE0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d4d01 = dFE0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d4d20 = dFE0d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d4d02 = dFE0d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d5d00 = dFE0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d5d10 = dFE0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d5d01 = dFE0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d5d20 = dFE0d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d5d02 = dFE0d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d6d00 = dFE0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d6d10 = dFE0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d6d01 = dFE0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d6d20 = dFE0d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d6d02 = dFE0d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d7d00 = dFE0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d7d10 = dFE0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d7d01 = dFE0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d7d20 = dFE0d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d7d02 = dFE0d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d8d00 = dFE0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d8d10 = dFE0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d8d01 = dFE0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d8d20 = dFE0d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d8d02 = dFE0d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE0d9d00 = dFE0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d9d10 = dFE0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d9d01 = dFE0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d9d20 = dFE0d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE0d9d02 = dFE0d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 0;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE0d0d00*a0d0[a] + dFE0d0d10*a0d1[a] + dFE0d0d20*a0d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d1d00*a1d0[a] + dFE0d1d10*a1d1[a] + dFE0d1d20*a1d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d2d00*a2d0[a] + dFE0d2d10*a2d1[a] + dFE0d2d20*a2d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d3d00*a3d0[a] + dFE0d3d10*a3d1[a] + dFE0d3d20*a3d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d4d00*a4d0[a] + dFE0d4d10*a4d1[a] + dFE0d4d20*a4d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d5d00*a5d0[a] + dFE0d5d10*a5d1[a] + dFE0d5d20*a5d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d7d00*a7d0[a] + dFE0d7d10*a7d1[a] + dFE0d7d20*a7d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d9d00*a9d0[a] + dFE0d9d10*a9d1[a] + dFE0d9d20*a9d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE0d0d01*b0d1[b] + dFE0d0d02*b0d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d1d01*b1d1[b] + dFE0d1d02*b1d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d2d01*b2d1[b] + dFE0d2d02*b2d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d3d01*b3d1[b] + dFE0d3d02*b3d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE0d4d01*b4d1[b] + dFE0d4d02*b4d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE0);
			Dxvec[k +pn] = dFE0d0d00*u0d00dk + dFE0d0d10*u0d10dk + dFE0d0d20*u0d20dk + dFE0d1d00*u1d00dk + dFE0d1d10*u1d10dk + dFE0d1d20*u1d20dk + dFE0d2d00*u2d00dk + dFE0d2d10*u2d10dk + dFE0d2d20*u2d20dk + dFE0d3d00*u3d00dk + dFE0d3d10*u3d10dk + dFE0d3d20*u3d20dk + dFE0d4d00*u4d00dk + dFE0d4d10*u4d10dk + dFE0d4d20*u4d20dk + dFE0d5d00*u5d00dk + dFE0d5d10*u5d10dk + dFE0d5d20*u5d20dk + dFE0d6d00*u6d00dk + dFE0d6d10*u6d10dk + dFE0d6d20*u6d20dk + dFE0d7d00*u7d00dk + dFE0d7d10*u7d10dk + dFE0d7d20*u7d20dk + dFE0d8d00*u8d00dk + dFE0d8d10*u8d10dk + dFE0d8d20*u8d20dk + dFE0d9d00*u9d00dk + dFE0d9d10*u9d10dk + dFE0d9d20*u9d20dk;
			Dyvec[k +pn] = dFE0d0d00*u0d00dk + dFE0d0d01*u0d01dk + dFE0d0d02*u0d02dk + dFE0d1d00*u1d00dk + dFE0d1d01*u1d01dk + dFE0d1d02*u1d02dk + dFE0d2d00*u2d00dk + dFE0d2d01*u2d01dk + dFE0d2d02*u2d02dk + dFE0d3d00*u3d00dk + dFE0d3d01*u3d01dk + dFE0d3d02*u3d02dk + dFE0d4d00*u4d00dk + dFE0d4d01*u4d01dk + dFE0d4d02*u4d02dk + dFE0d5d00*u5d00dk + dFE0d5d01*u5d01dk + dFE0d5d02*u5d02dk + dFE0d6d00*u6d00dk + dFE0d6d01*u6d01dk + dFE0d6d02*u6d02dk + dFE0d7d00*u7d00dk + dFE0d7d01*u7d01dk + dFE0d7d02*u7d02dk + dFE0d8d00*u8d00dk + dFE0d8d01*u8d01dk + dFE0d8d02*u8d02dk + dFE0d9d00*u9d00dk + dFE0d9d01*u9d01dk + dFE0d9d02*u9d02dk;
			
			//FE1
			FE1       =       FE1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d0d00 = dFE1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d0d10 = dFE1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d0d01 = dFE1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d0d20 = dFE1d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d0d02 = dFE1d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d1d00 = dFE1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d1d10 = dFE1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d1d01 = dFE1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d1d20 = dFE1d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d1d02 = dFE1d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d2d00 = dFE1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d2d10 = dFE1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d2d01 = dFE1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d2d20 = dFE1d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d2d02 = dFE1d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d3d00 = dFE1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d3d10 = dFE1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d3d01 = dFE1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d3d20 = dFE1d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d3d02 = dFE1d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d4d00 = dFE1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d4d10 = dFE1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d4d01 = dFE1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d4d20 = dFE1d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d4d02 = dFE1d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d5d00 = dFE1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d5d10 = dFE1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d5d01 = dFE1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d5d20 = dFE1d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d5d02 = dFE1d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d6d00 = dFE1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d6d10 = dFE1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d6d01 = dFE1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d6d20 = dFE1d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d6d02 = dFE1d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d7d00 = dFE1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d7d10 = dFE1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d7d01 = dFE1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d7d20 = dFE1d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d7d02 = dFE1d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d8d00 = dFE1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d8d10 = dFE1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d8d01 = dFE1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d8d20 = dFE1d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d8d02 = dFE1d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE1d9d00 = dFE1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d9d10 = dFE1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d9d01 = dFE1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d9d20 = dFE1d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE1d9d02 = dFE1d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 1;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE1d0d00*a0d0[a] + dFE1d0d10*a0d1[a] + dFE1d0d20*a0d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d1d00*a1d0[a] + dFE1d1d10*a1d1[a] + dFE1d1d20*a1d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d2d00*a2d0[a] + dFE1d2d10*a2d1[a] + dFE1d2d20*a2d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d3d00*a3d0[a] + dFE1d3d10*a3d1[a] + dFE1d3d20*a3d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d4d00*a4d0[a] + dFE1d4d10*a4d1[a] + dFE1d4d20*a4d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d5d00*a5d0[a] + dFE1d5d10*a5d1[a] + dFE1d5d20*a5d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d9d00*a9d0[a] + dFE1d9d10*a9d1[a] + dFE1d9d20*a9d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE1d0d01*b0d1[b] + dFE1d0d02*b0d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d1d01*b1d1[b] + dFE1d1d02*b1d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d2d01*b2d1[b] + dFE1d2d02*b2d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d3d01*b3d1[b] + dFE1d3d02*b3d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE1d4d01*b4d1[b] + dFE1d4d02*b4d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE1);
			Dxvec[k +pn] = dFE1d0d00*u0d00dk + dFE1d0d10*u0d10dk + dFE1d0d20*u0d20dk + dFE1d1d00*u1d00dk + dFE1d1d10*u1d10dk + dFE1d1d20*u1d20dk + dFE1d2d00*u2d00dk + dFE1d2d10*u2d10dk + dFE1d2d20*u2d20dk + dFE1d3d00*u3d00dk + dFE1d3d10*u3d10dk + dFE1d3d20*u3d20dk + dFE1d4d00*u4d00dk + dFE1d4d10*u4d10dk + dFE1d4d20*u4d20dk + dFE1d5d00*u5d00dk + dFE1d5d10*u5d10dk + dFE1d5d20*u5d20dk + dFE1d6d00*u6d00dk + dFE1d6d10*u6d10dk + dFE1d6d20*u6d20dk + dFE1d7d00*u7d00dk + dFE1d7d10*u7d10dk + dFE1d7d20*u7d20dk + dFE1d8d00*u8d00dk + dFE1d8d10*u8d10dk + dFE1d8d20*u8d20dk + dFE1d9d00*u9d00dk + dFE1d9d10*u9d10dk + dFE1d9d20*u9d20dk;
			Dyvec[k +pn] = dFE1d0d00*u0d00dk + dFE1d0d01*u0d01dk + dFE1d0d02*u0d02dk + dFE1d1d00*u1d00dk + dFE1d1d01*u1d01dk + dFE1d1d02*u1d02dk + dFE1d2d00*u2d00dk + dFE1d2d01*u2d01dk + dFE1d2d02*u2d02dk + dFE1d3d00*u3d00dk + dFE1d3d01*u3d01dk + dFE1d3d02*u3d02dk + dFE1d4d00*u4d00dk + dFE1d4d01*u4d01dk + dFE1d4d02*u4d02dk + dFE1d5d00*u5d00dk + dFE1d5d01*u5d01dk + dFE1d5d02*u5d02dk + dFE1d6d00*u6d00dk + dFE1d6d01*u6d01dk + dFE1d6d02*u6d02dk + dFE1d7d00*u7d00dk + dFE1d7d01*u7d01dk + dFE1d7d02*u7d02dk + dFE1d8d00*u8d00dk + dFE1d8d01*u8d01dk + dFE1d8d02*u8d02dk + dFE1d9d00*u9d00dk + dFE1d9d01*u9d01dk + dFE1d9d02*u9d02dk;
			
			//FE2
			FE2       =       FE2out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d0d00 = dFE2d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d0d10 = dFE2d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d0d01 = dFE2d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d0d20 = dFE2d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d0d02 = dFE2d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d1d00 = dFE2d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d1d10 = dFE2d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d1d01 = dFE2d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d1d20 = dFE2d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d1d02 = dFE2d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d2d00 = dFE2d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d2d10 = dFE2d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d2d01 = dFE2d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d2d20 = dFE2d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d2d02 = dFE2d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d3d00 = dFE2d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d3d10 = dFE2d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d3d01 = dFE2d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d3d20 = dFE2d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d3d02 = dFE2d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d4d00 = dFE2d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d4d10 = dFE2d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d4d01 = dFE2d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d4d20 = dFE2d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d4d02 = dFE2d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d5d00 = dFE2d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d5d10 = dFE2d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d5d01 = dFE2d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d5d20 = dFE2d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d5d02 = dFE2d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d6d00 = dFE2d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d6d10 = dFE2d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d6d01 = dFE2d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d6d20 = dFE2d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d6d02 = dFE2d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d7d00 = dFE2d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d7d10 = dFE2d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d7d01 = dFE2d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d7d20 = dFE2d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d7d02 = dFE2d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d8d00 = dFE2d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d8d10 = dFE2d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d8d01 = dFE2d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d8d20 = dFE2d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d8d02 = dFE2d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE2d9d00 = dFE2d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d9d10 = dFE2d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d9d01 = dFE2d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d9d20 = dFE2d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE2d9d02 = dFE2d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 2;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE2d0d00*a0d0[a] + dFE2d0d10*a0d1[a] + dFE2d0d20*a0d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d1d00*a1d0[a] + dFE2d1d10*a1d1[a] + dFE2d1d20*a1d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d2d00*a2d0[a] + dFE2d2d10*a2d1[a] + dFE2d2d20*a2d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d3d00*a3d0[a] + dFE2d3d10*a3d1[a] + dFE2d3d20*a3d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d4d00*a4d0[a] + dFE2d4d10*a4d1[a] + dFE2d4d20*a4d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE2d0d01*b0d1[b] + dFE2d0d02*b0d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d2d01*b2d1[b] + dFE2d2d02*b2d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d3d01*b3d1[b] + dFE2d3d02*b3d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE2d4d01*b4d1[b] + dFE2d4d02*b4d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE2);
			Dxvec[k +pn] = dFE2d0d00*u0d00dk + dFE2d0d10*u0d10dk + dFE2d0d20*u0d20dk + dFE2d1d00*u1d00dk + dFE2d1d10*u1d10dk + dFE2d1d20*u1d20dk + dFE2d2d00*u2d00dk + dFE2d2d10*u2d10dk + dFE2d2d20*u2d20dk + dFE2d3d00*u3d00dk + dFE2d3d10*u3d10dk + dFE2d3d20*u3d20dk + dFE2d4d00*u4d00dk + dFE2d4d10*u4d10dk + dFE2d4d20*u4d20dk + dFE2d5d00*u5d00dk + dFE2d5d10*u5d10dk + dFE2d5d20*u5d20dk + dFE2d6d00*u6d00dk + dFE2d6d10*u6d10dk + dFE2d6d20*u6d20dk + dFE2d7d00*u7d00dk + dFE2d7d10*u7d10dk + dFE2d7d20*u7d20dk + dFE2d8d00*u8d00dk + dFE2d8d10*u8d10dk + dFE2d8d20*u8d20dk + dFE2d9d00*u9d00dk + dFE2d9d10*u9d10dk + dFE2d9d20*u9d20dk;
			Dyvec[k +pn] = dFE2d0d00*u0d00dk + dFE2d0d01*u0d01dk + dFE2d0d02*u0d02dk + dFE2d1d00*u1d00dk + dFE2d1d01*u1d01dk + dFE2d1d02*u1d02dk + dFE2d2d00*u2d00dk + dFE2d2d01*u2d01dk + dFE2d2d02*u2d02dk + dFE2d3d00*u3d00dk + dFE2d3d01*u3d01dk + dFE2d3d02*u3d02dk + dFE2d4d00*u4d00dk + dFE2d4d01*u4d01dk + dFE2d4d02*u4d02dk + dFE2d5d00*u5d00dk + dFE2d5d01*u5d01dk + dFE2d5d02*u5d02dk + dFE2d6d00*u6d00dk + dFE2d6d01*u6d01dk + dFE2d6d02*u6d02dk + dFE2d7d00*u7d00dk + dFE2d7d01*u7d01dk + dFE2d7d02*u7d02dk + dFE2d8d00*u8d00dk + dFE2d8d01*u8d01dk + dFE2d8d02*u8d02dk + dFE2d9d00*u9d00dk + dFE2d9d01*u9d01dk + dFE2d9d02*u9d02dk;
			
			//FE3
			FE3       =       FE3out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d0d00 = dFE3d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d0d10 = dFE3d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d0d01 = dFE3d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d0d20 = dFE3d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d0d02 = dFE3d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d1d00 = dFE3d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d1d10 = dFE3d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d1d01 = dFE3d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d1d20 = dFE3d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d1d02 = dFE3d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d2d00 = dFE3d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d2d10 = dFE3d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d2d01 = dFE3d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d2d20 = dFE3d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d2d02 = dFE3d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d3d00 = dFE3d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d3d10 = dFE3d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d3d01 = dFE3d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d3d20 = dFE3d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d3d02 = dFE3d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d4d00 = dFE3d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d4d10 = dFE3d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d4d01 = dFE3d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d4d20 = dFE3d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d4d02 = dFE3d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d5d00 = dFE3d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d5d10 = dFE3d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d5d01 = dFE3d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d5d20 = dFE3d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d5d02 = dFE3d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d6d00 = dFE3d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d6d10 = dFE3d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d6d01 = dFE3d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d6d20 = dFE3d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d6d02 = dFE3d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d7d00 = dFE3d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d7d10 = dFE3d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d7d01 = dFE3d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d7d20 = dFE3d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d7d02 = dFE3d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d8d00 = dFE3d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d8d10 = dFE3d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d8d01 = dFE3d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d8d20 = dFE3d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d8d02 = dFE3d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE3d9d00 = dFE3d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d9d10 = dFE3d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d9d01 = dFE3d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d9d20 = dFE3d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE3d9d02 = dFE3d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 3;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE3d0d00*a0d0[a] + dFE3d0d10*a0d1[a] + dFE3d0d20*a0d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d1d00*a1d0[a] + dFE3d1d10*a1d1[a] + dFE3d1d20*a1d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d2d00*a2d0[a] + dFE3d2d10*a2d1[a] + dFE3d2d20*a2d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d3d00*a3d0[a] + dFE3d3d10*a3d1[a] + dFE3d3d20*a3d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d4d00*a4d0[a] + dFE3d4d10*a4d1[a] + dFE3d4d20*a4d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d8d00*a8d0[a] + dFE3d8d10*a8d1[a] + dFE3d8d20*a8d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d9d00*a9d0[a] + dFE3d9d10*a9d1[a] + dFE3d9d20*a9d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE3d0d01*b0d1[b] + dFE3d0d02*b0d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d1d01*b1d1[b] + dFE3d1d02*b1d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d2d01*b2d1[b] + dFE3d2d02*b2d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d3d01*b3d1[b] + dFE3d3d02*b3d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE3d4d01*b4d1[b] + dFE3d4d02*b4d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE3);
			Dxvec[k +pn] = dFE3d0d00*u0d00dk + dFE3d0d10*u0d10dk + dFE3d0d20*u0d20dk + dFE3d1d00*u1d00dk + dFE3d1d10*u1d10dk + dFE3d1d20*u1d20dk + dFE3d2d00*u2d00dk + dFE3d2d10*u2d10dk + dFE3d2d20*u2d20dk + dFE3d3d00*u3d00dk + dFE3d3d10*u3d10dk + dFE3d3d20*u3d20dk + dFE3d4d00*u4d00dk + dFE3d4d10*u4d10dk + dFE3d4d20*u4d20dk + dFE3d5d00*u5d00dk + dFE3d5d10*u5d10dk + dFE3d5d20*u5d20dk + dFE3d6d00*u6d00dk + dFE3d6d10*u6d10dk + dFE3d6d20*u6d20dk + dFE3d7d00*u7d00dk + dFE3d7d10*u7d10dk + dFE3d7d20*u7d20dk + dFE3d8d00*u8d00dk + dFE3d8d10*u8d10dk + dFE3d8d20*u8d20dk + dFE3d9d00*u9d00dk + dFE3d9d10*u9d10dk + dFE3d9d20*u9d20dk;
			Dyvec[k +pn] = dFE3d0d00*u0d00dk + dFE3d0d01*u0d01dk + dFE3d0d02*u0d02dk + dFE3d1d00*u1d00dk + dFE3d1d01*u1d01dk + dFE3d1d02*u1d02dk + dFE3d2d00*u2d00dk + dFE3d2d01*u2d01dk + dFE3d2d02*u2d02dk + dFE3d3d00*u3d00dk + dFE3d3d01*u3d01dk + dFE3d3d02*u3d02dk + dFE3d4d00*u4d00dk + dFE3d4d01*u4d01dk + dFE3d4d02*u4d02dk + dFE3d5d00*u5d00dk + dFE3d5d01*u5d01dk + dFE3d5d02*u5d02dk + dFE3d6d00*u6d00dk + dFE3d6d01*u6d01dk + dFE3d6d02*u6d02dk + dFE3d7d00*u7d00dk + dFE3d7d01*u7d01dk + dFE3d7d02*u7d02dk + dFE3d8d00*u8d00dk + dFE3d8d01*u8d01dk + dFE3d8d02*u8d02dk + dFE3d9d00*u9d00dk + dFE3d9d01*u9d01dk + dFE3d9d02*u9d02dk;
			
			//FE4
			FE4       =       FE4out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d0d00 = dFE4d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d0d10 = dFE4d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d0d01 = dFE4d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d0d20 = dFE4d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d0d02 = dFE4d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d1d00 = dFE4d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d1d10 = dFE4d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d1d01 = dFE4d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d1d20 = dFE4d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d1d02 = dFE4d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d2d00 = dFE4d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d2d10 = dFE4d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d2d01 = dFE4d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d2d20 = dFE4d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d2d02 = dFE4d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d3d00 = dFE4d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d3d10 = dFE4d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d3d01 = dFE4d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d3d20 = dFE4d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d3d02 = dFE4d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d4d00 = dFE4d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d4d10 = dFE4d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d4d01 = dFE4d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d4d20 = dFE4d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d4d02 = dFE4d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d5d00 = dFE4d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d5d10 = dFE4d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d5d01 = dFE4d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d5d20 = dFE4d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d5d02 = dFE4d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d6d00 = dFE4d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d6d10 = dFE4d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d6d01 = dFE4d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d6d20 = dFE4d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d6d02 = dFE4d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d7d00 = dFE4d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d7d10 = dFE4d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d7d01 = dFE4d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d7d20 = dFE4d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d7d02 = dFE4d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d8d00 = dFE4d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d8d10 = dFE4d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d8d01 = dFE4d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d8d20 = dFE4d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d8d02 = dFE4d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE4d9d00 = dFE4d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d9d10 = dFE4d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d9d01 = dFE4d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d9d20 = dFE4d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE4d9d02 = dFE4d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 4;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE4d0d00*a0d0[a] + dFE4d0d10*a0d1[a] + dFE4d0d20*a0d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d1d00*a1d0[a] + dFE4d1d10*a1d1[a] + dFE4d1d20*a1d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d2d00*a2d0[a] + dFE4d2d10*a2d1[a] + dFE4d2d20*a2d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d3d00*a3d0[a] + dFE4d3d10*a3d1[a] + dFE4d3d20*a3d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d4d00*a4d0[a] + dFE4d4d10*a4d1[a] + dFE4d4d20*a4d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d5d00*a5d0[a] + dFE4d5d10*a5d1[a] + dFE4d5d20*a5d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d7d00*a7d0[a] + dFE4d7d10*a7d1[a] + dFE4d7d20*a7d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d8d00*a8d0[a] + dFE4d8d10*a8d1[a] + dFE4d8d20*a8d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE4d0d01*b0d1[b] + dFE4d0d02*b0d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d1d01*b1d1[b] + dFE4d1d02*b1d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d2d01*b2d1[b] + dFE4d2d02*b2d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d3d01*b3d1[b] + dFE4d3d02*b3d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
				
				JacVal[nnztemp] = dFE4d4d01*b4d1[b] + dFE4d4d02*b4d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE4);
			Dxvec[k +pn] = dFE4d0d00*u0d00dk + dFE4d0d10*u0d10dk + dFE4d0d20*u0d20dk + dFE4d1d00*u1d00dk + dFE4d1d10*u1d10dk + dFE4d1d20*u1d20dk + dFE4d2d00*u2d00dk + dFE4d2d10*u2d10dk + dFE4d2d20*u2d20dk + dFE4d3d00*u3d00dk + dFE4d3d10*u3d10dk + dFE4d3d20*u3d20dk + dFE4d4d00*u4d00dk + dFE4d4d10*u4d10dk + dFE4d4d20*u4d20dk + dFE4d5d00*u5d00dk + dFE4d5d10*u5d10dk + dFE4d5d20*u5d20dk + dFE4d6d00*u6d00dk + dFE4d6d10*u6d10dk + dFE4d6d20*u6d20dk + dFE4d7d00*u7d00dk + dFE4d7d10*u7d10dk + dFE4d7d20*u7d20dk + dFE4d8d00*u8d00dk + dFE4d8d10*u8d10dk + dFE4d8d20*u8d20dk + dFE4d9d00*u9d00dk + dFE4d9d10*u9d10dk + dFE4d9d20*u9d20dk;
			Dyvec[k +pn] = dFE4d0d00*u0d00dk + dFE4d0d01*u0d01dk + dFE4d0d02*u0d02dk + dFE4d1d00*u1d00dk + dFE4d1d01*u1d01dk + dFE4d1d02*u1d02dk + dFE4d2d00*u2d00dk + dFE4d2d01*u2d01dk + dFE4d2d02*u2d02dk + dFE4d3d00*u3d00dk + dFE4d3d01*u3d01dk + dFE4d3d02*u3d02dk + dFE4d4d00*u4d00dk + dFE4d4d01*u4d01dk + dFE4d4d02*u4d02dk + dFE4d5d00*u5d00dk + dFE4d5d01*u5d01dk + dFE4d5d02*u5d02dk + dFE4d6d00*u6d00dk + dFE4d6d01*u6d01dk + dFE4d6d02*u6d02dk + dFE4d7d00*u7d00dk + dFE4d7d01*u7d01dk + dFE4d7d02*u7d02dk + dFE4d8d00*u8d00dk + dFE4d8d01*u8d01dk + dFE4d8d02*u8d02dk + dFE4d9d00*u9d00dk + dFE4d9d01*u9d01dk + dFE4d9d02*u9d02dk;
			
			//FE5
			FE5       =       FE5out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d0d00 = dFE5d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d0d10 = dFE5d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d0d01 = dFE5d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d0d20 = dFE5d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d0d02 = dFE5d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d1d00 = dFE5d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d1d10 = dFE5d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d1d01 = dFE5d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d1d20 = dFE5d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d1d02 = dFE5d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d2d00 = dFE5d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d2d10 = dFE5d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d2d01 = dFE5d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d2d20 = dFE5d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d2d02 = dFE5d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d3d00 = dFE5d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d3d10 = dFE5d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d3d01 = dFE5d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d3d20 = dFE5d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d3d02 = dFE5d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d4d00 = dFE5d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d4d10 = dFE5d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d4d01 = dFE5d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d4d20 = dFE5d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d4d02 = dFE5d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d5d00 = dFE5d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d5d10 = dFE5d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d5d01 = dFE5d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d5d20 = dFE5d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d5d02 = dFE5d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d6d00 = dFE5d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d6d10 = dFE5d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d6d01 = dFE5d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d6d20 = dFE5d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d6d02 = dFE5d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d7d00 = dFE5d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d7d10 = dFE5d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d7d01 = dFE5d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d7d20 = dFE5d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d7d02 = dFE5d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d8d00 = dFE5d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d8d10 = dFE5d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d8d01 = dFE5d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d8d20 = dFE5d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d8d02 = dFE5d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE5d9d00 = dFE5d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d9d10 = dFE5d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d9d01 = dFE5d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d9d20 = dFE5d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE5d9d02 = dFE5d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 5;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE5d5d00*a5d0[a] + dFE5d5d10*a5d1[a] + dFE5d5d20*a5d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE5d0d01*b0d1[b] + dFE5d0d02*b0d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE5);
			Dxvec[k +pn] = dFE5d0d00*u0d00dk + dFE5d0d10*u0d10dk + dFE5d0d20*u0d20dk + dFE5d1d00*u1d00dk + dFE5d1d10*u1d10dk + dFE5d1d20*u1d20dk + dFE5d2d00*u2d00dk + dFE5d2d10*u2d10dk + dFE5d2d20*u2d20dk + dFE5d3d00*u3d00dk + dFE5d3d10*u3d10dk + dFE5d3d20*u3d20dk + dFE5d4d00*u4d00dk + dFE5d4d10*u4d10dk + dFE5d4d20*u4d20dk + dFE5d5d00*u5d00dk + dFE5d5d10*u5d10dk + dFE5d5d20*u5d20dk + dFE5d6d00*u6d00dk + dFE5d6d10*u6d10dk + dFE5d6d20*u6d20dk + dFE5d7d00*u7d00dk + dFE5d7d10*u7d10dk + dFE5d7d20*u7d20dk + dFE5d8d00*u8d00dk + dFE5d8d10*u8d10dk + dFE5d8d20*u8d20dk + dFE5d9d00*u9d00dk + dFE5d9d10*u9d10dk + dFE5d9d20*u9d20dk;
			Dyvec[k +pn] = dFE5d0d00*u0d00dk + dFE5d0d01*u0d01dk + dFE5d0d02*u0d02dk + dFE5d1d00*u1d00dk + dFE5d1d01*u1d01dk + dFE5d1d02*u1d02dk + dFE5d2d00*u2d00dk + dFE5d2d01*u2d01dk + dFE5d2d02*u2d02dk + dFE5d3d00*u3d00dk + dFE5d3d01*u3d01dk + dFE5d3d02*u3d02dk + dFE5d4d00*u4d00dk + dFE5d4d01*u4d01dk + dFE5d4d02*u4d02dk + dFE5d5d00*u5d00dk + dFE5d5d01*u5d01dk + dFE5d5d02*u5d02dk + dFE5d6d00*u6d00dk + dFE5d6d01*u6d01dk + dFE5d6d02*u6d02dk + dFE5d7d00*u7d00dk + dFE5d7d01*u7d01dk + dFE5d7d02*u7d02dk + dFE5d8d00*u8d00dk + dFE5d8d01*u8d01dk + dFE5d8d02*u8d02dk + dFE5d9d00*u9d00dk + dFE5d9d01*u9d01dk + dFE5d9d02*u9d02dk;
			
			//FE6
			FE6       =       FE6out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d0d00 = dFE6d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d0d10 = dFE6d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d0d01 = dFE6d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d0d20 = dFE6d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d0d02 = dFE6d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d1d00 = dFE6d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d1d10 = dFE6d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d1d01 = dFE6d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d1d20 = dFE6d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d1d02 = dFE6d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d2d00 = dFE6d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d2d10 = dFE6d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d2d01 = dFE6d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d2d20 = dFE6d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d2d02 = dFE6d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d3d00 = dFE6d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d3d10 = dFE6d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d3d01 = dFE6d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d3d20 = dFE6d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d3d02 = dFE6d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d4d00 = dFE6d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d4d10 = dFE6d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d4d01 = dFE6d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d4d20 = dFE6d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d4d02 = dFE6d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d5d00 = dFE6d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d5d10 = dFE6d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d5d01 = dFE6d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d5d20 = dFE6d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d5d02 = dFE6d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d6d00 = dFE6d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d6d10 = dFE6d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d6d01 = dFE6d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d6d20 = dFE6d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d6d02 = dFE6d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d7d00 = dFE6d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d7d10 = dFE6d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d7d01 = dFE6d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d7d20 = dFE6d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d7d02 = dFE6d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d8d00 = dFE6d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d8d10 = dFE6d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d8d01 = dFE6d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d8d20 = dFE6d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d8d02 = dFE6d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE6d9d00 = dFE6d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d9d10 = dFE6d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d9d01 = dFE6d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d9d20 = dFE6d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE6d9d02 = dFE6d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 6;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE6d6d00*a6d0[a] + dFE6d6d10*a6d1[a] + dFE6d6d20*a6d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +6;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE6d1d01*b1d1[b] + dFE6d1d02*b1d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE6);
			Dxvec[k +pn] = dFE6d0d00*u0d00dk + dFE6d0d10*u0d10dk + dFE6d0d20*u0d20dk + dFE6d1d00*u1d00dk + dFE6d1d10*u1d10dk + dFE6d1d20*u1d20dk + dFE6d2d00*u2d00dk + dFE6d2d10*u2d10dk + dFE6d2d20*u2d20dk + dFE6d3d00*u3d00dk + dFE6d3d10*u3d10dk + dFE6d3d20*u3d20dk + dFE6d4d00*u4d00dk + dFE6d4d10*u4d10dk + dFE6d4d20*u4d20dk + dFE6d5d00*u5d00dk + dFE6d5d10*u5d10dk + dFE6d5d20*u5d20dk + dFE6d6d00*u6d00dk + dFE6d6d10*u6d10dk + dFE6d6d20*u6d20dk + dFE6d7d00*u7d00dk + dFE6d7d10*u7d10dk + dFE6d7d20*u7d20dk + dFE6d8d00*u8d00dk + dFE6d8d10*u8d10dk + dFE6d8d20*u8d20dk + dFE6d9d00*u9d00dk + dFE6d9d10*u9d10dk + dFE6d9d20*u9d20dk;
			Dyvec[k +pn] = dFE6d0d00*u0d00dk + dFE6d0d01*u0d01dk + dFE6d0d02*u0d02dk + dFE6d1d00*u1d00dk + dFE6d1d01*u1d01dk + dFE6d1d02*u1d02dk + dFE6d2d00*u2d00dk + dFE6d2d01*u2d01dk + dFE6d2d02*u2d02dk + dFE6d3d00*u3d00dk + dFE6d3d01*u3d01dk + dFE6d3d02*u3d02dk + dFE6d4d00*u4d00dk + dFE6d4d01*u4d01dk + dFE6d4d02*u4d02dk + dFE6d5d00*u5d00dk + dFE6d5d01*u5d01dk + dFE6d5d02*u5d02dk + dFE6d6d00*u6d00dk + dFE6d6d01*u6d01dk + dFE6d6d02*u6d02dk + dFE6d7d00*u7d00dk + dFE6d7d01*u7d01dk + dFE6d7d02*u7d02dk + dFE6d8d00*u8d00dk + dFE6d8d01*u8d01dk + dFE6d8d02*u8d02dk + dFE6d9d00*u9d00dk + dFE6d9d01*u9d01dk + dFE6d9d02*u9d02dk;
			
			//FE7
			FE7       =       FE7out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d0d00 = dFE7d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d0d10 = dFE7d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d0d01 = dFE7d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d0d20 = dFE7d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d0d02 = dFE7d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d1d00 = dFE7d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d1d10 = dFE7d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d1d01 = dFE7d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d1d20 = dFE7d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d1d02 = dFE7d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d2d00 = dFE7d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d2d10 = dFE7d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d2d01 = dFE7d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d2d20 = dFE7d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d2d02 = dFE7d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d3d00 = dFE7d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d3d10 = dFE7d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d3d01 = dFE7d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d3d20 = dFE7d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d3d02 = dFE7d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d4d00 = dFE7d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d4d10 = dFE7d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d4d01 = dFE7d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d4d20 = dFE7d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d4d02 = dFE7d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d5d00 = dFE7d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d5d10 = dFE7d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d5d01 = dFE7d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d5d20 = dFE7d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d5d02 = dFE7d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d6d00 = dFE7d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d6d10 = dFE7d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d6d01 = dFE7d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d6d20 = dFE7d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d6d02 = dFE7d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d7d00 = dFE7d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d7d10 = dFE7d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d7d01 = dFE7d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d7d20 = dFE7d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d7d02 = dFE7d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d8d00 = dFE7d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d8d10 = dFE7d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d8d01 = dFE7d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d8d20 = dFE7d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d8d02 = dFE7d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE7d9d00 = dFE7d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d9d10 = dFE7d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d9d01 = dFE7d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d9d20 = dFE7d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE7d9d02 = dFE7d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 7;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE7d7d00*a7d0[a] + dFE7d7d10*a7d1[a] + dFE7d7d20*a7d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE7d2d01*b2d1[b] + dFE7d2d02*b2d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE7);
			Dxvec[k +pn] = dFE7d0d00*u0d00dk + dFE7d0d10*u0d10dk + dFE7d0d20*u0d20dk + dFE7d1d00*u1d00dk + dFE7d1d10*u1d10dk + dFE7d1d20*u1d20dk + dFE7d2d00*u2d00dk + dFE7d2d10*u2d10dk + dFE7d2d20*u2d20dk + dFE7d3d00*u3d00dk + dFE7d3d10*u3d10dk + dFE7d3d20*u3d20dk + dFE7d4d00*u4d00dk + dFE7d4d10*u4d10dk + dFE7d4d20*u4d20dk + dFE7d5d00*u5d00dk + dFE7d5d10*u5d10dk + dFE7d5d20*u5d20dk + dFE7d6d00*u6d00dk + dFE7d6d10*u6d10dk + dFE7d6d20*u6d20dk + dFE7d7d00*u7d00dk + dFE7d7d10*u7d10dk + dFE7d7d20*u7d20dk + dFE7d8d00*u8d00dk + dFE7d8d10*u8d10dk + dFE7d8d20*u8d20dk + dFE7d9d00*u9d00dk + dFE7d9d10*u9d10dk + dFE7d9d20*u9d20dk;
			Dyvec[k +pn] = dFE7d0d00*u0d00dk + dFE7d0d01*u0d01dk + dFE7d0d02*u0d02dk + dFE7d1d00*u1d00dk + dFE7d1d01*u1d01dk + dFE7d1d02*u1d02dk + dFE7d2d00*u2d00dk + dFE7d2d01*u2d01dk + dFE7d2d02*u2d02dk + dFE7d3d00*u3d00dk + dFE7d3d01*u3d01dk + dFE7d3d02*u3d02dk + dFE7d4d00*u4d00dk + dFE7d4d01*u4d01dk + dFE7d4d02*u4d02dk + dFE7d5d00*u5d00dk + dFE7d5d01*u5d01dk + dFE7d5d02*u5d02dk + dFE7d6d00*u6d00dk + dFE7d6d01*u6d01dk + dFE7d6d02*u6d02dk + dFE7d7d00*u7d00dk + dFE7d7d01*u7d01dk + dFE7d7d02*u7d02dk + dFE7d8d00*u8d00dk + dFE7d8d01*u8d01dk + dFE7d8d02*u8d02dk + dFE7d9d00*u9d00dk + dFE7d9d01*u9d01dk + dFE7d9d02*u9d02dk;
			
			//FE8
			FE8       =       FE8out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d0d00 = dFE8d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d0d10 = dFE8d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d0d01 = dFE8d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d0d20 = dFE8d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d0d02 = dFE8d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d1d00 = dFE8d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d1d10 = dFE8d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d1d01 = dFE8d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d1d20 = dFE8d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d1d02 = dFE8d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d2d00 = dFE8d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d2d10 = dFE8d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d2d01 = dFE8d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d2d20 = dFE8d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d2d02 = dFE8d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d3d00 = dFE8d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d3d10 = dFE8d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d3d01 = dFE8d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d3d20 = dFE8d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d3d02 = dFE8d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d4d00 = dFE8d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d4d10 = dFE8d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d4d01 = dFE8d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d4d20 = dFE8d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d4d02 = dFE8d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d5d00 = dFE8d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d5d10 = dFE8d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d5d01 = dFE8d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d5d20 = dFE8d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d5d02 = dFE8d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d6d00 = dFE8d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d6d10 = dFE8d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d6d01 = dFE8d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d6d20 = dFE8d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d6d02 = dFE8d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d7d00 = dFE8d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d7d10 = dFE8d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d7d01 = dFE8d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d7d20 = dFE8d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d7d02 = dFE8d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d8d00 = dFE8d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d8d10 = dFE8d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d8d01 = dFE8d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d8d20 = dFE8d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d8d02 = dFE8d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE8d9d00 = dFE8d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d9d10 = dFE8d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d9d01 = dFE8d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d9d20 = dFE8d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE8d9d02 = dFE8d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 8;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE8d8d00*a8d0[a] + dFE8d8d10*a8d1[a] + dFE8d8d20*a8d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE8d3d01*b3d1[b] + dFE8d3d02*b3d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE8);
			Dxvec[k +pn] = dFE8d0d00*u0d00dk + dFE8d0d10*u0d10dk + dFE8d0d20*u0d20dk + dFE8d1d00*u1d00dk + dFE8d1d10*u1d10dk + dFE8d1d20*u1d20dk + dFE8d2d00*u2d00dk + dFE8d2d10*u2d10dk + dFE8d2d20*u2d20dk + dFE8d3d00*u3d00dk + dFE8d3d10*u3d10dk + dFE8d3d20*u3d20dk + dFE8d4d00*u4d00dk + dFE8d4d10*u4d10dk + dFE8d4d20*u4d20dk + dFE8d5d00*u5d00dk + dFE8d5d10*u5d10dk + dFE8d5d20*u5d20dk + dFE8d6d00*u6d00dk + dFE8d6d10*u6d10dk + dFE8d6d20*u6d20dk + dFE8d7d00*u7d00dk + dFE8d7d10*u7d10dk + dFE8d7d20*u7d20dk + dFE8d8d00*u8d00dk + dFE8d8d10*u8d10dk + dFE8d8d20*u8d20dk + dFE8d9d00*u9d00dk + dFE8d9d10*u9d10dk + dFE8d9d20*u9d20dk;
			Dyvec[k +pn] = dFE8d0d00*u0d00dk + dFE8d0d01*u0d01dk + dFE8d0d02*u0d02dk + dFE8d1d00*u1d00dk + dFE8d1d01*u1d01dk + dFE8d1d02*u1d02dk + dFE8d2d00*u2d00dk + dFE8d2d01*u2d01dk + dFE8d2d02*u2d02dk + dFE8d3d00*u3d00dk + dFE8d3d01*u3d01dk + dFE8d3d02*u3d02dk + dFE8d4d00*u4d00dk + dFE8d4d01*u4d01dk + dFE8d4d02*u4d02dk + dFE8d5d00*u5d00dk + dFE8d5d01*u5d01dk + dFE8d5d02*u5d02dk + dFE8d6d00*u6d00dk + dFE8d6d01*u6d01dk + dFE8d6d02*u6d02dk + dFE8d7d00*u7d00dk + dFE8d7d01*u7d01dk + dFE8d7d02*u7d02dk + dFE8d8d00*u8d00dk + dFE8d8d01*u8d01dk + dFE8d8d02*u8d02dk + dFE8d9d00*u9d00dk + dFE8d9d01*u9d01dk + dFE8d9d02*u9d02dk;
			
			//FE9
			FE9       =       FE9out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d0d00 = dFE9d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d0d10 = dFE9d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d0d01 = dFE9d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d0d20 = dFE9d0d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d0d02 = dFE9d0d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d1d00 = dFE9d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d1d10 = dFE9d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d1d01 = dFE9d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d1d20 = dFE9d1d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d1d02 = dFE9d1d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d2d00 = dFE9d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d2d10 = dFE9d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d2d01 = dFE9d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d2d20 = dFE9d2d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d2d02 = dFE9d2d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d3d00 = dFE9d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d3d10 = dFE9d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d3d01 = dFE9d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d3d20 = dFE9d3d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d3d02 = dFE9d3d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d4d00 = dFE9d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d4d10 = dFE9d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d4d01 = dFE9d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d4d20 = dFE9d4d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d4d02 = dFE9d4d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d5d00 = dFE9d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d5d10 = dFE9d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d5d01 = dFE9d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d5d20 = dFE9d5d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d5d02 = dFE9d5d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d6d00 = dFE9d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d6d10 = dFE9d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d6d01 = dFE9d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d6d20 = dFE9d6d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d6d02 = dFE9d6d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d7d00 = dFE9d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d7d10 = dFE9d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d7d01 = dFE9d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d7d20 = dFE9d7d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d7d02 = dFE9d7d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d8d00 = dFE9d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d8d10 = dFE9d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d8d01 = dFE9d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d8d20 = dFE9d8d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d8d02 = dFE9d8d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			dFE9d9d00 = dFE9d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d9d10 = dFE9d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d9d01 = dFE9d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d9d20 = dFE9d9d20out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			dFE9d9d02 = dFE9d9d02out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u0d20k, u0d02k, u1d00k, u1d10k, u1d01k, u1d20k, u1d02k, u2d00k, u2d10k, u2d01k, u2d20k, u2d02k, u3d00k, u3d10k, u3d01k, u3d20k, u3d02k, u4d00k, u4d10k, u4d01k, u4d20k, u4d02k, u5d00k, u5d10k, u5d01k, u5d20k, u5d02k, u6d00k, u6d10k, u6d01k, u6d20k, u6d02k, u7d00k, u7d10k, u7d01k, u7d20k, u7d02k, u8d00k, u8d10k, u8d01k, u8d20k, u8d02k, u9d00k, u9d10k, u9d01k, u9d20k, u9d02k);
			
			//-------------------------------------------------------
			pn = 9;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dFE9d9d00*a9d0[a] + dFE9d9d10*a9d1[a] + dFE9d9d20*a9d2[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dFE9d4d01*b4d1[b] + dFE9d4d02*b4d2[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
				
			}
			
			bvec[k +pn]  = sgn*(FE9);
			Dxvec[k +pn] = dFE9d0d00*u0d00dk + dFE9d0d10*u0d10dk + dFE9d0d20*u0d20dk + dFE9d1d00*u1d00dk + dFE9d1d10*u1d10dk + dFE9d1d20*u1d20dk + dFE9d2d00*u2d00dk + dFE9d2d10*u2d10dk + dFE9d2d20*u2d20dk + dFE9d3d00*u3d00dk + dFE9d3d10*u3d10dk + dFE9d3d20*u3d20dk + dFE9d4d00*u4d00dk + dFE9d4d10*u4d10dk + dFE9d4d20*u4d20dk + dFE9d5d00*u5d00dk + dFE9d5d10*u5d10dk + dFE9d5d20*u5d20dk + dFE9d6d00*u6d00dk + dFE9d6d10*u6d10dk + dFE9d6d20*u6d20dk + dFE9d7d00*u7d00dk + dFE9d7d10*u7d10dk + dFE9d7d20*u7d20dk + dFE9d8d00*u8d00dk + dFE9d8d10*u8d10dk + dFE9d8d20*u8d20dk + dFE9d9d00*u9d00dk + dFE9d9d10*u9d10dk + dFE9d9d20*u9d20dk;
			Dyvec[k +pn] = dFE9d0d00*u0d00dk + dFE9d0d01*u0d01dk + dFE9d0d02*u0d02dk + dFE9d1d00*u1d00dk + dFE9d1d01*u1d01dk + dFE9d1d02*u1d02dk + dFE9d2d00*u2d00dk + dFE9d2d01*u2d01dk + dFE9d2d02*u2d02dk + dFE9d3d00*u3d00dk + dFE9d3d01*u3d01dk + dFE9d3d02*u3d02dk + dFE9d4d00*u4d00dk + dFE9d4d01*u4d01dk + dFE9d4d02*u4d02dk + dFE9d5d00*u5d00dk + dFE9d5d01*u5d01dk + dFE9d5d02*u5d02dk + dFE9d6d00*u6d00dk + dFE9d6d01*u6d01dk + dFE9d6d02*u6d02dk + dFE9d7d00*u7d00dk + dFE9d7d01*u7d01dk + dFE9d7d02*u7d02dk + dFE9d8d00*u8d00dk + dFE9d8d01*u8d01dk + dFE9d8d02*u8d02dk + dFE9d9d00*u9d00dk + dFE9d9d01*u9d01dk + dFE9d9d02*u9d02dk;
			
		}
	}
	
	//Cleanup
	free(sx);
	free(sy);
	
	free(a0d0);
	free(a0d1);
	free(a0d2);
	free(b0d0);
	free(b0d1);
	free(b0d2);
	
	free(a1d0);
	free(a1d1);
	free(a1d2);
	free(b1d0);
	free(b1d1);
	free(b1d2);
	
	free(a2d0);
	free(a2d1);
	free(a2d2);
	free(b2d0);
	free(b2d1);
	free(b2d2);
	
	free(a3d0);
	free(a3d1);
	free(a3d2);
	free(b3d0);
	free(b3d1);
	free(b3d2);
	
	free(a4d0);
	free(a4d1);
	free(a4d2);
	free(b4d0);
	free(b4d1);
	free(b4d2);
	
	free(a5d0);
	free(a5d1);
	free(a5d2);
	free(b5d0);
	free(b5d1);
	free(b5d2);
	
	free(a6d0);
	free(a6d1);
	free(a6d2);
	free(b6d0);
	free(b6d1);
	free(b6d2);
	
	free(a7d0);
	free(a7d1);
	free(a7d2);
	free(b7d0);
	free(b7d1);
	free(b7d2);
	
	free(a8d0);
	free(a8d1);
	free(a8d2);
	free(b8d0);
	free(b8d1);
	free(b8d2);
	
	free(a9d0);
	free(a9d1);
	free(a9d2);
	free(b9d0);
	free(b9d1);
	free(b9d2);
	
	
	return nnztemp;
}

int BVP_interp(struct param_type *params, double x[], double y[], double Psi[], int nnew, int mnew, double xnew[], double ynew[], double Psinew[])
{
	int i,j,k,l, status;
	const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;
	const int nnzJac = (*params).nnzparam;
	const double r_H = (*params).r_Hparam;
	const double M = (*params).Mparam;
	const double alpha = (*params).alphaparam;
	const double beta = (*params).betaparam;
	const size_t slen = r+4;
	int pn;
	const double sgn = -1.0;
	int a, b; //a for x, b for y(theta)
	int I, onexside, twoxside;
	int J, oneyside, twoyside;
	double xk, yk;
	const int dimx = 1;
	const int dimy = n;
	int nnztemp;
	
	int *sx;
	sx = (int*)calloc(slen, sizeof(int));
	int *sy;
	sy = (int*)calloc(slen, sizeof(int));
	
	//----- Finite difference stencil -----
	double *a0d0, *a0d1, *a0d2;
	a0d0 = (double*)calloc(slen, sizeof(double));
	a0d1 = (double*)calloc(slen, sizeof(double));
	a0d2 = (double*)calloc(slen, sizeof(double));
	double *b0d0, *b0d1, *b0d2;
	b0d0 = (double*)calloc(slen, sizeof(double));
	b0d1 = (double*)calloc(slen, sizeof(double));
	b0d2 = (double*)calloc(slen, sizeof(double));
	
	double *a1d0, *a1d1, *a1d2;
	a1d0 = (double*)calloc(slen, sizeof(double));
	a1d1 = (double*)calloc(slen, sizeof(double));
	a1d2 = (double*)calloc(slen, sizeof(double));
	double *b1d0, *b1d1, *b1d2;
	b1d0 = (double*)calloc(slen, sizeof(double));
	b1d1 = (double*)calloc(slen, sizeof(double));
	b1d2 = (double*)calloc(slen, sizeof(double));
	
	double *a2d0, *a2d1, *a2d2;
	a2d0 = (double*)calloc(slen, sizeof(double));
	a2d1 = (double*)calloc(slen, sizeof(double));
	a2d2 = (double*)calloc(slen, sizeof(double));
	double *b2d0, *b2d1, *b2d2;
	b2d0 = (double*)calloc(slen, sizeof(double));
	b2d1 = (double*)calloc(slen, sizeof(double));
	b2d2 = (double*)calloc(slen, sizeof(double));
	
	double *a3d0, *a3d1, *a3d2;
	a3d0 = (double*)calloc(slen, sizeof(double));
	a3d1 = (double*)calloc(slen, sizeof(double));
	a3d2 = (double*)calloc(slen, sizeof(double));
	double *b3d0, *b3d1, *b3d2;
	b3d0 = (double*)calloc(slen, sizeof(double));
	b3d1 = (double*)calloc(slen, sizeof(double));
	b3d2 = (double*)calloc(slen, sizeof(double));
	
	double *a4d0, *a4d1, *a4d2;
	a4d0 = (double*)calloc(slen, sizeof(double));
	a4d1 = (double*)calloc(slen, sizeof(double));
	a4d2 = (double*)calloc(slen, sizeof(double));
	double *b4d0, *b4d1, *b4d2;
	b4d0 = (double*)calloc(slen, sizeof(double));
	b4d1 = (double*)calloc(slen, sizeof(double));
	b4d2 = (double*)calloc(slen, sizeof(double));
	
	double *a5d0, *a5d1, *a5d2;
	a5d0 = (double*)calloc(slen, sizeof(double));
	a5d1 = (double*)calloc(slen, sizeof(double));
	a5d2 = (double*)calloc(slen, sizeof(double));
	double *b5d0, *b5d1, *b5d2;
	b5d0 = (double*)calloc(slen, sizeof(double));
	b5d1 = (double*)calloc(slen, sizeof(double));
	b5d2 = (double*)calloc(slen, sizeof(double));
	
	double *a6d0, *a6d1, *a6d2;
	a6d0 = (double*)calloc(slen, sizeof(double));
	a6d1 = (double*)calloc(slen, sizeof(double));
	a6d2 = (double*)calloc(slen, sizeof(double));
	double *b6d0, *b6d1, *b6d2;
	b6d0 = (double*)calloc(slen, sizeof(double));
	b6d1 = (double*)calloc(slen, sizeof(double));
	b6d2 = (double*)calloc(slen, sizeof(double));
	
	double *a7d0, *a7d1, *a7d2;
	a7d0 = (double*)calloc(slen, sizeof(double));
	a7d1 = (double*)calloc(slen, sizeof(double));
	a7d2 = (double*)calloc(slen, sizeof(double));
	double *b7d0, *b7d1, *b7d2;
	b7d0 = (double*)calloc(slen, sizeof(double));
	b7d1 = (double*)calloc(slen, sizeof(double));
	b7d2 = (double*)calloc(slen, sizeof(double));
	
	double *a8d0, *a8d1, *a8d2;
	a8d0 = (double*)calloc(slen, sizeof(double));
	a8d1 = (double*)calloc(slen, sizeof(double));
	a8d2 = (double*)calloc(slen, sizeof(double));
	double *b8d0, *b8d1, *b8d2;
	b8d0 = (double*)calloc(slen, sizeof(double));
	b8d1 = (double*)calloc(slen, sizeof(double));
	b8d2 = (double*)calloc(slen, sizeof(double));
	
	double *a9d0, *a9d1, *a9d2;
	a9d0 = (double*)calloc(slen, sizeof(double));
	a9d1 = (double*)calloc(slen, sizeof(double));
	a9d2 = (double*)calloc(slen, sizeof(double));
	double *b9d0, *b9d1, *b9d2;
	b9d0 = (double*)calloc(slen, sizeof(double));
	b9d1 = (double*)calloc(slen, sizeof(double));
	b9d2 = (double*)calloc(slen, sizeof(double));
	
	//----- Field derivatives initialization -----
	double u0d00k,  u0d10k,  u0d20k,  u0d01k,  u0d02k;
	double u0d00dk, u0d10dk, u0d20dk, u0d01dk, u0d02dk;
	
	double u1d00k,  u1d10k,  u1d20k,  u1d01k,  u1d02k;
	double u1d00dk, u1d10dk, u1d20dk, u1d01dk, u1d02dk;
	
	double u2d00k,  u2d10k,  u2d20k,  u2d01k,  u2d02k;
	double u2d00dk, u2d10dk, u2d20dk, u2d01dk, u2d02dk;
	
	double u3d00k,  u3d10k,  u3d20k,  u3d01k,  u3d02k;
	double u3d00dk, u3d10dk, u3d20dk, u3d01dk, u3d02dk;
	
	double u4d00k,  u4d10k,  u4d20k,  u4d01k,  u4d02k;
	double u4d00dk, u4d10dk, u4d20dk, u4d01dk, u4d02dk;
	
	double u5d00k,  u5d10k,  u5d20k,  u5d01k,  u5d02k;
	double u5d00dk, u5d10dk, u5d20dk, u5d01dk, u5d02dk;
	
	double u6d00k,  u6d10k,  u6d20k,  u6d01k,  u6d02k;
	double u6d00dk, u6d10dk, u6d20dk, u6d01dk, u6d02dk;
	
	double u7d00k,  u7d10k,  u7d20k,  u7d01k,  u7d02k;
	double u7d00dk, u7d10dk, u7d20dk, u7d01dk, u7d02dk;
	
	double u8d00k,  u8d10k,  u8d20k,  u8d01k,  u8d02k;
	double u8d00dk, u8d10dk, u8d20dk, u8d01dk, u8d02dk;
	
	double u9d00k,  u9d10k,  u9d20k,  u9d01k,  u9d02k;
	double u9d00dk, u9d10dk, u9d20dk, u9d01dk, u9d02dk;
	
	//:::::::::::::::::::: Main Loop ::::::::::::::::::::
	for (i = 0; i < mnew; ++i)
	{
		for (j = 0; j < nnew; ++j)
		{
			k = (i*nnew+j)*p;
			xk = xnew[j];
			yk = ynew[i];
			
			I = 0; //Determine I location in old y grid
			for (l = 0; l < m; ++l)
			{
				if (y[l] >= yk)
				{
					I = l;
					break;
				}
			}
			J = 0; //Determine J location in old x grid
			for (l = 0; l < n; ++l)
			{
				if (x[l] >= xk)
				{
					J = l;
					break;
				}
			}
			
			//Zero derivatives of each field
			onexside = 0; twoxside = 0;
			oneyside = 0; twoyside = 0;
			u0d00k  = 0.0, u0d10k  = 0.0, u0d20k  = 0.0, u0d01k  = 0.0, u0d02k  = 0.0;
			u0d00dk = 0.0, u0d10dk = 0.0, u0d20dk = 0.0, u0d01dk = 0.0, u0d02dk = 0.0;
			
			u1d00k  = 0.0, u1d10k  = 0.0, u1d20k  = 0.0, u1d01k  = 0.0, u1d02k  = 0.0;
			u1d00dk = 0.0, u1d10dk = 0.0, u1d20dk = 0.0, u1d01dk = 0.0, u1d02dk = 0.0;
			
			u2d00k  = 0.0, u2d10k  = 0.0, u2d20k  = 0.0, u2d01k  = 0.0, u2d02k  = 0.0;
			u2d00dk = 0.0, u2d10dk = 0.0, u2d20dk = 0.0, u2d01dk = 0.0, u2d02dk = 0.0;
			
			u3d00k  = 0.0, u3d10k  = 0.0, u3d20k  = 0.0, u3d01k  = 0.0, u3d02k  = 0.0;
			u3d00dk = 0.0, u3d10dk = 0.0, u3d20dk = 0.0, u3d01dk = 0.0, u3d02dk = 0.0;
			
			u4d00k  = 0.0, u4d10k  = 0.0, u4d20k  = 0.0, u4d01k  = 0.0, u4d02k  = 0.0;
			u4d00dk = 0.0, u4d10dk = 0.0, u4d20dk = 0.0, u4d01dk = 0.0, u4d02dk = 0.0;
			
			u5d00k  = 0.0, u5d10k  = 0.0, u5d20k  = 0.0, u5d01k  = 0.0, u5d02k  = 0.0;
			u5d00dk = 0.0, u5d10dk = 0.0, u5d20dk = 0.0, u5d01dk = 0.0, u5d02dk = 0.0;
			
			u6d00k  = 0.0, u6d10k  = 0.0, u6d20k  = 0.0, u6d01k  = 0.0, u6d02k  = 0.0;
			u6d00dk = 0.0, u6d10dk = 0.0, u6d20dk = 0.0, u6d01dk = 0.0, u6d02dk = 0.0;
			
			u7d00k  = 0.0, u7d10k  = 0.0, u7d20k  = 0.0, u7d01k  = 0.0, u7d02k  = 0.0;
			u7d00dk = 0.0, u7d10dk = 0.0, u7d20dk = 0.0, u7d01dk = 0.0, u7d02dk = 0.0;
			
			u8d00k  = 0.0, u8d10k  = 0.0, u8d20k  = 0.0, u8d01k  = 0.0, u8d02k  = 0.0;
			u8d00dk = 0.0, u8d10dk = 0.0, u8d20dk = 0.0, u8d01dk = 0.0, u8d02dk = 0.0;
			
			u9d00k  = 0.0, u9d10k  = 0.0, u9d20k  = 0.0, u9d01k  = 0.0, u9d02k  = 0.0;
			u9d00dk = 0.0, u9d10dk = 0.0, u9d20dk = 0.0, u9d01dk = 0.0, u9d02dk = 0.0;
			
			//Compute stencil
			status = steninit(sy, &oneyside, &twoyside, I, r, m);
			status = steninit(sx, &onexside, &twoxside, J, r, n);
			
			//Compute Newton polynomial representation for each field
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+0, Psi, 0, b0d0, b0d1, b0d2, &u0d00dk, &u0d01dk, &u0d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+0, Psi, 0, a0d0, a0d1, a0d2, &u0d00dk, &u0d10dk, &u0d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+1, Psi, 1, b1d0, b1d1, b1d2, &u1d00dk, &u1d01dk, &u1d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+1, Psi, 1, a1d0, a1d1, a1d2, &u1d00dk, &u1d10dk, &u1d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+2, Psi, 2, b2d0, b2d1, b2d2, &u2d00dk, &u2d01dk, &u2d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+2, Psi, 2, a2d0, a2d1, a2d2, &u2d00dk, &u2d10dk, &u2d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+3, Psi, 3, b3d0, b3d1, b3d2, &u3d00dk, &u3d01dk, &u3d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+3, Psi, 3, a3d0, a3d1, a3d2, &u3d00dk, &u3d10dk, &u3d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+4, Psi, 4, b4d0, b4d1, b4d2, &u4d00dk, &u4d01dk, &u4d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+4, Psi, 4, a4d0, a4d1, a4d2, &u4d00dk, &u4d10dk, &u4d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+5, Psi, 5, b5d0, b5d1, b5d2, &u5d00dk, &u5d01dk, &u5d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+5, Psi, 5, a5d0, a5d1, a5d2, &u5d00dk, &u5d10dk, &u5d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+6, Psi, 6, b6d0, b6d1, b6d2, &u6d00dk, &u6d01dk, &u6d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+6, Psi, 6, a6d0, a6d1, a6d2, &u6d00dk, &u6d10dk, &u6d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+7, Psi, 7, b7d0, b7d1, b7d2, &u7d00dk, &u7d01dk, &u7d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+7, Psi, 7, a7d0, a7d1, a7d2, &u7d00dk, &u7d10dk, &u7d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+8, Psi, 8, b8d0, b8d1, b8d2, &u8d00dk, &u8d01dk, &u8d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+8, Psi, 8, a8d0, a8d1, a8d2, &u8d00dk, &u8d10dk, &u8d20dk, xk);
			
			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+9, Psi, 9, b9d0, b9d1, b9d2, &u9d00dk, &u9d01dk, &u9d02dk, yk);
			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+9, Psi, 9, a9d0, a9d1, a9d2, &u9d00dk, &u9d10dk, &u9d20dk, xk);
			
			//Calculate derivative of each field
			for (b = 0; b <= r+oneyside; ++b)
			{
				u0d01k += Psi[((I+sy[b])*n +J)*p +0]*b0d1[b];
				u0d02k += Psi[((I+sy[b])*n +J)*p +0]*b0d2[b];
				
				u1d01k += Psi[((I+sy[b])*n +J)*p +1]*b1d1[b];
				u1d02k += Psi[((I+sy[b])*n +J)*p +1]*b1d2[b];
				
				u2d01k += Psi[((I+sy[b])*n +J)*p +2]*b2d1[b];
				u2d02k += Psi[((I+sy[b])*n +J)*p +2]*b2d2[b];
				
				u3d01k += Psi[((I+sy[b])*n +J)*p +3]*b3d1[b];
				u3d02k += Psi[((I+sy[b])*n +J)*p +3]*b3d2[b];
				
				u4d01k += Psi[((I+sy[b])*n +J)*p +4]*b4d1[b];
				u4d02k += Psi[((I+sy[b])*n +J)*p +4]*b4d2[b];
				
				u5d01k += Psi[((I+sy[b])*n +J)*p +5]*b5d1[b];
				u5d02k += Psi[((I+sy[b])*n +J)*p +5]*b5d2[b];
				
				u6d01k += Psi[((I+sy[b])*n +J)*p +6]*b6d1[b];
				u6d02k += Psi[((I+sy[b])*n +J)*p +6]*b6d2[b];
				
				u7d01k += Psi[((I+sy[b])*n +J)*p +7]*b7d1[b];
				u7d02k += Psi[((I+sy[b])*n +J)*p +7]*b7d2[b];
				
				u8d01k += Psi[((I+sy[b])*n +J)*p +8]*b8d1[b];
				u8d02k += Psi[((I+sy[b])*n +J)*p +8]*b8d2[b];
				
				u9d01k += Psi[((I+sy[b])*n +J)*p +9]*b9d1[b];
				u9d02k += Psi[((I+sy[b])*n +J)*p +9]*b9d2[b];
				
			}
			for (a = 0; a <= r+onexside; ++a)
			{
				u0d10k += Psi[(I*n +J+sx[a])*p +0]*a0d1[a];
				u0d20k += Psi[(I*n +J+sx[a])*p +0]*a0d2[a];
				
				u1d10k += Psi[(I*n +J+sx[a])*p +1]*a1d1[a];
				u1d20k += Psi[(I*n +J+sx[a])*p +1]*a1d2[a];
				
				u2d10k += Psi[(I*n +J+sx[a])*p +2]*a2d1[a];
				u2d20k += Psi[(I*n +J+sx[a])*p +2]*a2d2[a];
				
				u3d10k += Psi[(I*n +J+sx[a])*p +3]*a3d1[a];
				u3d20k += Psi[(I*n +J+sx[a])*p +3]*a3d2[a];
				
				u4d10k += Psi[(I*n +J+sx[a])*p +4]*a4d1[a];
				u4d20k += Psi[(I*n +J+sx[a])*p +4]*a4d2[a];
				
				u5d10k += Psi[(I*n +J+sx[a])*p +5]*a5d1[a];
				u5d20k += Psi[(I*n +J+sx[a])*p +5]*a5d2[a];
				
				u6d10k += Psi[(I*n +J+sx[a])*p +6]*a6d1[a];
				u6d20k += Psi[(I*n +J+sx[a])*p +6]*a6d2[a];
				
				u7d10k += Psi[(I*n +J+sx[a])*p +7]*a7d1[a];
				u7d20k += Psi[(I*n +J+sx[a])*p +7]*a7d2[a];
				
				u8d10k += Psi[(I*n +J+sx[a])*p +8]*a8d1[a];
				u8d20k += Psi[(I*n +J+sx[a])*p +8]*a8d2[a];
				
				u9d10k += Psi[(I*n +J+sx[a])*p +9]*a9d1[a];
				u9d20k += Psi[(I*n +J+sx[a])*p +9]*a9d2[a];
				
			}
			
			if (nnew == n)
			{
				for (b = 0; b <= r+oneyside; ++b)
				{
					u0d00k += Psi[((I+sy[b])*n +J)*p +0]*b0d0[b];
					u1d00k += Psi[((I+sy[b])*n +J)*p +1]*b1d0[b];
					u2d00k += Psi[((I+sy[b])*n +J)*p +2]*b2d0[b];
					u3d00k += Psi[((I+sy[b])*n +J)*p +3]*b3d0[b];
					u4d00k += Psi[((I+sy[b])*n +J)*p +4]*b4d0[b];
					u5d00k += Psi[((I+sy[b])*n +J)*p +5]*b5d0[b];
					u6d00k += Psi[((I+sy[b])*n +J)*p +6]*b6d0[b];
					u7d00k += Psi[((I+sy[b])*n +J)*p +7]*b7d0[b];
					u8d00k += Psi[((I+sy[b])*n +J)*p +8]*b8d0[b];
					u9d00k += Psi[((I+sy[b])*n +J)*p +9]*b9d0[b];
				}
			}
			else if (mnew == m)
			{
				for (a = 0; a <= r+onexside; ++a)
				{
					u0d00k += Psi[(I*n +J+sx[a])*p +0]*a0d0[a];
					u1d00k += Psi[(I*n +J+sx[a])*p +1]*a1d0[a];
					u2d00k += Psi[(I*n +J+sx[a])*p +2]*a2d0[a];
					u3d00k += Psi[(I*n +J+sx[a])*p +3]*a3d0[a];
					u4d00k += Psi[(I*n +J+sx[a])*p +4]*a4d0[a];
					u5d00k += Psi[(I*n +J+sx[a])*p +5]*a5d0[a];
					u6d00k += Psi[(I*n +J+sx[a])*p +6]*a6d0[a];
					u7d00k += Psi[(I*n +J+sx[a])*p +7]*a7d0[a];
					u8d00k += Psi[(I*n +J+sx[a])*p +8]*a8d0[a];
					u9d00k += Psi[(I*n +J+sx[a])*p +9]*a9d0[a];
				}
			}
			else
			{
				printf("ERROR: Interp cannot interpolate x and y simultaneously.\n");
				return -1;
			}
			
			//FE0
			pn = 0;
			Psinew[k +pn] = u0d00k;
			
			//FE1
			pn = 1;
			Psinew[k +pn] = u1d00k;
			
			//FE2
			pn = 2;
			Psinew[k +pn] = u2d00k;
			
			//FE3
			pn = 3;
			Psinew[k +pn] = u3d00k;
			
			//FE4
			pn = 4;
			Psinew[k +pn] = u4d00k;
			
			//FE5
			pn = 5;
			Psinew[k +pn] = u5d00k;
			
			//FE6
			pn = 6;
			Psinew[k +pn] = u6d00k;
			
			//FE7
			pn = 7;
			Psinew[k +pn] = u7d00k;
			
			//FE8
			pn = 8;
			Psinew[k +pn] = u8d00k;
			
			//FE9
			pn = 9;
			Psinew[k +pn] = u9d00k;
			
		}
	}
	
	//Cleanup
	free(sx);
	free(sy);
	
	free(a0d0);
	free(a0d1);
	free(a0d2);
	free(b0d0);
	free(b0d1);
	free(b0d2);
	
	free(a1d0);
	free(a1d1);
	free(a1d2);
	free(b1d0);
	free(b1d1);
	free(b1d2);
	
	free(a2d0);
	free(a2d1);
	free(a2d2);
	free(b2d0);
	free(b2d1);
	free(b2d2);
	
	free(a3d0);
	free(a3d1);
	free(a3d2);
	free(b3d0);
	free(b3d1);
	free(b3d2);
	
	free(a4d0);
	free(a4d1);
	free(a4d2);
	free(b4d0);
	free(b4d1);
	free(b4d2);
	
	free(a5d0);
	free(a5d1);
	free(a5d2);
	free(b5d0);
	free(b5d1);
	free(b5d2);
	
	free(a6d0);
	free(a6d1);
	free(a6d2);
	free(b6d0);
	free(b6d1);
	free(b6d2);
	
	free(a7d0);
	free(a7d1);
	free(a7d2);
	free(b7d0);
	free(b7d1);
	free(b7d2);
	
	free(a8d0);
	free(a8d1);
	free(a8d2);
	free(b8d0);
	free(b8d1);
	free(b8d2);
	
	free(a9d0);
	free(a9d1);
	free(a9d2);
	free(b9d0);
	free(b9d1);
	free(b9d2);
	
	
	return 0;
}



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::: Function Definitions :::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int BVP_NPcalc(const int p, const int r, const int rn, const int I, int s[], double x[], const int dim, const int IJ, double y[], const int pn, double Q0_a[], double Q1_a[], double Q2_a[], double *d0ydk, double *d1ydk, double *d2ydk, const double xk)
{
	int k,t,q,l,a,b,c;
	double s_ba;
	double P0_b;
	double P1_b;
	double P2_b;
	
	double temp;
	
	for (a = 0; a <= rn; ++a) //[2] loop
	{
		Q0_a[a] = 0.0;
		Q1_a[a] = 0.0;
		Q2_a[a] = 0.0;
	
		for (b = a; b <= rn; ++b)
		{
			s_ba = 1.0;
			for (c = 0; c <= b; ++c) //[3] loop
			{
				if (c != a)
				{
					s_ba /= (x[I+s[a]]-x[I+s[c]]); //[3]
				}
			}
			
			//P
			P0_b  = 0.0;
			temp = 1.0;
			for (c = 0; c <= b-1; ++c) //[4] loop
			{
				temp *= (xk - x[I+s[c]]); //[4]
			}
			P0_b += temp;
			
			//P'
			P1_b = 0.0;
			for (q = 0; q <= b-1; ++q) //[6] loop
			{
				temp = 1.0;
				for (c = 0; c <= b-1; ++c)
				{
					if (c != q)
					{
						temp *= (xk - x[I+s[c]]); //[6]
					}
				}
				P1_b += temp; //[5]
			}
			
			//P''
			P2_b = 0.0;
			for (t = 0; t <= b-1; ++t) //[7] loop
			{
				for (q = 0; q <= b-1; ++q)
				{
					if (q != t)
					{
						temp = 1.0;
						for (c = 0; c <= b-1; ++c)
						{
							if ((c != q) && (c != t))
							{
								temp *= (xk - x[I+s[c]]); //[7]
							}
						}
						P2_b += temp; //[8]
					}
				}
			}
			
			if (b > r)
			{
				*d0ydk += y[IJ + s[a]*p*dim]*s_ba*P0_b; //[1]
				*d1ydk += y[IJ + s[a]*p*dim]*s_ba*P1_b; //[1]
				*d2ydk += y[IJ + s[a]*p*dim]*s_ba*P2_b; //[1]
			}
			else
			{
				Q0_a[a] += s_ba*P0_b; //[2]
				Q1_a[a] += s_ba*P1_b; //[2]
				Q2_a[a] += s_ba*P2_b; //[2]
			}
		}
	}
	
	
	return 0;
}

int steninit(int s[], int *oneside, int *twoside, const int k, const int r, const int N)
{
	int l;
	s[0] = 0; //[1]
	s[(r+2)+1] = 0; //[1]
	for (l = 0; l < (r+2)/2; ++l) //[2]
	{
		s[2*l+1] = l+1;
		s[2*l+2] = -(l+1);
	}
	if (k < r/2) //[3]
	{
		*twoside = 0;
		*oneside = 1;
		for (l = 2*k+2; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] + 1;
		}
	}
	else if (k < (r+2)/2) //[4]
	{
		*twoside = 1;
		*oneside = 0;
		for (l = 2*k+2; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] + 1;
		}
	}
	else if ((N-1) - k < r/2) //[5]
	{
		*twoside = 0;
		*oneside = 1;
		for (l = 2*((N-1)-k)+1; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] - 1;
		}
	}
	else if ((N-1) - k < (r+2)/2) //[6]
	{
		*twoside = 1;
		*oneside = 0;
		for (l = 2*((N-1)-k)+1; l <= (r+2)+1; ++l)
		{
			s[l] = s[l-1] - 1;
		}
	}
	/*
	printf("k = %03i. oneside = %i. twoside = %i. s[] = ",k,*oneside,*twoside); //[7]
	for (l = 0; l < (r+2)+1; ++l)
	{
		printf(" %2i,",s[l]);
	}
	printf(" %2i.\n",s[(r+2)+1]);
	*/
	
	
	return 0;
}




