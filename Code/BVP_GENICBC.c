#include "BVP_header.h"

int BVP_GENBC(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[])
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
	
	//----- Boundary partials -----
	//BC X0
	double BC0X0;
	double dBC0X0d0d00, dBC0X0d0d10, dBC0X0d0d01;
	double dBC0X0d1d00, dBC0X0d1d10, dBC0X0d1d01;
	double dBC0X0d2d00, dBC0X0d2d10, dBC0X0d2d01;
	double dBC0X0d3d00, dBC0X0d3d10, dBC0X0d3d01;
	double dBC0X0d4d00, dBC0X0d4d10, dBC0X0d4d01;
	double dBC0X0d5d00, dBC0X0d5d10, dBC0X0d5d01;
	double dBC0X0d6d00, dBC0X0d6d10, dBC0X0d6d01;
	double dBC0X0d7d00, dBC0X0d7d10, dBC0X0d7d01;
	double dBC0X0d8d00, dBC0X0d8d10, dBC0X0d8d01;
	double dBC0X0d9d00, dBC0X0d9d10, dBC0X0d9d01;
	
	double BC1X0;
	double dBC1X0d0d00, dBC1X0d0d10, dBC1X0d0d01;
	double dBC1X0d1d00, dBC1X0d1d10, dBC1X0d1d01;
	double dBC1X0d2d00, dBC1X0d2d10, dBC1X0d2d01;
	double dBC1X0d3d00, dBC1X0d3d10, dBC1X0d3d01;
	double dBC1X0d4d00, dBC1X0d4d10, dBC1X0d4d01;
	double dBC1X0d5d00, dBC1X0d5d10, dBC1X0d5d01;
	double dBC1X0d6d00, dBC1X0d6d10, dBC1X0d6d01;
	double dBC1X0d7d00, dBC1X0d7d10, dBC1X0d7d01;
	double dBC1X0d8d00, dBC1X0d8d10, dBC1X0d8d01;
	double dBC1X0d9d00, dBC1X0d9d10, dBC1X0d9d01;
	
	double BC2X0;
	double dBC2X0d0d00, dBC2X0d0d10, dBC2X0d0d01;
	double dBC2X0d1d00, dBC2X0d1d10, dBC2X0d1d01;
	double dBC2X0d2d00, dBC2X0d2d10, dBC2X0d2d01;
	double dBC2X0d3d00, dBC2X0d3d10, dBC2X0d3d01;
	double dBC2X0d4d00, dBC2X0d4d10, dBC2X0d4d01;
	double dBC2X0d5d00, dBC2X0d5d10, dBC2X0d5d01;
	double dBC2X0d6d00, dBC2X0d6d10, dBC2X0d6d01;
	double dBC2X0d7d00, dBC2X0d7d10, dBC2X0d7d01;
	double dBC2X0d8d00, dBC2X0d8d10, dBC2X0d8d01;
	double dBC2X0d9d00, dBC2X0d9d10, dBC2X0d9d01;
	
	double BC3X0;
	double dBC3X0d0d00, dBC3X0d0d10, dBC3X0d0d01;
	double dBC3X0d1d00, dBC3X0d1d10, dBC3X0d1d01;
	double dBC3X0d2d00, dBC3X0d2d10, dBC3X0d2d01;
	double dBC3X0d3d00, dBC3X0d3d10, dBC3X0d3d01;
	double dBC3X0d4d00, dBC3X0d4d10, dBC3X0d4d01;
	double dBC3X0d5d00, dBC3X0d5d10, dBC3X0d5d01;
	double dBC3X0d6d00, dBC3X0d6d10, dBC3X0d6d01;
	double dBC3X0d7d00, dBC3X0d7d10, dBC3X0d7d01;
	double dBC3X0d8d00, dBC3X0d8d10, dBC3X0d8d01;
	double dBC3X0d9d00, dBC3X0d9d10, dBC3X0d9d01;
	
	double BC4X0;
	double dBC4X0d0d00, dBC4X0d0d10, dBC4X0d0d01;
	double dBC4X0d1d00, dBC4X0d1d10, dBC4X0d1d01;
	double dBC4X0d2d00, dBC4X0d2d10, dBC4X0d2d01;
	double dBC4X0d3d00, dBC4X0d3d10, dBC4X0d3d01;
	double dBC4X0d4d00, dBC4X0d4d10, dBC4X0d4d01;
	double dBC4X0d5d00, dBC4X0d5d10, dBC4X0d5d01;
	double dBC4X0d6d00, dBC4X0d6d10, dBC4X0d6d01;
	double dBC4X0d7d00, dBC4X0d7d10, dBC4X0d7d01;
	double dBC4X0d8d00, dBC4X0d8d10, dBC4X0d8d01;
	double dBC4X0d9d00, dBC4X0d9d10, dBC4X0d9d01;
	
	double BC5X0;
	double dBC5X0d0d00, dBC5X0d0d10, dBC5X0d0d01;
	double dBC5X0d1d00, dBC5X0d1d10, dBC5X0d1d01;
	double dBC5X0d2d00, dBC5X0d2d10, dBC5X0d2d01;
	double dBC5X0d3d00, dBC5X0d3d10, dBC5X0d3d01;
	double dBC5X0d4d00, dBC5X0d4d10, dBC5X0d4d01;
	double dBC5X0d5d00, dBC5X0d5d10, dBC5X0d5d01;
	double dBC5X0d6d00, dBC5X0d6d10, dBC5X0d6d01;
	double dBC5X0d7d00, dBC5X0d7d10, dBC5X0d7d01;
	double dBC5X0d8d00, dBC5X0d8d10, dBC5X0d8d01;
	double dBC5X0d9d00, dBC5X0d9d10, dBC5X0d9d01;
	
	double BC6X0;
	double dBC6X0d0d00, dBC6X0d0d10, dBC6X0d0d01;
	double dBC6X0d1d00, dBC6X0d1d10, dBC6X0d1d01;
	double dBC6X0d2d00, dBC6X0d2d10, dBC6X0d2d01;
	double dBC6X0d3d00, dBC6X0d3d10, dBC6X0d3d01;
	double dBC6X0d4d00, dBC6X0d4d10, dBC6X0d4d01;
	double dBC6X0d5d00, dBC6X0d5d10, dBC6X0d5d01;
	double dBC6X0d6d00, dBC6X0d6d10, dBC6X0d6d01;
	double dBC6X0d7d00, dBC6X0d7d10, dBC6X0d7d01;
	double dBC6X0d8d00, dBC6X0d8d10, dBC6X0d8d01;
	double dBC6X0d9d00, dBC6X0d9d10, dBC6X0d9d01;
	
	double BC7X0;
	double dBC7X0d0d00, dBC7X0d0d10, dBC7X0d0d01;
	double dBC7X0d1d00, dBC7X0d1d10, dBC7X0d1d01;
	double dBC7X0d2d00, dBC7X0d2d10, dBC7X0d2d01;
	double dBC7X0d3d00, dBC7X0d3d10, dBC7X0d3d01;
	double dBC7X0d4d00, dBC7X0d4d10, dBC7X0d4d01;
	double dBC7X0d5d00, dBC7X0d5d10, dBC7X0d5d01;
	double dBC7X0d6d00, dBC7X0d6d10, dBC7X0d6d01;
	double dBC7X0d7d00, dBC7X0d7d10, dBC7X0d7d01;
	double dBC7X0d8d00, dBC7X0d8d10, dBC7X0d8d01;
	double dBC7X0d9d00, dBC7X0d9d10, dBC7X0d9d01;
	
	double BC8X0;
	double dBC8X0d0d00, dBC8X0d0d10, dBC8X0d0d01;
	double dBC8X0d1d00, dBC8X0d1d10, dBC8X0d1d01;
	double dBC8X0d2d00, dBC8X0d2d10, dBC8X0d2d01;
	double dBC8X0d3d00, dBC8X0d3d10, dBC8X0d3d01;
	double dBC8X0d4d00, dBC8X0d4d10, dBC8X0d4d01;
	double dBC8X0d5d00, dBC8X0d5d10, dBC8X0d5d01;
	double dBC8X0d6d00, dBC8X0d6d10, dBC8X0d6d01;
	double dBC8X0d7d00, dBC8X0d7d10, dBC8X0d7d01;
	double dBC8X0d8d00, dBC8X0d8d10, dBC8X0d8d01;
	double dBC8X0d9d00, dBC8X0d9d10, dBC8X0d9d01;
	
	double BC9X0;
	double dBC9X0d0d00, dBC9X0d0d10, dBC9X0d0d01;
	double dBC9X0d1d00, dBC9X0d1d10, dBC9X0d1d01;
	double dBC9X0d2d00, dBC9X0d2d10, dBC9X0d2d01;
	double dBC9X0d3d00, dBC9X0d3d10, dBC9X0d3d01;
	double dBC9X0d4d00, dBC9X0d4d10, dBC9X0d4d01;
	double dBC9X0d5d00, dBC9X0d5d10, dBC9X0d5d01;
	double dBC9X0d6d00, dBC9X0d6d10, dBC9X0d6d01;
	double dBC9X0d7d00, dBC9X0d7d10, dBC9X0d7d01;
	double dBC9X0d8d00, dBC9X0d8d10, dBC9X0d8d01;
	double dBC9X0d9d00, dBC9X0d9d10, dBC9X0d9d01;
	
	//BC X1
	double BC0X1;
	double dBC0X1d0d00, dBC0X1d0d10, dBC0X1d0d01;
	double dBC0X1d1d00, dBC0X1d1d10, dBC0X1d1d01;
	double dBC0X1d2d00, dBC0X1d2d10, dBC0X1d2d01;
	double dBC0X1d3d00, dBC0X1d3d10, dBC0X1d3d01;
	double dBC0X1d4d00, dBC0X1d4d10, dBC0X1d4d01;
	double dBC0X1d5d00, dBC0X1d5d10, dBC0X1d5d01;
	double dBC0X1d6d00, dBC0X1d6d10, dBC0X1d6d01;
	double dBC0X1d7d00, dBC0X1d7d10, dBC0X1d7d01;
	double dBC0X1d8d00, dBC0X1d8d10, dBC0X1d8d01;
	double dBC0X1d9d00, dBC0X1d9d10, dBC0X1d9d01;
	
	double BC1X1;
	double dBC1X1d0d00, dBC1X1d0d10, dBC1X1d0d01;
	double dBC1X1d1d00, dBC1X1d1d10, dBC1X1d1d01;
	double dBC1X1d2d00, dBC1X1d2d10, dBC1X1d2d01;
	double dBC1X1d3d00, dBC1X1d3d10, dBC1X1d3d01;
	double dBC1X1d4d00, dBC1X1d4d10, dBC1X1d4d01;
	double dBC1X1d5d00, dBC1X1d5d10, dBC1X1d5d01;
	double dBC1X1d6d00, dBC1X1d6d10, dBC1X1d6d01;
	double dBC1X1d7d00, dBC1X1d7d10, dBC1X1d7d01;
	double dBC1X1d8d00, dBC1X1d8d10, dBC1X1d8d01;
	double dBC1X1d9d00, dBC1X1d9d10, dBC1X1d9d01;
	
	double BC2X1;
	double dBC2X1d0d00, dBC2X1d0d10, dBC2X1d0d01;
	double dBC2X1d1d00, dBC2X1d1d10, dBC2X1d1d01;
	double dBC2X1d2d00, dBC2X1d2d10, dBC2X1d2d01;
	double dBC2X1d3d00, dBC2X1d3d10, dBC2X1d3d01;
	double dBC2X1d4d00, dBC2X1d4d10, dBC2X1d4d01;
	double dBC2X1d5d00, dBC2X1d5d10, dBC2X1d5d01;
	double dBC2X1d6d00, dBC2X1d6d10, dBC2X1d6d01;
	double dBC2X1d7d00, dBC2X1d7d10, dBC2X1d7d01;
	double dBC2X1d8d00, dBC2X1d8d10, dBC2X1d8d01;
	double dBC2X1d9d00, dBC2X1d9d10, dBC2X1d9d01;
	
	double BC3X1;
	double dBC3X1d0d00, dBC3X1d0d10, dBC3X1d0d01;
	double dBC3X1d1d00, dBC3X1d1d10, dBC3X1d1d01;
	double dBC3X1d2d00, dBC3X1d2d10, dBC3X1d2d01;
	double dBC3X1d3d00, dBC3X1d3d10, dBC3X1d3d01;
	double dBC3X1d4d00, dBC3X1d4d10, dBC3X1d4d01;
	double dBC3X1d5d00, dBC3X1d5d10, dBC3X1d5d01;
	double dBC3X1d6d00, dBC3X1d6d10, dBC3X1d6d01;
	double dBC3X1d7d00, dBC3X1d7d10, dBC3X1d7d01;
	double dBC3X1d8d00, dBC3X1d8d10, dBC3X1d8d01;
	double dBC3X1d9d00, dBC3X1d9d10, dBC3X1d9d01;
	
	double BC4X1;
	double dBC4X1d0d00, dBC4X1d0d10, dBC4X1d0d01;
	double dBC4X1d1d00, dBC4X1d1d10, dBC4X1d1d01;
	double dBC4X1d2d00, dBC4X1d2d10, dBC4X1d2d01;
	double dBC4X1d3d00, dBC4X1d3d10, dBC4X1d3d01;
	double dBC4X1d4d00, dBC4X1d4d10, dBC4X1d4d01;
	double dBC4X1d5d00, dBC4X1d5d10, dBC4X1d5d01;
	double dBC4X1d6d00, dBC4X1d6d10, dBC4X1d6d01;
	double dBC4X1d7d00, dBC4X1d7d10, dBC4X1d7d01;
	double dBC4X1d8d00, dBC4X1d8d10, dBC4X1d8d01;
	double dBC4X1d9d00, dBC4X1d9d10, dBC4X1d9d01;
	
	double BC5X1;
	double dBC5X1d0d00, dBC5X1d0d10, dBC5X1d0d01;
	double dBC5X1d1d00, dBC5X1d1d10, dBC5X1d1d01;
	double dBC5X1d2d00, dBC5X1d2d10, dBC5X1d2d01;
	double dBC5X1d3d00, dBC5X1d3d10, dBC5X1d3d01;
	double dBC5X1d4d00, dBC5X1d4d10, dBC5X1d4d01;
	double dBC5X1d5d00, dBC5X1d5d10, dBC5X1d5d01;
	double dBC5X1d6d00, dBC5X1d6d10, dBC5X1d6d01;
	double dBC5X1d7d00, dBC5X1d7d10, dBC5X1d7d01;
	double dBC5X1d8d00, dBC5X1d8d10, dBC5X1d8d01;
	double dBC5X1d9d00, dBC5X1d9d10, dBC5X1d9d01;
	
	double BC6X1;
	double dBC6X1d0d00, dBC6X1d0d10, dBC6X1d0d01;
	double dBC6X1d1d00, dBC6X1d1d10, dBC6X1d1d01;
	double dBC6X1d2d00, dBC6X1d2d10, dBC6X1d2d01;
	double dBC6X1d3d00, dBC6X1d3d10, dBC6X1d3d01;
	double dBC6X1d4d00, dBC6X1d4d10, dBC6X1d4d01;
	double dBC6X1d5d00, dBC6X1d5d10, dBC6X1d5d01;
	double dBC6X1d6d00, dBC6X1d6d10, dBC6X1d6d01;
	double dBC6X1d7d00, dBC6X1d7d10, dBC6X1d7d01;
	double dBC6X1d8d00, dBC6X1d8d10, dBC6X1d8d01;
	double dBC6X1d9d00, dBC6X1d9d10, dBC6X1d9d01;
	
	double BC7X1;
	double dBC7X1d0d00, dBC7X1d0d10, dBC7X1d0d01;
	double dBC7X1d1d00, dBC7X1d1d10, dBC7X1d1d01;
	double dBC7X1d2d00, dBC7X1d2d10, dBC7X1d2d01;
	double dBC7X1d3d00, dBC7X1d3d10, dBC7X1d3d01;
	double dBC7X1d4d00, dBC7X1d4d10, dBC7X1d4d01;
	double dBC7X1d5d00, dBC7X1d5d10, dBC7X1d5d01;
	double dBC7X1d6d00, dBC7X1d6d10, dBC7X1d6d01;
	double dBC7X1d7d00, dBC7X1d7d10, dBC7X1d7d01;
	double dBC7X1d8d00, dBC7X1d8d10, dBC7X1d8d01;
	double dBC7X1d9d00, dBC7X1d9d10, dBC7X1d9d01;
	
	double BC8X1;
	double dBC8X1d0d00, dBC8X1d0d10, dBC8X1d0d01;
	double dBC8X1d1d00, dBC8X1d1d10, dBC8X1d1d01;
	double dBC8X1d2d00, dBC8X1d2d10, dBC8X1d2d01;
	double dBC8X1d3d00, dBC8X1d3d10, dBC8X1d3d01;
	double dBC8X1d4d00, dBC8X1d4d10, dBC8X1d4d01;
	double dBC8X1d5d00, dBC8X1d5d10, dBC8X1d5d01;
	double dBC8X1d6d00, dBC8X1d6d10, dBC8X1d6d01;
	double dBC8X1d7d00, dBC8X1d7d10, dBC8X1d7d01;
	double dBC8X1d8d00, dBC8X1d8d10, dBC8X1d8d01;
	double dBC8X1d9d00, dBC8X1d9d10, dBC8X1d9d01;
	
	double BC9X1;
	double dBC9X1d0d00, dBC9X1d0d10, dBC9X1d0d01;
	double dBC9X1d1d00, dBC9X1d1d10, dBC9X1d1d01;
	double dBC9X1d2d00, dBC9X1d2d10, dBC9X1d2d01;
	double dBC9X1d3d00, dBC9X1d3d10, dBC9X1d3d01;
	double dBC9X1d4d00, dBC9X1d4d10, dBC9X1d4d01;
	double dBC9X1d5d00, dBC9X1d5d10, dBC9X1d5d01;
	double dBC9X1d6d00, dBC9X1d6d10, dBC9X1d6d01;
	double dBC9X1d7d00, dBC9X1d7d10, dBC9X1d7d01;
	double dBC9X1d8d00, dBC9X1d8d10, dBC9X1d8d01;
	double dBC9X1d9d00, dBC9X1d9d10, dBC9X1d9d01;
	
	//BC Y1
	double BC0Y1;
	double dBC0Y1d0d00, dBC0Y1d0d10, dBC0Y1d0d01;
	double dBC0Y1d1d00, dBC0Y1d1d10, dBC0Y1d1d01;
	double dBC0Y1d2d00, dBC0Y1d2d10, dBC0Y1d2d01;
	double dBC0Y1d3d00, dBC0Y1d3d10, dBC0Y1d3d01;
	double dBC0Y1d4d00, dBC0Y1d4d10, dBC0Y1d4d01;
	double dBC0Y1d5d00, dBC0Y1d5d10, dBC0Y1d5d01;
	double dBC0Y1d6d00, dBC0Y1d6d10, dBC0Y1d6d01;
	double dBC0Y1d7d00, dBC0Y1d7d10, dBC0Y1d7d01;
	double dBC0Y1d8d00, dBC0Y1d8d10, dBC0Y1d8d01;
	double dBC0Y1d9d00, dBC0Y1d9d10, dBC0Y1d9d01;
	
	double BC1Y1;
	double dBC1Y1d0d00, dBC1Y1d0d10, dBC1Y1d0d01;
	double dBC1Y1d1d00, dBC1Y1d1d10, dBC1Y1d1d01;
	double dBC1Y1d2d00, dBC1Y1d2d10, dBC1Y1d2d01;
	double dBC1Y1d3d00, dBC1Y1d3d10, dBC1Y1d3d01;
	double dBC1Y1d4d00, dBC1Y1d4d10, dBC1Y1d4d01;
	double dBC1Y1d5d00, dBC1Y1d5d10, dBC1Y1d5d01;
	double dBC1Y1d6d00, dBC1Y1d6d10, dBC1Y1d6d01;
	double dBC1Y1d7d00, dBC1Y1d7d10, dBC1Y1d7d01;
	double dBC1Y1d8d00, dBC1Y1d8d10, dBC1Y1d8d01;
	double dBC1Y1d9d00, dBC1Y1d9d10, dBC1Y1d9d01;
	
	double BC2Y1;
	double dBC2Y1d0d00, dBC2Y1d0d10, dBC2Y1d0d01;
	double dBC2Y1d1d00, dBC2Y1d1d10, dBC2Y1d1d01;
	double dBC2Y1d2d00, dBC2Y1d2d10, dBC2Y1d2d01;
	double dBC2Y1d3d00, dBC2Y1d3d10, dBC2Y1d3d01;
	double dBC2Y1d4d00, dBC2Y1d4d10, dBC2Y1d4d01;
	double dBC2Y1d5d00, dBC2Y1d5d10, dBC2Y1d5d01;
	double dBC2Y1d6d00, dBC2Y1d6d10, dBC2Y1d6d01;
	double dBC2Y1d7d00, dBC2Y1d7d10, dBC2Y1d7d01;
	double dBC2Y1d8d00, dBC2Y1d8d10, dBC2Y1d8d01;
	double dBC2Y1d9d00, dBC2Y1d9d10, dBC2Y1d9d01;
	
	double BC3Y1;
	double dBC3Y1d0d00, dBC3Y1d0d10, dBC3Y1d0d01;
	double dBC3Y1d1d00, dBC3Y1d1d10, dBC3Y1d1d01;
	double dBC3Y1d2d00, dBC3Y1d2d10, dBC3Y1d2d01;
	double dBC3Y1d3d00, dBC3Y1d3d10, dBC3Y1d3d01;
	double dBC3Y1d4d00, dBC3Y1d4d10, dBC3Y1d4d01;
	double dBC3Y1d5d00, dBC3Y1d5d10, dBC3Y1d5d01;
	double dBC3Y1d6d00, dBC3Y1d6d10, dBC3Y1d6d01;
	double dBC3Y1d7d00, dBC3Y1d7d10, dBC3Y1d7d01;
	double dBC3Y1d8d00, dBC3Y1d8d10, dBC3Y1d8d01;
	double dBC3Y1d9d00, dBC3Y1d9d10, dBC3Y1d9d01;
	
	double BC4Y1;
	double dBC4Y1d0d00, dBC4Y1d0d10, dBC4Y1d0d01;
	double dBC4Y1d1d00, dBC4Y1d1d10, dBC4Y1d1d01;
	double dBC4Y1d2d00, dBC4Y1d2d10, dBC4Y1d2d01;
	double dBC4Y1d3d00, dBC4Y1d3d10, dBC4Y1d3d01;
	double dBC4Y1d4d00, dBC4Y1d4d10, dBC4Y1d4d01;
	double dBC4Y1d5d00, dBC4Y1d5d10, dBC4Y1d5d01;
	double dBC4Y1d6d00, dBC4Y1d6d10, dBC4Y1d6d01;
	double dBC4Y1d7d00, dBC4Y1d7d10, dBC4Y1d7d01;
	double dBC4Y1d8d00, dBC4Y1d8d10, dBC4Y1d8d01;
	double dBC4Y1d9d00, dBC4Y1d9d10, dBC4Y1d9d01;
	
	double BC5Y1;
	double dBC5Y1d0d00, dBC5Y1d0d10, dBC5Y1d0d01;
	double dBC5Y1d1d00, dBC5Y1d1d10, dBC5Y1d1d01;
	double dBC5Y1d2d00, dBC5Y1d2d10, dBC5Y1d2d01;
	double dBC5Y1d3d00, dBC5Y1d3d10, dBC5Y1d3d01;
	double dBC5Y1d4d00, dBC5Y1d4d10, dBC5Y1d4d01;
	double dBC5Y1d5d00, dBC5Y1d5d10, dBC5Y1d5d01;
	double dBC5Y1d6d00, dBC5Y1d6d10, dBC5Y1d6d01;
	double dBC5Y1d7d00, dBC5Y1d7d10, dBC5Y1d7d01;
	double dBC5Y1d8d00, dBC5Y1d8d10, dBC5Y1d8d01;
	double dBC5Y1d9d00, dBC5Y1d9d10, dBC5Y1d9d01;
	
	double BC6Y1;
	double dBC6Y1d0d00, dBC6Y1d0d10, dBC6Y1d0d01;
	double dBC6Y1d1d00, dBC6Y1d1d10, dBC6Y1d1d01;
	double dBC6Y1d2d00, dBC6Y1d2d10, dBC6Y1d2d01;
	double dBC6Y1d3d00, dBC6Y1d3d10, dBC6Y1d3d01;
	double dBC6Y1d4d00, dBC6Y1d4d10, dBC6Y1d4d01;
	double dBC6Y1d5d00, dBC6Y1d5d10, dBC6Y1d5d01;
	double dBC6Y1d6d00, dBC6Y1d6d10, dBC6Y1d6d01;
	double dBC6Y1d7d00, dBC6Y1d7d10, dBC6Y1d7d01;
	double dBC6Y1d8d00, dBC6Y1d8d10, dBC6Y1d8d01;
	double dBC6Y1d9d00, dBC6Y1d9d10, dBC6Y1d9d01;
	
	double BC7Y1;
	double dBC7Y1d0d00, dBC7Y1d0d10, dBC7Y1d0d01;
	double dBC7Y1d1d00, dBC7Y1d1d10, dBC7Y1d1d01;
	double dBC7Y1d2d00, dBC7Y1d2d10, dBC7Y1d2d01;
	double dBC7Y1d3d00, dBC7Y1d3d10, dBC7Y1d3d01;
	double dBC7Y1d4d00, dBC7Y1d4d10, dBC7Y1d4d01;
	double dBC7Y1d5d00, dBC7Y1d5d10, dBC7Y1d5d01;
	double dBC7Y1d6d00, dBC7Y1d6d10, dBC7Y1d6d01;
	double dBC7Y1d7d00, dBC7Y1d7d10, dBC7Y1d7d01;
	double dBC7Y1d8d00, dBC7Y1d8d10, dBC7Y1d8d01;
	double dBC7Y1d9d00, dBC7Y1d9d10, dBC7Y1d9d01;
	
	double BC8Y1;
	double dBC8Y1d0d00, dBC8Y1d0d10, dBC8Y1d0d01;
	double dBC8Y1d1d00, dBC8Y1d1d10, dBC8Y1d1d01;
	double dBC8Y1d2d00, dBC8Y1d2d10, dBC8Y1d2d01;
	double dBC8Y1d3d00, dBC8Y1d3d10, dBC8Y1d3d01;
	double dBC8Y1d4d00, dBC8Y1d4d10, dBC8Y1d4d01;
	double dBC8Y1d5d00, dBC8Y1d5d10, dBC8Y1d5d01;
	double dBC8Y1d6d00, dBC8Y1d6d10, dBC8Y1d6d01;
	double dBC8Y1d7d00, dBC8Y1d7d10, dBC8Y1d7d01;
	double dBC8Y1d8d00, dBC8Y1d8d10, dBC8Y1d8d01;
	double dBC8Y1d9d00, dBC8Y1d9d10, dBC8Y1d9d01;
	
	double BC9Y1;
	double dBC9Y1d0d00, dBC9Y1d0d10, dBC9Y1d0d01;
	double dBC9Y1d1d00, dBC9Y1d1d10, dBC9Y1d1d01;
	double dBC9Y1d2d00, dBC9Y1d2d10, dBC9Y1d2d01;
	double dBC9Y1d3d00, dBC9Y1d3d10, dBC9Y1d3d01;
	double dBC9Y1d4d00, dBC9Y1d4d10, dBC9Y1d4d01;
	double dBC9Y1d5d00, dBC9Y1d5d10, dBC9Y1d5d01;
	double dBC9Y1d6d00, dBC9Y1d6d10, dBC9Y1d6d01;
	double dBC9Y1d7d00, dBC9Y1d7d10, dBC9Y1d7d01;
	double dBC9Y1d8d00, dBC9Y1d8d10, dBC9Y1d8d01;
	double dBC9Y1d9d00, dBC9Y1d9d10, dBC9Y1d9d01;
	
	//BC Y0
	double BC0Y0;
	double dBC0Y0d0d00, dBC0Y0d0d10, dBC0Y0d0d01;
	double dBC0Y0d1d00, dBC0Y0d1d10, dBC0Y0d1d01;
	double dBC0Y0d2d00, dBC0Y0d2d10, dBC0Y0d2d01;
	double dBC0Y0d3d00, dBC0Y0d3d10, dBC0Y0d3d01;
	double dBC0Y0d4d00, dBC0Y0d4d10, dBC0Y0d4d01;
	double dBC0Y0d5d00, dBC0Y0d5d10, dBC0Y0d5d01;
	double dBC0Y0d6d00, dBC0Y0d6d10, dBC0Y0d6d01;
	double dBC0Y0d7d00, dBC0Y0d7d10, dBC0Y0d7d01;
	double dBC0Y0d8d00, dBC0Y0d8d10, dBC0Y0d8d01;
	double dBC0Y0d9d00, dBC0Y0d9d10, dBC0Y0d9d01;
	
	double BC1Y0;
	double dBC1Y0d0d00, dBC1Y0d0d10, dBC1Y0d0d01;
	double dBC1Y0d1d00, dBC1Y0d1d10, dBC1Y0d1d01;
	double dBC1Y0d2d00, dBC1Y0d2d10, dBC1Y0d2d01;
	double dBC1Y0d3d00, dBC1Y0d3d10, dBC1Y0d3d01;
	double dBC1Y0d4d00, dBC1Y0d4d10, dBC1Y0d4d01;
	double dBC1Y0d5d00, dBC1Y0d5d10, dBC1Y0d5d01;
	double dBC1Y0d6d00, dBC1Y0d6d10, dBC1Y0d6d01;
	double dBC1Y0d7d00, dBC1Y0d7d10, dBC1Y0d7d01;
	double dBC1Y0d8d00, dBC1Y0d8d10, dBC1Y0d8d01;
	double dBC1Y0d9d00, dBC1Y0d9d10, dBC1Y0d9d01;
	
	double BC2Y0;
	double dBC2Y0d0d00, dBC2Y0d0d10, dBC2Y0d0d01;
	double dBC2Y0d1d00, dBC2Y0d1d10, dBC2Y0d1d01;
	double dBC2Y0d2d00, dBC2Y0d2d10, dBC2Y0d2d01;
	double dBC2Y0d3d00, dBC2Y0d3d10, dBC2Y0d3d01;
	double dBC2Y0d4d00, dBC2Y0d4d10, dBC2Y0d4d01;
	double dBC2Y0d5d00, dBC2Y0d5d10, dBC2Y0d5d01;
	double dBC2Y0d6d00, dBC2Y0d6d10, dBC2Y0d6d01;
	double dBC2Y0d7d00, dBC2Y0d7d10, dBC2Y0d7d01;
	double dBC2Y0d8d00, dBC2Y0d8d10, dBC2Y0d8d01;
	double dBC2Y0d9d00, dBC2Y0d9d10, dBC2Y0d9d01;
	
	double BC3Y0;
	double dBC3Y0d0d00, dBC3Y0d0d10, dBC3Y0d0d01;
	double dBC3Y0d1d00, dBC3Y0d1d10, dBC3Y0d1d01;
	double dBC3Y0d2d00, dBC3Y0d2d10, dBC3Y0d2d01;
	double dBC3Y0d3d00, dBC3Y0d3d10, dBC3Y0d3d01;
	double dBC3Y0d4d00, dBC3Y0d4d10, dBC3Y0d4d01;
	double dBC3Y0d5d00, dBC3Y0d5d10, dBC3Y0d5d01;
	double dBC3Y0d6d00, dBC3Y0d6d10, dBC3Y0d6d01;
	double dBC3Y0d7d00, dBC3Y0d7d10, dBC3Y0d7d01;
	double dBC3Y0d8d00, dBC3Y0d8d10, dBC3Y0d8d01;
	double dBC3Y0d9d00, dBC3Y0d9d10, dBC3Y0d9d01;
	
	double BC4Y0;
	double dBC4Y0d0d00, dBC4Y0d0d10, dBC4Y0d0d01;
	double dBC4Y0d1d00, dBC4Y0d1d10, dBC4Y0d1d01;
	double dBC4Y0d2d00, dBC4Y0d2d10, dBC4Y0d2d01;
	double dBC4Y0d3d00, dBC4Y0d3d10, dBC4Y0d3d01;
	double dBC4Y0d4d00, dBC4Y0d4d10, dBC4Y0d4d01;
	double dBC4Y0d5d00, dBC4Y0d5d10, dBC4Y0d5d01;
	double dBC4Y0d6d00, dBC4Y0d6d10, dBC4Y0d6d01;
	double dBC4Y0d7d00, dBC4Y0d7d10, dBC4Y0d7d01;
	double dBC4Y0d8d00, dBC4Y0d8d10, dBC4Y0d8d01;
	double dBC4Y0d9d00, dBC4Y0d9d10, dBC4Y0d9d01;
	
	double BC5Y0;
	double dBC5Y0d0d00, dBC5Y0d0d10, dBC5Y0d0d01;
	double dBC5Y0d1d00, dBC5Y0d1d10, dBC5Y0d1d01;
	double dBC5Y0d2d00, dBC5Y0d2d10, dBC5Y0d2d01;
	double dBC5Y0d3d00, dBC5Y0d3d10, dBC5Y0d3d01;
	double dBC5Y0d4d00, dBC5Y0d4d10, dBC5Y0d4d01;
	double dBC5Y0d5d00, dBC5Y0d5d10, dBC5Y0d5d01;
	double dBC5Y0d6d00, dBC5Y0d6d10, dBC5Y0d6d01;
	double dBC5Y0d7d00, dBC5Y0d7d10, dBC5Y0d7d01;
	double dBC5Y0d8d00, dBC5Y0d8d10, dBC5Y0d8d01;
	double dBC5Y0d9d00, dBC5Y0d9d10, dBC5Y0d9d01;
	
	double BC6Y0;
	double dBC6Y0d0d00, dBC6Y0d0d10, dBC6Y0d0d01;
	double dBC6Y0d1d00, dBC6Y0d1d10, dBC6Y0d1d01;
	double dBC6Y0d2d00, dBC6Y0d2d10, dBC6Y0d2d01;
	double dBC6Y0d3d00, dBC6Y0d3d10, dBC6Y0d3d01;
	double dBC6Y0d4d00, dBC6Y0d4d10, dBC6Y0d4d01;
	double dBC6Y0d5d00, dBC6Y0d5d10, dBC6Y0d5d01;
	double dBC6Y0d6d00, dBC6Y0d6d10, dBC6Y0d6d01;
	double dBC6Y0d7d00, dBC6Y0d7d10, dBC6Y0d7d01;
	double dBC6Y0d8d00, dBC6Y0d8d10, dBC6Y0d8d01;
	double dBC6Y0d9d00, dBC6Y0d9d10, dBC6Y0d9d01;
	
	double BC7Y0;
	double dBC7Y0d0d00, dBC7Y0d0d10, dBC7Y0d0d01;
	double dBC7Y0d1d00, dBC7Y0d1d10, dBC7Y0d1d01;
	double dBC7Y0d2d00, dBC7Y0d2d10, dBC7Y0d2d01;
	double dBC7Y0d3d00, dBC7Y0d3d10, dBC7Y0d3d01;
	double dBC7Y0d4d00, dBC7Y0d4d10, dBC7Y0d4d01;
	double dBC7Y0d5d00, dBC7Y0d5d10, dBC7Y0d5d01;
	double dBC7Y0d6d00, dBC7Y0d6d10, dBC7Y0d6d01;
	double dBC7Y0d7d00, dBC7Y0d7d10, dBC7Y0d7d01;
	double dBC7Y0d8d00, dBC7Y0d8d10, dBC7Y0d8d01;
	double dBC7Y0d9d00, dBC7Y0d9d10, dBC7Y0d9d01;
	
	double BC8Y0;
	double dBC8Y0d0d00, dBC8Y0d0d10, dBC8Y0d0d01;
	double dBC8Y0d1d00, dBC8Y0d1d10, dBC8Y0d1d01;
	double dBC8Y0d2d00, dBC8Y0d2d10, dBC8Y0d2d01;
	double dBC8Y0d3d00, dBC8Y0d3d10, dBC8Y0d3d01;
	double dBC8Y0d4d00, dBC8Y0d4d10, dBC8Y0d4d01;
	double dBC8Y0d5d00, dBC8Y0d5d10, dBC8Y0d5d01;
	double dBC8Y0d6d00, dBC8Y0d6d10, dBC8Y0d6d01;
	double dBC8Y0d7d00, dBC8Y0d7d10, dBC8Y0d7d01;
	double dBC8Y0d8d00, dBC8Y0d8d10, dBC8Y0d8d01;
	double dBC8Y0d9d00, dBC8Y0d9d10, dBC8Y0d9d01;
	
	double BC9Y0;
	double dBC9Y0d0d00, dBC9Y0d0d10, dBC9Y0d0d01;
	double dBC9Y0d1d00, dBC9Y0d1d10, dBC9Y0d1d01;
	double dBC9Y0d2d00, dBC9Y0d2d10, dBC9Y0d2d01;
	double dBC9Y0d3d00, dBC9Y0d3d10, dBC9Y0d3d01;
	double dBC9Y0d4d00, dBC9Y0d4d10, dBC9Y0d4d01;
	double dBC9Y0d5d00, dBC9Y0d5d10, dBC9Y0d5d01;
	double dBC9Y0d6d00, dBC9Y0d6d10, dBC9Y0d6d01;
	double dBC9Y0d7d00, dBC9Y0d7d10, dBC9Y0d7d01;
	double dBC9Y0d8d00, dBC9Y0d8d10, dBC9Y0d8d01;
	double dBC9Y0d9d00, dBC9Y0d9d10, dBC9Y0d9d01;
	
	nnztemp = 0;
	
	
	//:::::::::::::::::::: X0 ::::::::::::::::::::
	for (i = 0; i < m; ++i)
	{
		j = 0;
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
			
			//----- BC0 -----
			BC0X0       =       BC0X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d0d00 = dBC0X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d0d10 = dBC0X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d0d01 = dBC0X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d1d00 = dBC0X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d1d10 = dBC0X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d1d01 = dBC0X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d2d00 = dBC0X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d2d10 = dBC0X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d2d01 = dBC0X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d3d00 = dBC0X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d3d10 = dBC0X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d3d01 = dBC0X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d4d00 = dBC0X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d4d10 = dBC0X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d4d01 = dBC0X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d5d00 = dBC0X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d5d10 = dBC0X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d5d01 = dBC0X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d6d00 = dBC0X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d6d10 = dBC0X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d6d01 = dBC0X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d7d00 = dBC0X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d7d10 = dBC0X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d7d01 = dBC0X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d8d00 = dBC0X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d8d10 = dBC0X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d8d01 = dBC0X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X0d9d00 = dBC0X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d9d10 = dBC0X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X0d9d01 = dBC0X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 0;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC0X0d0d00*a0d0[a] + dBC0X0d0d10*a0d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC0X0);
			Dxvec[k +pn] = dBC0X0d0d10*u0d10dk + dBC0X0d0d00*u0d00dk + dBC0X0d1d10*u1d10dk + dBC0X0d1d00*u1d00dk + dBC0X0d2d10*u2d10dk + dBC0X0d2d00*u2d00dk + dBC0X0d3d10*u3d10dk + dBC0X0d3d00*u3d00dk + dBC0X0d4d10*u4d10dk + dBC0X0d4d00*u4d00dk + dBC0X0d5d10*u5d10dk + dBC0X0d5d00*u5d00dk + dBC0X0d6d10*u6d10dk + dBC0X0d6d00*u6d00dk + dBC0X0d7d10*u7d10dk + dBC0X0d7d00*u7d00dk + dBC0X0d8d10*u8d10dk + dBC0X0d8d00*u8d00dk + dBC0X0d9d10*u9d10dk + dBC0X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC0X0d0d01*u0d01dk + dBC0X0d0d00*u0d00dk + dBC0X0d1d01*u1d01dk + dBC0X0d1d00*u1d00dk + dBC0X0d2d01*u2d01dk + dBC0X0d2d00*u2d00dk + dBC0X0d3d01*u3d01dk + dBC0X0d3d00*u3d00dk + dBC0X0d4d01*u4d01dk + dBC0X0d4d00*u4d00dk + dBC0X0d5d01*u5d01dk + dBC0X0d5d00*u5d00dk + dBC0X0d6d01*u6d01dk + dBC0X0d6d00*u6d00dk + dBC0X0d7d01*u7d01dk + dBC0X0d7d00*u7d00dk + dBC0X0d8d01*u8d01dk + dBC0X0d8d00*u8d00dk + dBC0X0d9d01*u9d01dk + dBC0X0d9d00*u9d00dk;
			
			//----- BC1 -----
			BC1X0       =       BC1X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d0d00 = dBC1X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d0d10 = dBC1X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d0d01 = dBC1X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d1d00 = dBC1X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d1d10 = dBC1X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d1d01 = dBC1X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d2d00 = dBC1X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d2d10 = dBC1X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d2d01 = dBC1X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d3d00 = dBC1X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d3d10 = dBC1X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d3d01 = dBC1X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d4d00 = dBC1X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d4d10 = dBC1X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d4d01 = dBC1X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d5d00 = dBC1X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d5d10 = dBC1X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d5d01 = dBC1X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d6d00 = dBC1X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d6d10 = dBC1X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d6d01 = dBC1X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d7d00 = dBC1X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d7d10 = dBC1X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d7d01 = dBC1X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d8d00 = dBC1X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d8d10 = dBC1X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d8d01 = dBC1X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X0d9d00 = dBC1X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d9d10 = dBC1X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X0d9d01 = dBC1X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 1;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC1X0d1d00*a1d0[a] + dBC1X0d1d10*a1d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC1X0);
			Dxvec[k +pn] = dBC1X0d0d10*u0d10dk + dBC1X0d0d00*u0d00dk + dBC1X0d1d10*u1d10dk + dBC1X0d1d00*u1d00dk + dBC1X0d2d10*u2d10dk + dBC1X0d2d00*u2d00dk + dBC1X0d3d10*u3d10dk + dBC1X0d3d00*u3d00dk + dBC1X0d4d10*u4d10dk + dBC1X0d4d00*u4d00dk + dBC1X0d5d10*u5d10dk + dBC1X0d5d00*u5d00dk + dBC1X0d6d10*u6d10dk + dBC1X0d6d00*u6d00dk + dBC1X0d7d10*u7d10dk + dBC1X0d7d00*u7d00dk + dBC1X0d8d10*u8d10dk + dBC1X0d8d00*u8d00dk + dBC1X0d9d10*u9d10dk + dBC1X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC1X0d0d01*u0d01dk + dBC1X0d0d00*u0d00dk + dBC1X0d1d01*u1d01dk + dBC1X0d1d00*u1d00dk + dBC1X0d2d01*u2d01dk + dBC1X0d2d00*u2d00dk + dBC1X0d3d01*u3d01dk + dBC1X0d3d00*u3d00dk + dBC1X0d4d01*u4d01dk + dBC1X0d4d00*u4d00dk + dBC1X0d5d01*u5d01dk + dBC1X0d5d00*u5d00dk + dBC1X0d6d01*u6d01dk + dBC1X0d6d00*u6d00dk + dBC1X0d7d01*u7d01dk + dBC1X0d7d00*u7d00dk + dBC1X0d8d01*u8d01dk + dBC1X0d8d00*u8d00dk + dBC1X0d9d01*u9d01dk + dBC1X0d9d00*u9d00dk;
			
			//----- BC2 -----
			BC2X0       =       BC2X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d0d00 = dBC2X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d0d10 = dBC2X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d0d01 = dBC2X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d1d00 = dBC2X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d1d10 = dBC2X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d1d01 = dBC2X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d2d00 = dBC2X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d2d10 = dBC2X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d2d01 = dBC2X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d3d00 = dBC2X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d3d10 = dBC2X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d3d01 = dBC2X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d4d00 = dBC2X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d4d10 = dBC2X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d4d01 = dBC2X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d5d00 = dBC2X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d5d10 = dBC2X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d5d01 = dBC2X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d6d00 = dBC2X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d6d10 = dBC2X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d6d01 = dBC2X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d7d00 = dBC2X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d7d10 = dBC2X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d7d01 = dBC2X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d8d00 = dBC2X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d8d10 = dBC2X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d8d01 = dBC2X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X0d9d00 = dBC2X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d9d10 = dBC2X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X0d9d01 = dBC2X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 2;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC2X0d2d00*a2d0[a] + dBC2X0d2d10*a2d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC2X0);
			Dxvec[k +pn] = dBC2X0d0d10*u0d10dk + dBC2X0d0d00*u0d00dk + dBC2X0d1d10*u1d10dk + dBC2X0d1d00*u1d00dk + dBC2X0d2d10*u2d10dk + dBC2X0d2d00*u2d00dk + dBC2X0d3d10*u3d10dk + dBC2X0d3d00*u3d00dk + dBC2X0d4d10*u4d10dk + dBC2X0d4d00*u4d00dk + dBC2X0d5d10*u5d10dk + dBC2X0d5d00*u5d00dk + dBC2X0d6d10*u6d10dk + dBC2X0d6d00*u6d00dk + dBC2X0d7d10*u7d10dk + dBC2X0d7d00*u7d00dk + dBC2X0d8d10*u8d10dk + dBC2X0d8d00*u8d00dk + dBC2X0d9d10*u9d10dk + dBC2X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC2X0d0d01*u0d01dk + dBC2X0d0d00*u0d00dk + dBC2X0d1d01*u1d01dk + dBC2X0d1d00*u1d00dk + dBC2X0d2d01*u2d01dk + dBC2X0d2d00*u2d00dk + dBC2X0d3d01*u3d01dk + dBC2X0d3d00*u3d00dk + dBC2X0d4d01*u4d01dk + dBC2X0d4d00*u4d00dk + dBC2X0d5d01*u5d01dk + dBC2X0d5d00*u5d00dk + dBC2X0d6d01*u6d01dk + dBC2X0d6d00*u6d00dk + dBC2X0d7d01*u7d01dk + dBC2X0d7d00*u7d00dk + dBC2X0d8d01*u8d01dk + dBC2X0d8d00*u8d00dk + dBC2X0d9d01*u9d01dk + dBC2X0d9d00*u9d00dk;
			
			//----- BC3 -----
			BC3X0       =       BC3X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d0d00 = dBC3X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d0d10 = dBC3X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d0d01 = dBC3X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d1d00 = dBC3X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d1d10 = dBC3X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d1d01 = dBC3X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d2d00 = dBC3X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d2d10 = dBC3X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d2d01 = dBC3X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d3d00 = dBC3X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d3d10 = dBC3X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d3d01 = dBC3X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d4d00 = dBC3X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d4d10 = dBC3X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d4d01 = dBC3X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d5d00 = dBC3X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d5d10 = dBC3X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d5d01 = dBC3X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d6d00 = dBC3X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d6d10 = dBC3X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d6d01 = dBC3X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d7d00 = dBC3X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d7d10 = dBC3X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d7d01 = dBC3X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d8d00 = dBC3X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d8d10 = dBC3X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d8d01 = dBC3X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X0d9d00 = dBC3X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d9d10 = dBC3X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X0d9d01 = dBC3X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 3;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC3X0d3d00*a3d0[a] + dBC3X0d3d10*a3d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC3X0);
			Dxvec[k +pn] = dBC3X0d0d10*u0d10dk + dBC3X0d0d00*u0d00dk + dBC3X0d1d10*u1d10dk + dBC3X0d1d00*u1d00dk + dBC3X0d2d10*u2d10dk + dBC3X0d2d00*u2d00dk + dBC3X0d3d10*u3d10dk + dBC3X0d3d00*u3d00dk + dBC3X0d4d10*u4d10dk + dBC3X0d4d00*u4d00dk + dBC3X0d5d10*u5d10dk + dBC3X0d5d00*u5d00dk + dBC3X0d6d10*u6d10dk + dBC3X0d6d00*u6d00dk + dBC3X0d7d10*u7d10dk + dBC3X0d7d00*u7d00dk + dBC3X0d8d10*u8d10dk + dBC3X0d8d00*u8d00dk + dBC3X0d9d10*u9d10dk + dBC3X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC3X0d0d01*u0d01dk + dBC3X0d0d00*u0d00dk + dBC3X0d1d01*u1d01dk + dBC3X0d1d00*u1d00dk + dBC3X0d2d01*u2d01dk + dBC3X0d2d00*u2d00dk + dBC3X0d3d01*u3d01dk + dBC3X0d3d00*u3d00dk + dBC3X0d4d01*u4d01dk + dBC3X0d4d00*u4d00dk + dBC3X0d5d01*u5d01dk + dBC3X0d5d00*u5d00dk + dBC3X0d6d01*u6d01dk + dBC3X0d6d00*u6d00dk + dBC3X0d7d01*u7d01dk + dBC3X0d7d00*u7d00dk + dBC3X0d8d01*u8d01dk + dBC3X0d8d00*u8d00dk + dBC3X0d9d01*u9d01dk + dBC3X0d9d00*u9d00dk;
			
			//----- BC4 -----
			BC4X0       =       BC4X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d0d00 = dBC4X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d0d10 = dBC4X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d0d01 = dBC4X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d1d00 = dBC4X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d1d10 = dBC4X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d1d01 = dBC4X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d2d00 = dBC4X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d2d10 = dBC4X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d2d01 = dBC4X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d3d00 = dBC4X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d3d10 = dBC4X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d3d01 = dBC4X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d4d00 = dBC4X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d4d10 = dBC4X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d4d01 = dBC4X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d5d00 = dBC4X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d5d10 = dBC4X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d5d01 = dBC4X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d6d00 = dBC4X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d6d10 = dBC4X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d6d01 = dBC4X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d7d00 = dBC4X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d7d10 = dBC4X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d7d01 = dBC4X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d8d00 = dBC4X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d8d10 = dBC4X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d8d01 = dBC4X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X0d9d00 = dBC4X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d9d10 = dBC4X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X0d9d01 = dBC4X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 4;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC4X0d4d00*a4d0[a] + dBC4X0d4d10*a4d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC4X0);
			Dxvec[k +pn] = dBC4X0d0d10*u0d10dk + dBC4X0d0d00*u0d00dk + dBC4X0d1d10*u1d10dk + dBC4X0d1d00*u1d00dk + dBC4X0d2d10*u2d10dk + dBC4X0d2d00*u2d00dk + dBC4X0d3d10*u3d10dk + dBC4X0d3d00*u3d00dk + dBC4X0d4d10*u4d10dk + dBC4X0d4d00*u4d00dk + dBC4X0d5d10*u5d10dk + dBC4X0d5d00*u5d00dk + dBC4X0d6d10*u6d10dk + dBC4X0d6d00*u6d00dk + dBC4X0d7d10*u7d10dk + dBC4X0d7d00*u7d00dk + dBC4X0d8d10*u8d10dk + dBC4X0d8d00*u8d00dk + dBC4X0d9d10*u9d10dk + dBC4X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC4X0d0d01*u0d01dk + dBC4X0d0d00*u0d00dk + dBC4X0d1d01*u1d01dk + dBC4X0d1d00*u1d00dk + dBC4X0d2d01*u2d01dk + dBC4X0d2d00*u2d00dk + dBC4X0d3d01*u3d01dk + dBC4X0d3d00*u3d00dk + dBC4X0d4d01*u4d01dk + dBC4X0d4d00*u4d00dk + dBC4X0d5d01*u5d01dk + dBC4X0d5d00*u5d00dk + dBC4X0d6d01*u6d01dk + dBC4X0d6d00*u6d00dk + dBC4X0d7d01*u7d01dk + dBC4X0d7d00*u7d00dk + dBC4X0d8d01*u8d01dk + dBC4X0d8d00*u8d00dk + dBC4X0d9d01*u9d01dk + dBC4X0d9d00*u9d00dk;
			
			//----- BC5 -----
			BC5X0       =       BC5X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d0d00 = dBC5X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d0d10 = dBC5X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d0d01 = dBC5X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d1d00 = dBC5X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d1d10 = dBC5X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d1d01 = dBC5X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d2d00 = dBC5X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d2d10 = dBC5X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d2d01 = dBC5X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d3d00 = dBC5X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d3d10 = dBC5X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d3d01 = dBC5X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d4d00 = dBC5X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d4d10 = dBC5X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d4d01 = dBC5X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d5d00 = dBC5X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d5d10 = dBC5X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d5d01 = dBC5X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d6d00 = dBC5X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d6d10 = dBC5X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d6d01 = dBC5X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d7d00 = dBC5X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d7d10 = dBC5X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d7d01 = dBC5X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d8d00 = dBC5X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d8d10 = dBC5X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d8d01 = dBC5X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X0d9d00 = dBC5X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d9d10 = dBC5X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X0d9d01 = dBC5X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 5;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC5X0d5d00*a5d0[a] + dBC5X0d5d10*a5d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC5X0d0d01*b0d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC5X0);
			Dxvec[k +pn] = dBC5X0d0d10*u0d10dk + dBC5X0d0d00*u0d00dk + dBC5X0d1d10*u1d10dk + dBC5X0d1d00*u1d00dk + dBC5X0d2d10*u2d10dk + dBC5X0d2d00*u2d00dk + dBC5X0d3d10*u3d10dk + dBC5X0d3d00*u3d00dk + dBC5X0d4d10*u4d10dk + dBC5X0d4d00*u4d00dk + dBC5X0d5d10*u5d10dk + dBC5X0d5d00*u5d00dk + dBC5X0d6d10*u6d10dk + dBC5X0d6d00*u6d00dk + dBC5X0d7d10*u7d10dk + dBC5X0d7d00*u7d00dk + dBC5X0d8d10*u8d10dk + dBC5X0d8d00*u8d00dk + dBC5X0d9d10*u9d10dk + dBC5X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC5X0d0d01*u0d01dk + dBC5X0d0d00*u0d00dk + dBC5X0d1d01*u1d01dk + dBC5X0d1d00*u1d00dk + dBC5X0d2d01*u2d01dk + dBC5X0d2d00*u2d00dk + dBC5X0d3d01*u3d01dk + dBC5X0d3d00*u3d00dk + dBC5X0d4d01*u4d01dk + dBC5X0d4d00*u4d00dk + dBC5X0d5d01*u5d01dk + dBC5X0d5d00*u5d00dk + dBC5X0d6d01*u6d01dk + dBC5X0d6d00*u6d00dk + dBC5X0d7d01*u7d01dk + dBC5X0d7d00*u7d00dk + dBC5X0d8d01*u8d01dk + dBC5X0d8d00*u8d00dk + dBC5X0d9d01*u9d01dk + dBC5X0d9d00*u9d00dk;
			
			//----- BC6 -----
			BC6X0       =       BC6X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d0d00 = dBC6X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d0d10 = dBC6X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d0d01 = dBC6X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d1d00 = dBC6X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d1d10 = dBC6X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d1d01 = dBC6X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d2d00 = dBC6X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d2d10 = dBC6X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d2d01 = dBC6X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d3d00 = dBC6X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d3d10 = dBC6X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d3d01 = dBC6X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d4d00 = dBC6X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d4d10 = dBC6X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d4d01 = dBC6X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d5d00 = dBC6X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d5d10 = dBC6X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d5d01 = dBC6X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d6d00 = dBC6X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d6d10 = dBC6X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d6d01 = dBC6X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d7d00 = dBC6X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d7d10 = dBC6X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d7d01 = dBC6X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d8d00 = dBC6X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d8d10 = dBC6X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d8d01 = dBC6X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X0d9d00 = dBC6X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d9d10 = dBC6X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X0d9d01 = dBC6X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 6;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC6X0d6d00*a6d0[a] + dBC6X0d6d10*a6d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +6;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC6X0d1d01*b1d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC6X0);
			Dxvec[k +pn] = dBC6X0d0d10*u0d10dk + dBC6X0d0d00*u0d00dk + dBC6X0d1d10*u1d10dk + dBC6X0d1d00*u1d00dk + dBC6X0d2d10*u2d10dk + dBC6X0d2d00*u2d00dk + dBC6X0d3d10*u3d10dk + dBC6X0d3d00*u3d00dk + dBC6X0d4d10*u4d10dk + dBC6X0d4d00*u4d00dk + dBC6X0d5d10*u5d10dk + dBC6X0d5d00*u5d00dk + dBC6X0d6d10*u6d10dk + dBC6X0d6d00*u6d00dk + dBC6X0d7d10*u7d10dk + dBC6X0d7d00*u7d00dk + dBC6X0d8d10*u8d10dk + dBC6X0d8d00*u8d00dk + dBC6X0d9d10*u9d10dk + dBC6X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC6X0d0d01*u0d01dk + dBC6X0d0d00*u0d00dk + dBC6X0d1d01*u1d01dk + dBC6X0d1d00*u1d00dk + dBC6X0d2d01*u2d01dk + dBC6X0d2d00*u2d00dk + dBC6X0d3d01*u3d01dk + dBC6X0d3d00*u3d00dk + dBC6X0d4d01*u4d01dk + dBC6X0d4d00*u4d00dk + dBC6X0d5d01*u5d01dk + dBC6X0d5d00*u5d00dk + dBC6X0d6d01*u6d01dk + dBC6X0d6d00*u6d00dk + dBC6X0d7d01*u7d01dk + dBC6X0d7d00*u7d00dk + dBC6X0d8d01*u8d01dk + dBC6X0d8d00*u8d00dk + dBC6X0d9d01*u9d01dk + dBC6X0d9d00*u9d00dk;
			
			//----- BC7 -----
			BC7X0       =       BC7X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d0d00 = dBC7X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d0d10 = dBC7X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d0d01 = dBC7X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d1d00 = dBC7X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d1d10 = dBC7X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d1d01 = dBC7X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d2d00 = dBC7X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d2d10 = dBC7X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d2d01 = dBC7X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d3d00 = dBC7X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d3d10 = dBC7X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d3d01 = dBC7X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d4d00 = dBC7X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d4d10 = dBC7X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d4d01 = dBC7X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d5d00 = dBC7X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d5d10 = dBC7X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d5d01 = dBC7X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d6d00 = dBC7X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d6d10 = dBC7X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d6d01 = dBC7X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d7d00 = dBC7X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d7d10 = dBC7X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d7d01 = dBC7X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d8d00 = dBC7X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d8d10 = dBC7X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d8d01 = dBC7X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X0d9d00 = dBC7X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d9d10 = dBC7X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X0d9d01 = dBC7X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 7;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC7X0d7d00*a7d0[a] + dBC7X0d7d10*a7d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC7X0d2d01*b2d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC7X0);
			Dxvec[k +pn] = dBC7X0d0d10*u0d10dk + dBC7X0d0d00*u0d00dk + dBC7X0d1d10*u1d10dk + dBC7X0d1d00*u1d00dk + dBC7X0d2d10*u2d10dk + dBC7X0d2d00*u2d00dk + dBC7X0d3d10*u3d10dk + dBC7X0d3d00*u3d00dk + dBC7X0d4d10*u4d10dk + dBC7X0d4d00*u4d00dk + dBC7X0d5d10*u5d10dk + dBC7X0d5d00*u5d00dk + dBC7X0d6d10*u6d10dk + dBC7X0d6d00*u6d00dk + dBC7X0d7d10*u7d10dk + dBC7X0d7d00*u7d00dk + dBC7X0d8d10*u8d10dk + dBC7X0d8d00*u8d00dk + dBC7X0d9d10*u9d10dk + dBC7X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC7X0d0d01*u0d01dk + dBC7X0d0d00*u0d00dk + dBC7X0d1d01*u1d01dk + dBC7X0d1d00*u1d00dk + dBC7X0d2d01*u2d01dk + dBC7X0d2d00*u2d00dk + dBC7X0d3d01*u3d01dk + dBC7X0d3d00*u3d00dk + dBC7X0d4d01*u4d01dk + dBC7X0d4d00*u4d00dk + dBC7X0d5d01*u5d01dk + dBC7X0d5d00*u5d00dk + dBC7X0d6d01*u6d01dk + dBC7X0d6d00*u6d00dk + dBC7X0d7d01*u7d01dk + dBC7X0d7d00*u7d00dk + dBC7X0d8d01*u8d01dk + dBC7X0d8d00*u8d00dk + dBC7X0d9d01*u9d01dk + dBC7X0d9d00*u9d00dk;
			
			//----- BC8 -----
			BC8X0       =       BC8X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d0d00 = dBC8X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d0d10 = dBC8X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d0d01 = dBC8X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d1d00 = dBC8X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d1d10 = dBC8X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d1d01 = dBC8X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d2d00 = dBC8X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d2d10 = dBC8X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d2d01 = dBC8X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d3d00 = dBC8X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d3d10 = dBC8X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d3d01 = dBC8X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d4d00 = dBC8X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d4d10 = dBC8X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d4d01 = dBC8X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d5d00 = dBC8X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d5d10 = dBC8X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d5d01 = dBC8X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d6d00 = dBC8X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d6d10 = dBC8X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d6d01 = dBC8X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d7d00 = dBC8X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d7d10 = dBC8X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d7d01 = dBC8X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d8d00 = dBC8X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d8d10 = dBC8X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d8d01 = dBC8X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X0d9d00 = dBC8X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d9d10 = dBC8X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X0d9d01 = dBC8X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 8;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC8X0d8d00*a8d0[a] + dBC8X0d8d10*a8d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC8X0d3d01*b3d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC8X0);
			Dxvec[k +pn] = dBC8X0d0d10*u0d10dk + dBC8X0d0d00*u0d00dk + dBC8X0d1d10*u1d10dk + dBC8X0d1d00*u1d00dk + dBC8X0d2d10*u2d10dk + dBC8X0d2d00*u2d00dk + dBC8X0d3d10*u3d10dk + dBC8X0d3d00*u3d00dk + dBC8X0d4d10*u4d10dk + dBC8X0d4d00*u4d00dk + dBC8X0d5d10*u5d10dk + dBC8X0d5d00*u5d00dk + dBC8X0d6d10*u6d10dk + dBC8X0d6d00*u6d00dk + dBC8X0d7d10*u7d10dk + dBC8X0d7d00*u7d00dk + dBC8X0d8d10*u8d10dk + dBC8X0d8d00*u8d00dk + dBC8X0d9d10*u9d10dk + dBC8X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC8X0d0d01*u0d01dk + dBC8X0d0d00*u0d00dk + dBC8X0d1d01*u1d01dk + dBC8X0d1d00*u1d00dk + dBC8X0d2d01*u2d01dk + dBC8X0d2d00*u2d00dk + dBC8X0d3d01*u3d01dk + dBC8X0d3d00*u3d00dk + dBC8X0d4d01*u4d01dk + dBC8X0d4d00*u4d00dk + dBC8X0d5d01*u5d01dk + dBC8X0d5d00*u5d00dk + dBC8X0d6d01*u6d01dk + dBC8X0d6d00*u6d00dk + dBC8X0d7d01*u7d01dk + dBC8X0d7d00*u7d00dk + dBC8X0d8d01*u8d01dk + dBC8X0d8d00*u8d00dk + dBC8X0d9d01*u9d01dk + dBC8X0d9d00*u9d00dk;
			
			//----- BC9 -----
			BC9X0       =       BC9X0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d0d00 = dBC9X0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d0d10 = dBC9X0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d0d01 = dBC9X0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d1d00 = dBC9X0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d1d10 = dBC9X0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d1d01 = dBC9X0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d2d00 = dBC9X0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d2d10 = dBC9X0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d2d01 = dBC9X0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d3d00 = dBC9X0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d3d10 = dBC9X0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d3d01 = dBC9X0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d4d00 = dBC9X0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d4d10 = dBC9X0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d4d01 = dBC9X0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d5d00 = dBC9X0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d5d10 = dBC9X0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d5d01 = dBC9X0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d6d00 = dBC9X0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d6d10 = dBC9X0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d6d01 = dBC9X0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d7d00 = dBC9X0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d7d10 = dBC9X0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d7d01 = dBC9X0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d8d00 = dBC9X0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d8d10 = dBC9X0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d8d01 = dBC9X0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X0d9d00 = dBC9X0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d9d10 = dBC9X0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X0d9d01 = dBC9X0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 9;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC9X0d9d00*a9d0[a] + dBC9X0d9d10*a9d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC9X0d4d01*b4d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC9X0);
			Dxvec[k +pn] = dBC9X0d0d10*u0d10dk + dBC9X0d0d00*u0d00dk + dBC9X0d1d10*u1d10dk + dBC9X0d1d00*u1d00dk + dBC9X0d2d10*u2d10dk + dBC9X0d2d00*u2d00dk + dBC9X0d3d10*u3d10dk + dBC9X0d3d00*u3d00dk + dBC9X0d4d10*u4d10dk + dBC9X0d4d00*u4d00dk + dBC9X0d5d10*u5d10dk + dBC9X0d5d00*u5d00dk + dBC9X0d6d10*u6d10dk + dBC9X0d6d00*u6d00dk + dBC9X0d7d10*u7d10dk + dBC9X0d7d00*u7d00dk + dBC9X0d8d10*u8d10dk + dBC9X0d8d00*u8d00dk + dBC9X0d9d10*u9d10dk + dBC9X0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC9X0d0d01*u0d01dk + dBC9X0d0d00*u0d00dk + dBC9X0d1d01*u1d01dk + dBC9X0d1d00*u1d00dk + dBC9X0d2d01*u2d01dk + dBC9X0d2d00*u2d00dk + dBC9X0d3d01*u3d01dk + dBC9X0d3d00*u3d00dk + dBC9X0d4d01*u4d01dk + dBC9X0d4d00*u4d00dk + dBC9X0d5d01*u5d01dk + dBC9X0d5d00*u5d00dk + dBC9X0d6d01*u6d01dk + dBC9X0d6d00*u6d00dk + dBC9X0d7d01*u7d01dk + dBC9X0d7d00*u7d00dk + dBC9X0d8d01*u8d01dk + dBC9X0d8d00*u8d00dk + dBC9X0d9d01*u9d01dk + dBC9X0d9d00*u9d00dk;
			
		}
	}
	
	//:::::::::::::::::::: X1 ::::::::::::::::::::
	for (i = 0; i < m; ++i)
	{
		j = n-1;
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
			
			//----- BC0 -----
			BC0X1       =       BC0X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d0d00 = dBC0X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d0d10 = dBC0X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d0d01 = dBC0X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d1d00 = dBC0X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d1d10 = dBC0X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d1d01 = dBC0X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d2d00 = dBC0X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d2d10 = dBC0X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d2d01 = dBC0X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d3d00 = dBC0X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d3d10 = dBC0X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d3d01 = dBC0X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d4d00 = dBC0X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d4d10 = dBC0X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d4d01 = dBC0X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d5d00 = dBC0X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d5d10 = dBC0X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d5d01 = dBC0X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d6d00 = dBC0X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d6d10 = dBC0X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d6d01 = dBC0X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d7d00 = dBC0X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d7d10 = dBC0X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d7d01 = dBC0X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d8d00 = dBC0X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d8d10 = dBC0X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d8d01 = dBC0X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0X1d9d00 = dBC0X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d9d10 = dBC0X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0X1d9d01 = dBC0X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 0;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC0X1d0d00*a0d0[a] + dBC0X1d0d10*a0d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +0;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC0X1);
			Dxvec[k +pn] = dBC0X1d0d10*u0d10dk + dBC0X1d0d00*u0d00dk + dBC0X1d1d10*u1d10dk + dBC0X1d1d00*u1d00dk + dBC0X1d2d10*u2d10dk + dBC0X1d2d00*u2d00dk + dBC0X1d3d10*u3d10dk + dBC0X1d3d00*u3d00dk + dBC0X1d4d10*u4d10dk + dBC0X1d4d00*u4d00dk + dBC0X1d5d10*u5d10dk + dBC0X1d5d00*u5d00dk + dBC0X1d6d10*u6d10dk + dBC0X1d6d00*u6d00dk + dBC0X1d7d10*u7d10dk + dBC0X1d7d00*u7d00dk + dBC0X1d8d10*u8d10dk + dBC0X1d8d00*u8d00dk + dBC0X1d9d10*u9d10dk + dBC0X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC0X1d0d01*u0d01dk + dBC0X1d0d00*u0d00dk + dBC0X1d1d01*u1d01dk + dBC0X1d1d00*u1d00dk + dBC0X1d2d01*u2d01dk + dBC0X1d2d00*u2d00dk + dBC0X1d3d01*u3d01dk + dBC0X1d3d00*u3d00dk + dBC0X1d4d01*u4d01dk + dBC0X1d4d00*u4d00dk + dBC0X1d5d01*u5d01dk + dBC0X1d5d00*u5d00dk + dBC0X1d6d01*u6d01dk + dBC0X1d6d00*u6d00dk + dBC0X1d7d01*u7d01dk + dBC0X1d7d00*u7d00dk + dBC0X1d8d01*u8d01dk + dBC0X1d8d00*u8d00dk + dBC0X1d9d01*u9d01dk + dBC0X1d9d00*u9d00dk;
			
			//----- BC1 -----
			BC1X1       =       BC1X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d0d00 = dBC1X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d0d10 = dBC1X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d0d01 = dBC1X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d1d00 = dBC1X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d1d10 = dBC1X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d1d01 = dBC1X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d2d00 = dBC1X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d2d10 = dBC1X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d2d01 = dBC1X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d3d00 = dBC1X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d3d10 = dBC1X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d3d01 = dBC1X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d4d00 = dBC1X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d4d10 = dBC1X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d4d01 = dBC1X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d5d00 = dBC1X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d5d10 = dBC1X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d5d01 = dBC1X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d6d00 = dBC1X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d6d10 = dBC1X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d6d01 = dBC1X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d7d00 = dBC1X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d7d10 = dBC1X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d7d01 = dBC1X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d8d00 = dBC1X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d8d10 = dBC1X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d8d01 = dBC1X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1X1d9d00 = dBC1X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d9d10 = dBC1X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1X1d9d01 = dBC1X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 1;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC1X1d1d00*a1d0[a] + dBC1X1d1d10*a1d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +1;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC1X1);
			Dxvec[k +pn] = dBC1X1d0d10*u0d10dk + dBC1X1d0d00*u0d00dk + dBC1X1d1d10*u1d10dk + dBC1X1d1d00*u1d00dk + dBC1X1d2d10*u2d10dk + dBC1X1d2d00*u2d00dk + dBC1X1d3d10*u3d10dk + dBC1X1d3d00*u3d00dk + dBC1X1d4d10*u4d10dk + dBC1X1d4d00*u4d00dk + dBC1X1d5d10*u5d10dk + dBC1X1d5d00*u5d00dk + dBC1X1d6d10*u6d10dk + dBC1X1d6d00*u6d00dk + dBC1X1d7d10*u7d10dk + dBC1X1d7d00*u7d00dk + dBC1X1d8d10*u8d10dk + dBC1X1d8d00*u8d00dk + dBC1X1d9d10*u9d10dk + dBC1X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC1X1d0d01*u0d01dk + dBC1X1d0d00*u0d00dk + dBC1X1d1d01*u1d01dk + dBC1X1d1d00*u1d00dk + dBC1X1d2d01*u2d01dk + dBC1X1d2d00*u2d00dk + dBC1X1d3d01*u3d01dk + dBC1X1d3d00*u3d00dk + dBC1X1d4d01*u4d01dk + dBC1X1d4d00*u4d00dk + dBC1X1d5d01*u5d01dk + dBC1X1d5d00*u5d00dk + dBC1X1d6d01*u6d01dk + dBC1X1d6d00*u6d00dk + dBC1X1d7d01*u7d01dk + dBC1X1d7d00*u7d00dk + dBC1X1d8d01*u8d01dk + dBC1X1d8d00*u8d00dk + dBC1X1d9d01*u9d01dk + dBC1X1d9d00*u9d00dk;
			
			//----- BC2 -----
			BC2X1       =       BC2X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d0d00 = dBC2X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d0d10 = dBC2X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d0d01 = dBC2X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d1d00 = dBC2X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d1d10 = dBC2X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d1d01 = dBC2X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d2d00 = dBC2X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d2d10 = dBC2X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d2d01 = dBC2X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d3d00 = dBC2X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d3d10 = dBC2X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d3d01 = dBC2X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d4d00 = dBC2X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d4d10 = dBC2X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d4d01 = dBC2X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d5d00 = dBC2X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d5d10 = dBC2X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d5d01 = dBC2X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d6d00 = dBC2X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d6d10 = dBC2X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d6d01 = dBC2X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d7d00 = dBC2X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d7d10 = dBC2X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d7d01 = dBC2X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d8d00 = dBC2X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d8d10 = dBC2X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d8d01 = dBC2X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2X1d9d00 = dBC2X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d9d10 = dBC2X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2X1d9d01 = dBC2X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 2;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC2X1d2d00*a2d0[a] + dBC2X1d2d10*a2d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +2;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC2X1);
			Dxvec[k +pn] = dBC2X1d0d10*u0d10dk + dBC2X1d0d00*u0d00dk + dBC2X1d1d10*u1d10dk + dBC2X1d1d00*u1d00dk + dBC2X1d2d10*u2d10dk + dBC2X1d2d00*u2d00dk + dBC2X1d3d10*u3d10dk + dBC2X1d3d00*u3d00dk + dBC2X1d4d10*u4d10dk + dBC2X1d4d00*u4d00dk + dBC2X1d5d10*u5d10dk + dBC2X1d5d00*u5d00dk + dBC2X1d6d10*u6d10dk + dBC2X1d6d00*u6d00dk + dBC2X1d7d10*u7d10dk + dBC2X1d7d00*u7d00dk + dBC2X1d8d10*u8d10dk + dBC2X1d8d00*u8d00dk + dBC2X1d9d10*u9d10dk + dBC2X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC2X1d0d01*u0d01dk + dBC2X1d0d00*u0d00dk + dBC2X1d1d01*u1d01dk + dBC2X1d1d00*u1d00dk + dBC2X1d2d01*u2d01dk + dBC2X1d2d00*u2d00dk + dBC2X1d3d01*u3d01dk + dBC2X1d3d00*u3d00dk + dBC2X1d4d01*u4d01dk + dBC2X1d4d00*u4d00dk + dBC2X1d5d01*u5d01dk + dBC2X1d5d00*u5d00dk + dBC2X1d6d01*u6d01dk + dBC2X1d6d00*u6d00dk + dBC2X1d7d01*u7d01dk + dBC2X1d7d00*u7d00dk + dBC2X1d8d01*u8d01dk + dBC2X1d8d00*u8d00dk + dBC2X1d9d01*u9d01dk + dBC2X1d9d00*u9d00dk;
			
			//----- BC3 -----
			BC3X1       =       BC3X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d0d00 = dBC3X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d0d10 = dBC3X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d0d01 = dBC3X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d1d00 = dBC3X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d1d10 = dBC3X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d1d01 = dBC3X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d2d00 = dBC3X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d2d10 = dBC3X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d2d01 = dBC3X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d3d00 = dBC3X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d3d10 = dBC3X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d3d01 = dBC3X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d4d00 = dBC3X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d4d10 = dBC3X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d4d01 = dBC3X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d5d00 = dBC3X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d5d10 = dBC3X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d5d01 = dBC3X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d6d00 = dBC3X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d6d10 = dBC3X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d6d01 = dBC3X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d7d00 = dBC3X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d7d10 = dBC3X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d7d01 = dBC3X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d8d00 = dBC3X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d8d10 = dBC3X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d8d01 = dBC3X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3X1d9d00 = dBC3X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d9d10 = dBC3X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3X1d9d01 = dBC3X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 3;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC3X1d3d00*a3d0[a] + dBC3X1d3d10*a3d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +3;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC3X1);
			Dxvec[k +pn] = dBC3X1d0d10*u0d10dk + dBC3X1d0d00*u0d00dk + dBC3X1d1d10*u1d10dk + dBC3X1d1d00*u1d00dk + dBC3X1d2d10*u2d10dk + dBC3X1d2d00*u2d00dk + dBC3X1d3d10*u3d10dk + dBC3X1d3d00*u3d00dk + dBC3X1d4d10*u4d10dk + dBC3X1d4d00*u4d00dk + dBC3X1d5d10*u5d10dk + dBC3X1d5d00*u5d00dk + dBC3X1d6d10*u6d10dk + dBC3X1d6d00*u6d00dk + dBC3X1d7d10*u7d10dk + dBC3X1d7d00*u7d00dk + dBC3X1d8d10*u8d10dk + dBC3X1d8d00*u8d00dk + dBC3X1d9d10*u9d10dk + dBC3X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC3X1d0d01*u0d01dk + dBC3X1d0d00*u0d00dk + dBC3X1d1d01*u1d01dk + dBC3X1d1d00*u1d00dk + dBC3X1d2d01*u2d01dk + dBC3X1d2d00*u2d00dk + dBC3X1d3d01*u3d01dk + dBC3X1d3d00*u3d00dk + dBC3X1d4d01*u4d01dk + dBC3X1d4d00*u4d00dk + dBC3X1d5d01*u5d01dk + dBC3X1d5d00*u5d00dk + dBC3X1d6d01*u6d01dk + dBC3X1d6d00*u6d00dk + dBC3X1d7d01*u7d01dk + dBC3X1d7d00*u7d00dk + dBC3X1d8d01*u8d01dk + dBC3X1d8d00*u8d00dk + dBC3X1d9d01*u9d01dk + dBC3X1d9d00*u9d00dk;
			
			//----- BC4 -----
			BC4X1       =       BC4X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d0d00 = dBC4X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d0d10 = dBC4X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d0d01 = dBC4X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d1d00 = dBC4X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d1d10 = dBC4X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d1d01 = dBC4X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d2d00 = dBC4X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d2d10 = dBC4X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d2d01 = dBC4X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d3d00 = dBC4X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d3d10 = dBC4X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d3d01 = dBC4X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d4d00 = dBC4X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d4d10 = dBC4X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d4d01 = dBC4X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d5d00 = dBC4X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d5d10 = dBC4X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d5d01 = dBC4X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d6d00 = dBC4X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d6d10 = dBC4X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d6d01 = dBC4X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d7d00 = dBC4X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d7d10 = dBC4X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d7d01 = dBC4X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d8d00 = dBC4X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d8d10 = dBC4X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d8d01 = dBC4X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4X1d9d00 = dBC4X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d9d10 = dBC4X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4X1d9d01 = dBC4X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 4;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC4X1d4d00*a4d0[a] + dBC4X1d4d10*a4d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +4;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
			}
			
			bvec[k +pn]  = sgn*(BC4X1);
			Dxvec[k +pn] = dBC4X1d0d10*u0d10dk + dBC4X1d0d00*u0d00dk + dBC4X1d1d10*u1d10dk + dBC4X1d1d00*u1d00dk + dBC4X1d2d10*u2d10dk + dBC4X1d2d00*u2d00dk + dBC4X1d3d10*u3d10dk + dBC4X1d3d00*u3d00dk + dBC4X1d4d10*u4d10dk + dBC4X1d4d00*u4d00dk + dBC4X1d5d10*u5d10dk + dBC4X1d5d00*u5d00dk + dBC4X1d6d10*u6d10dk + dBC4X1d6d00*u6d00dk + dBC4X1d7d10*u7d10dk + dBC4X1d7d00*u7d00dk + dBC4X1d8d10*u8d10dk + dBC4X1d8d00*u8d00dk + dBC4X1d9d10*u9d10dk + dBC4X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC4X1d0d01*u0d01dk + dBC4X1d0d00*u0d00dk + dBC4X1d1d01*u1d01dk + dBC4X1d1d00*u1d00dk + dBC4X1d2d01*u2d01dk + dBC4X1d2d00*u2d00dk + dBC4X1d3d01*u3d01dk + dBC4X1d3d00*u3d00dk + dBC4X1d4d01*u4d01dk + dBC4X1d4d00*u4d00dk + dBC4X1d5d01*u5d01dk + dBC4X1d5d00*u5d00dk + dBC4X1d6d01*u6d01dk + dBC4X1d6d00*u6d00dk + dBC4X1d7d01*u7d01dk + dBC4X1d7d00*u7d00dk + dBC4X1d8d01*u8d01dk + dBC4X1d8d00*u8d00dk + dBC4X1d9d01*u9d01dk + dBC4X1d9d00*u9d00dk;
			
			//----- BC5 -----
			BC5X1       =       BC5X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d0d00 = dBC5X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d0d10 = dBC5X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d0d01 = dBC5X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d1d00 = dBC5X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d1d10 = dBC5X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d1d01 = dBC5X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d2d00 = dBC5X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d2d10 = dBC5X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d2d01 = dBC5X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d3d00 = dBC5X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d3d10 = dBC5X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d3d01 = dBC5X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d4d00 = dBC5X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d4d10 = dBC5X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d4d01 = dBC5X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d5d00 = dBC5X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d5d10 = dBC5X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d5d01 = dBC5X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d6d00 = dBC5X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d6d10 = dBC5X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d6d01 = dBC5X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d7d00 = dBC5X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d7d10 = dBC5X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d7d01 = dBC5X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d8d00 = dBC5X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d8d10 = dBC5X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d8d01 = dBC5X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5X1d9d00 = dBC5X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d9d10 = dBC5X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5X1d9d01 = dBC5X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 5;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC5X1d5d00*a5d0[a] + dBC5X1d5d10*a5d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC5X1d0d01*b0d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC5X1);
			Dxvec[k +pn] = dBC5X1d0d10*u0d10dk + dBC5X1d0d00*u0d00dk + dBC5X1d1d10*u1d10dk + dBC5X1d1d00*u1d00dk + dBC5X1d2d10*u2d10dk + dBC5X1d2d00*u2d00dk + dBC5X1d3d10*u3d10dk + dBC5X1d3d00*u3d00dk + dBC5X1d4d10*u4d10dk + dBC5X1d4d00*u4d00dk + dBC5X1d5d10*u5d10dk + dBC5X1d5d00*u5d00dk + dBC5X1d6d10*u6d10dk + dBC5X1d6d00*u6d00dk + dBC5X1d7d10*u7d10dk + dBC5X1d7d00*u7d00dk + dBC5X1d8d10*u8d10dk + dBC5X1d8d00*u8d00dk + dBC5X1d9d10*u9d10dk + dBC5X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC5X1d0d01*u0d01dk + dBC5X1d0d00*u0d00dk + dBC5X1d1d01*u1d01dk + dBC5X1d1d00*u1d00dk + dBC5X1d2d01*u2d01dk + dBC5X1d2d00*u2d00dk + dBC5X1d3d01*u3d01dk + dBC5X1d3d00*u3d00dk + dBC5X1d4d01*u4d01dk + dBC5X1d4d00*u4d00dk + dBC5X1d5d01*u5d01dk + dBC5X1d5d00*u5d00dk + dBC5X1d6d01*u6d01dk + dBC5X1d6d00*u6d00dk + dBC5X1d7d01*u7d01dk + dBC5X1d7d00*u7d00dk + dBC5X1d8d01*u8d01dk + dBC5X1d8d00*u8d00dk + dBC5X1d9d01*u9d01dk + dBC5X1d9d00*u9d00dk;
			
			//----- BC6 -----
			BC6X1       =       BC6X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d0d00 = dBC6X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d0d10 = dBC6X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d0d01 = dBC6X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d1d00 = dBC6X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d1d10 = dBC6X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d1d01 = dBC6X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d2d00 = dBC6X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d2d10 = dBC6X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d2d01 = dBC6X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d3d00 = dBC6X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d3d10 = dBC6X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d3d01 = dBC6X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d4d00 = dBC6X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d4d10 = dBC6X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d4d01 = dBC6X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d5d00 = dBC6X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d5d10 = dBC6X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d5d01 = dBC6X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d6d00 = dBC6X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d6d10 = dBC6X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d6d01 = dBC6X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d7d00 = dBC6X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d7d10 = dBC6X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d7d01 = dBC6X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d8d00 = dBC6X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d8d10 = dBC6X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d8d01 = dBC6X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6X1d9d00 = dBC6X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d9d10 = dBC6X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6X1d9d01 = dBC6X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 6;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC6X1d6d00*a6d0[a] + dBC6X1d6d10*a6d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +6;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC6X1d1d01*b1d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC6X1);
			Dxvec[k +pn] = dBC6X1d0d10*u0d10dk + dBC6X1d0d00*u0d00dk + dBC6X1d1d10*u1d10dk + dBC6X1d1d00*u1d00dk + dBC6X1d2d10*u2d10dk + dBC6X1d2d00*u2d00dk + dBC6X1d3d10*u3d10dk + dBC6X1d3d00*u3d00dk + dBC6X1d4d10*u4d10dk + dBC6X1d4d00*u4d00dk + dBC6X1d5d10*u5d10dk + dBC6X1d5d00*u5d00dk + dBC6X1d6d10*u6d10dk + dBC6X1d6d00*u6d00dk + dBC6X1d7d10*u7d10dk + dBC6X1d7d00*u7d00dk + dBC6X1d8d10*u8d10dk + dBC6X1d8d00*u8d00dk + dBC6X1d9d10*u9d10dk + dBC6X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC6X1d0d01*u0d01dk + dBC6X1d0d00*u0d00dk + dBC6X1d1d01*u1d01dk + dBC6X1d1d00*u1d00dk + dBC6X1d2d01*u2d01dk + dBC6X1d2d00*u2d00dk + dBC6X1d3d01*u3d01dk + dBC6X1d3d00*u3d00dk + dBC6X1d4d01*u4d01dk + dBC6X1d4d00*u4d00dk + dBC6X1d5d01*u5d01dk + dBC6X1d5d00*u5d00dk + dBC6X1d6d01*u6d01dk + dBC6X1d6d00*u6d00dk + dBC6X1d7d01*u7d01dk + dBC6X1d7d00*u7d00dk + dBC6X1d8d01*u8d01dk + dBC6X1d8d00*u8d00dk + dBC6X1d9d01*u9d01dk + dBC6X1d9d00*u9d00dk;
			
			//----- BC7 -----
			BC7X1       =       BC7X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d0d00 = dBC7X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d0d10 = dBC7X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d0d01 = dBC7X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d1d00 = dBC7X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d1d10 = dBC7X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d1d01 = dBC7X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d2d00 = dBC7X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d2d10 = dBC7X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d2d01 = dBC7X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d3d00 = dBC7X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d3d10 = dBC7X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d3d01 = dBC7X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d4d00 = dBC7X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d4d10 = dBC7X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d4d01 = dBC7X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d5d00 = dBC7X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d5d10 = dBC7X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d5d01 = dBC7X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d6d00 = dBC7X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d6d10 = dBC7X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d6d01 = dBC7X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d7d00 = dBC7X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d7d10 = dBC7X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d7d01 = dBC7X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d8d00 = dBC7X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d8d10 = dBC7X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d8d01 = dBC7X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7X1d9d00 = dBC7X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d9d10 = dBC7X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7X1d9d01 = dBC7X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 7;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC7X1d7d00*a7d0[a] + dBC7X1d7d10*a7d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC7X1d2d01*b2d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC7X1);
			Dxvec[k +pn] = dBC7X1d0d10*u0d10dk + dBC7X1d0d00*u0d00dk + dBC7X1d1d10*u1d10dk + dBC7X1d1d00*u1d00dk + dBC7X1d2d10*u2d10dk + dBC7X1d2d00*u2d00dk + dBC7X1d3d10*u3d10dk + dBC7X1d3d00*u3d00dk + dBC7X1d4d10*u4d10dk + dBC7X1d4d00*u4d00dk + dBC7X1d5d10*u5d10dk + dBC7X1d5d00*u5d00dk + dBC7X1d6d10*u6d10dk + dBC7X1d6d00*u6d00dk + dBC7X1d7d10*u7d10dk + dBC7X1d7d00*u7d00dk + dBC7X1d8d10*u8d10dk + dBC7X1d8d00*u8d00dk + dBC7X1d9d10*u9d10dk + dBC7X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC7X1d0d01*u0d01dk + dBC7X1d0d00*u0d00dk + dBC7X1d1d01*u1d01dk + dBC7X1d1d00*u1d00dk + dBC7X1d2d01*u2d01dk + dBC7X1d2d00*u2d00dk + dBC7X1d3d01*u3d01dk + dBC7X1d3d00*u3d00dk + dBC7X1d4d01*u4d01dk + dBC7X1d4d00*u4d00dk + dBC7X1d5d01*u5d01dk + dBC7X1d5d00*u5d00dk + dBC7X1d6d01*u6d01dk + dBC7X1d6d00*u6d00dk + dBC7X1d7d01*u7d01dk + dBC7X1d7d00*u7d00dk + dBC7X1d8d01*u8d01dk + dBC7X1d8d00*u8d00dk + dBC7X1d9d01*u9d01dk + dBC7X1d9d00*u9d00dk;
			
			//----- BC8 -----
			BC8X1       =       BC8X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d0d00 = dBC8X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d0d10 = dBC8X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d0d01 = dBC8X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d1d00 = dBC8X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d1d10 = dBC8X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d1d01 = dBC8X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d2d00 = dBC8X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d2d10 = dBC8X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d2d01 = dBC8X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d3d00 = dBC8X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d3d10 = dBC8X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d3d01 = dBC8X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d4d00 = dBC8X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d4d10 = dBC8X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d4d01 = dBC8X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d5d00 = dBC8X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d5d10 = dBC8X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d5d01 = dBC8X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d6d00 = dBC8X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d6d10 = dBC8X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d6d01 = dBC8X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d7d00 = dBC8X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d7d10 = dBC8X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d7d01 = dBC8X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d8d00 = dBC8X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d8d10 = dBC8X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d8d01 = dBC8X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8X1d9d00 = dBC8X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d9d10 = dBC8X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8X1d9d01 = dBC8X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 8;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC8X1d8d00*a8d0[a] + dBC8X1d8d10*a8d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC8X1d3d01*b3d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC8X1);
			Dxvec[k +pn] = dBC8X1d0d10*u0d10dk + dBC8X1d0d00*u0d00dk + dBC8X1d1d10*u1d10dk + dBC8X1d1d00*u1d00dk + dBC8X1d2d10*u2d10dk + dBC8X1d2d00*u2d00dk + dBC8X1d3d10*u3d10dk + dBC8X1d3d00*u3d00dk + dBC8X1d4d10*u4d10dk + dBC8X1d4d00*u4d00dk + dBC8X1d5d10*u5d10dk + dBC8X1d5d00*u5d00dk + dBC8X1d6d10*u6d10dk + dBC8X1d6d00*u6d00dk + dBC8X1d7d10*u7d10dk + dBC8X1d7d00*u7d00dk + dBC8X1d8d10*u8d10dk + dBC8X1d8d00*u8d00dk + dBC8X1d9d10*u9d10dk + dBC8X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC8X1d0d01*u0d01dk + dBC8X1d0d00*u0d00dk + dBC8X1d1d01*u1d01dk + dBC8X1d1d00*u1d00dk + dBC8X1d2d01*u2d01dk + dBC8X1d2d00*u2d00dk + dBC8X1d3d01*u3d01dk + dBC8X1d3d00*u3d00dk + dBC8X1d4d01*u4d01dk + dBC8X1d4d00*u4d00dk + dBC8X1d5d01*u5d01dk + dBC8X1d5d00*u5d00dk + dBC8X1d6d01*u6d01dk + dBC8X1d6d00*u6d00dk + dBC8X1d7d01*u7d01dk + dBC8X1d7d00*u7d00dk + dBC8X1d8d01*u8d01dk + dBC8X1d8d00*u8d00dk + dBC8X1d9d01*u9d01dk + dBC8X1d9d00*u9d00dk;
			
			//----- BC9 -----
			BC9X1       =       BC9X1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d0d00 = dBC9X1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d0d10 = dBC9X1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d0d01 = dBC9X1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d1d00 = dBC9X1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d1d10 = dBC9X1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d1d01 = dBC9X1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d2d00 = dBC9X1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d2d10 = dBC9X1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d2d01 = dBC9X1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d3d00 = dBC9X1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d3d10 = dBC9X1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d3d01 = dBC9X1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d4d00 = dBC9X1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d4d10 = dBC9X1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d4d01 = dBC9X1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d5d00 = dBC9X1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d5d10 = dBC9X1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d5d01 = dBC9X1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d6d00 = dBC9X1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d6d10 = dBC9X1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d6d01 = dBC9X1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d7d00 = dBC9X1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d7d10 = dBC9X1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d7d01 = dBC9X1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d8d00 = dBC9X1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d8d10 = dBC9X1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d8d01 = dBC9X1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9X1d9d00 = dBC9X1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d9d10 = dBC9X1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9X1d9d01 = dBC9X1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 9;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC9X1d9d00*a9d0[a] + dBC9X1d9d10*a9d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC9X1d4d01*b4d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC9X1);
			Dxvec[k +pn] = dBC9X1d0d10*u0d10dk + dBC9X1d0d00*u0d00dk + dBC9X1d1d10*u1d10dk + dBC9X1d1d00*u1d00dk + dBC9X1d2d10*u2d10dk + dBC9X1d2d00*u2d00dk + dBC9X1d3d10*u3d10dk + dBC9X1d3d00*u3d00dk + dBC9X1d4d10*u4d10dk + dBC9X1d4d00*u4d00dk + dBC9X1d5d10*u5d10dk + dBC9X1d5d00*u5d00dk + dBC9X1d6d10*u6d10dk + dBC9X1d6d00*u6d00dk + dBC9X1d7d10*u7d10dk + dBC9X1d7d00*u7d00dk + dBC9X1d8d10*u8d10dk + dBC9X1d8d00*u8d00dk + dBC9X1d9d10*u9d10dk + dBC9X1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC9X1d0d01*u0d01dk + dBC9X1d0d00*u0d00dk + dBC9X1d1d01*u1d01dk + dBC9X1d1d00*u1d00dk + dBC9X1d2d01*u2d01dk + dBC9X1d2d00*u2d00dk + dBC9X1d3d01*u3d01dk + dBC9X1d3d00*u3d00dk + dBC9X1d4d01*u4d01dk + dBC9X1d4d00*u4d00dk + dBC9X1d5d01*u5d01dk + dBC9X1d5d00*u5d00dk + dBC9X1d6d01*u6d01dk + dBC9X1d6d00*u6d00dk + dBC9X1d7d01*u7d01dk + dBC9X1d7d00*u7d00dk + dBC9X1d8d01*u8d01dk + dBC9X1d8d00*u8d00dk + dBC9X1d9d01*u9d01dk + dBC9X1d9d00*u9d00dk;
			
		}
	}
	
	//:::::::::::::::::::: Y1 ::::::::::::::::::::
	i = 0;
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
			
			//----- BC0 -----
			BC0Y1       =       BC0Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d0d00 = dBC0Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d0d10 = dBC0Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d0d01 = dBC0Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d1d00 = dBC0Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d1d10 = dBC0Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d1d01 = dBC0Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d2d00 = dBC0Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d2d10 = dBC0Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d2d01 = dBC0Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d3d00 = dBC0Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d3d10 = dBC0Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d3d01 = dBC0Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d4d00 = dBC0Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d4d10 = dBC0Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d4d01 = dBC0Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d5d00 = dBC0Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d5d10 = dBC0Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d5d01 = dBC0Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d6d00 = dBC0Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d6d10 = dBC0Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d6d01 = dBC0Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d7d00 = dBC0Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d7d10 = dBC0Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d7d01 = dBC0Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d8d00 = dBC0Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d8d10 = dBC0Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d8d01 = dBC0Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y1d9d00 = dBC0Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d9d10 = dBC0Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y1d9d01 = dBC0Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 0;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC0Y1d0d01*b0d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC0Y1);
			Dxvec[k +pn] = dBC0Y1d0d10*u0d10dk + dBC0Y1d0d00*u0d00dk + dBC0Y1d1d10*u1d10dk + dBC0Y1d1d00*u1d00dk + dBC0Y1d2d10*u2d10dk + dBC0Y1d2d00*u2d00dk + dBC0Y1d3d10*u3d10dk + dBC0Y1d3d00*u3d00dk + dBC0Y1d4d10*u4d10dk + dBC0Y1d4d00*u4d00dk + dBC0Y1d5d10*u5d10dk + dBC0Y1d5d00*u5d00dk + dBC0Y1d6d10*u6d10dk + dBC0Y1d6d00*u6d00dk + dBC0Y1d7d10*u7d10dk + dBC0Y1d7d00*u7d00dk + dBC0Y1d8d10*u8d10dk + dBC0Y1d8d00*u8d00dk + dBC0Y1d9d10*u9d10dk + dBC0Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC0Y1d0d01*u0d01dk + dBC0Y1d0d00*u0d00dk + dBC0Y1d1d01*u1d01dk + dBC0Y1d1d00*u1d00dk + dBC0Y1d2d01*u2d01dk + dBC0Y1d2d00*u2d00dk + dBC0Y1d3d01*u3d01dk + dBC0Y1d3d00*u3d00dk + dBC0Y1d4d01*u4d01dk + dBC0Y1d4d00*u4d00dk + dBC0Y1d5d01*u5d01dk + dBC0Y1d5d00*u5d00dk + dBC0Y1d6d01*u6d01dk + dBC0Y1d6d00*u6d00dk + dBC0Y1d7d01*u7d01dk + dBC0Y1d7d00*u7d00dk + dBC0Y1d8d01*u8d01dk + dBC0Y1d8d00*u8d00dk + dBC0Y1d9d01*u9d01dk + dBC0Y1d9d00*u9d00dk;
			
			//----- BC1 -----
			BC1Y1       =       BC1Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d0d00 = dBC1Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d0d10 = dBC1Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d0d01 = dBC1Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d1d00 = dBC1Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d1d10 = dBC1Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d1d01 = dBC1Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d2d00 = dBC1Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d2d10 = dBC1Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d2d01 = dBC1Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d3d00 = dBC1Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d3d10 = dBC1Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d3d01 = dBC1Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d4d00 = dBC1Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d4d10 = dBC1Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d4d01 = dBC1Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d5d00 = dBC1Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d5d10 = dBC1Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d5d01 = dBC1Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d6d00 = dBC1Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d6d10 = dBC1Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d6d01 = dBC1Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d7d00 = dBC1Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d7d10 = dBC1Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d7d01 = dBC1Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d8d00 = dBC1Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d8d10 = dBC1Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d8d01 = dBC1Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y1d9d00 = dBC1Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d9d10 = dBC1Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y1d9d01 = dBC1Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 1;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC1Y1d1d01*b1d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC1Y1);
			Dxvec[k +pn] = dBC1Y1d0d10*u0d10dk + dBC1Y1d0d00*u0d00dk + dBC1Y1d1d10*u1d10dk + dBC1Y1d1d00*u1d00dk + dBC1Y1d2d10*u2d10dk + dBC1Y1d2d00*u2d00dk + dBC1Y1d3d10*u3d10dk + dBC1Y1d3d00*u3d00dk + dBC1Y1d4d10*u4d10dk + dBC1Y1d4d00*u4d00dk + dBC1Y1d5d10*u5d10dk + dBC1Y1d5d00*u5d00dk + dBC1Y1d6d10*u6d10dk + dBC1Y1d6d00*u6d00dk + dBC1Y1d7d10*u7d10dk + dBC1Y1d7d00*u7d00dk + dBC1Y1d8d10*u8d10dk + dBC1Y1d8d00*u8d00dk + dBC1Y1d9d10*u9d10dk + dBC1Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC1Y1d0d01*u0d01dk + dBC1Y1d0d00*u0d00dk + dBC1Y1d1d01*u1d01dk + dBC1Y1d1d00*u1d00dk + dBC1Y1d2d01*u2d01dk + dBC1Y1d2d00*u2d00dk + dBC1Y1d3d01*u3d01dk + dBC1Y1d3d00*u3d00dk + dBC1Y1d4d01*u4d01dk + dBC1Y1d4d00*u4d00dk + dBC1Y1d5d01*u5d01dk + dBC1Y1d5d00*u5d00dk + dBC1Y1d6d01*u6d01dk + dBC1Y1d6d00*u6d00dk + dBC1Y1d7d01*u7d01dk + dBC1Y1d7d00*u7d00dk + dBC1Y1d8d01*u8d01dk + dBC1Y1d8d00*u8d00dk + dBC1Y1d9d01*u9d01dk + dBC1Y1d9d00*u9d00dk;
			
			//----- BC2 -----
			BC2Y1       =       BC2Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d0d00 = dBC2Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d0d10 = dBC2Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d0d01 = dBC2Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d1d00 = dBC2Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d1d10 = dBC2Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d1d01 = dBC2Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d2d00 = dBC2Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d2d10 = dBC2Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d2d01 = dBC2Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d3d00 = dBC2Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d3d10 = dBC2Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d3d01 = dBC2Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d4d00 = dBC2Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d4d10 = dBC2Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d4d01 = dBC2Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d5d00 = dBC2Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d5d10 = dBC2Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d5d01 = dBC2Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d6d00 = dBC2Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d6d10 = dBC2Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d6d01 = dBC2Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d7d00 = dBC2Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d7d10 = dBC2Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d7d01 = dBC2Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d8d00 = dBC2Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d8d10 = dBC2Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d8d01 = dBC2Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y1d9d00 = dBC2Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d9d10 = dBC2Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y1d9d01 = dBC2Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 2;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC2Y1d2d01*b2d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC2Y1);
			Dxvec[k +pn] = dBC2Y1d0d10*u0d10dk + dBC2Y1d0d00*u0d00dk + dBC2Y1d1d10*u1d10dk + dBC2Y1d1d00*u1d00dk + dBC2Y1d2d10*u2d10dk + dBC2Y1d2d00*u2d00dk + dBC2Y1d3d10*u3d10dk + dBC2Y1d3d00*u3d00dk + dBC2Y1d4d10*u4d10dk + dBC2Y1d4d00*u4d00dk + dBC2Y1d5d10*u5d10dk + dBC2Y1d5d00*u5d00dk + dBC2Y1d6d10*u6d10dk + dBC2Y1d6d00*u6d00dk + dBC2Y1d7d10*u7d10dk + dBC2Y1d7d00*u7d00dk + dBC2Y1d8d10*u8d10dk + dBC2Y1d8d00*u8d00dk + dBC2Y1d9d10*u9d10dk + dBC2Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC2Y1d0d01*u0d01dk + dBC2Y1d0d00*u0d00dk + dBC2Y1d1d01*u1d01dk + dBC2Y1d1d00*u1d00dk + dBC2Y1d2d01*u2d01dk + dBC2Y1d2d00*u2d00dk + dBC2Y1d3d01*u3d01dk + dBC2Y1d3d00*u3d00dk + dBC2Y1d4d01*u4d01dk + dBC2Y1d4d00*u4d00dk + dBC2Y1d5d01*u5d01dk + dBC2Y1d5d00*u5d00dk + dBC2Y1d6d01*u6d01dk + dBC2Y1d6d00*u6d00dk + dBC2Y1d7d01*u7d01dk + dBC2Y1d7d00*u7d00dk + dBC2Y1d8d01*u8d01dk + dBC2Y1d8d00*u8d00dk + dBC2Y1d9d01*u9d01dk + dBC2Y1d9d00*u9d00dk;
			
			//----- BC3 -----
			BC3Y1       =       BC3Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d0d00 = dBC3Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d0d10 = dBC3Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d0d01 = dBC3Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d1d00 = dBC3Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d1d10 = dBC3Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d1d01 = dBC3Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d2d00 = dBC3Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d2d10 = dBC3Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d2d01 = dBC3Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d3d00 = dBC3Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d3d10 = dBC3Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d3d01 = dBC3Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d4d00 = dBC3Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d4d10 = dBC3Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d4d01 = dBC3Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d5d00 = dBC3Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d5d10 = dBC3Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d5d01 = dBC3Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d6d00 = dBC3Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d6d10 = dBC3Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d6d01 = dBC3Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d7d00 = dBC3Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d7d10 = dBC3Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d7d01 = dBC3Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d8d00 = dBC3Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d8d10 = dBC3Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d8d01 = dBC3Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y1d9d00 = dBC3Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d9d10 = dBC3Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y1d9d01 = dBC3Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 3;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC3Y1d3d01*b3d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC3Y1);
			Dxvec[k +pn] = dBC3Y1d0d10*u0d10dk + dBC3Y1d0d00*u0d00dk + dBC3Y1d1d10*u1d10dk + dBC3Y1d1d00*u1d00dk + dBC3Y1d2d10*u2d10dk + dBC3Y1d2d00*u2d00dk + dBC3Y1d3d10*u3d10dk + dBC3Y1d3d00*u3d00dk + dBC3Y1d4d10*u4d10dk + dBC3Y1d4d00*u4d00dk + dBC3Y1d5d10*u5d10dk + dBC3Y1d5d00*u5d00dk + dBC3Y1d6d10*u6d10dk + dBC3Y1d6d00*u6d00dk + dBC3Y1d7d10*u7d10dk + dBC3Y1d7d00*u7d00dk + dBC3Y1d8d10*u8d10dk + dBC3Y1d8d00*u8d00dk + dBC3Y1d9d10*u9d10dk + dBC3Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC3Y1d0d01*u0d01dk + dBC3Y1d0d00*u0d00dk + dBC3Y1d1d01*u1d01dk + dBC3Y1d1d00*u1d00dk + dBC3Y1d2d01*u2d01dk + dBC3Y1d2d00*u2d00dk + dBC3Y1d3d01*u3d01dk + dBC3Y1d3d00*u3d00dk + dBC3Y1d4d01*u4d01dk + dBC3Y1d4d00*u4d00dk + dBC3Y1d5d01*u5d01dk + dBC3Y1d5d00*u5d00dk + dBC3Y1d6d01*u6d01dk + dBC3Y1d6d00*u6d00dk + dBC3Y1d7d01*u7d01dk + dBC3Y1d7d00*u7d00dk + dBC3Y1d8d01*u8d01dk + dBC3Y1d8d00*u8d00dk + dBC3Y1d9d01*u9d01dk + dBC3Y1d9d00*u9d00dk;
			
			//----- BC4 -----
			BC4Y1       =       BC4Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d0d00 = dBC4Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d0d10 = dBC4Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d0d01 = dBC4Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d1d00 = dBC4Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d1d10 = dBC4Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d1d01 = dBC4Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d2d00 = dBC4Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d2d10 = dBC4Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d2d01 = dBC4Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d3d00 = dBC4Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d3d10 = dBC4Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d3d01 = dBC4Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d4d00 = dBC4Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d4d10 = dBC4Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d4d01 = dBC4Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d5d00 = dBC4Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d5d10 = dBC4Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d5d01 = dBC4Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d6d00 = dBC4Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d6d10 = dBC4Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d6d01 = dBC4Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d7d00 = dBC4Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d7d10 = dBC4Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d7d01 = dBC4Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d8d00 = dBC4Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d8d10 = dBC4Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d8d01 = dBC4Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y1d9d00 = dBC4Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d9d10 = dBC4Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y1d9d01 = dBC4Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 4;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC4Y1d4d01*b4d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC4Y1);
			Dxvec[k +pn] = dBC4Y1d0d10*u0d10dk + dBC4Y1d0d00*u0d00dk + dBC4Y1d1d10*u1d10dk + dBC4Y1d1d00*u1d00dk + dBC4Y1d2d10*u2d10dk + dBC4Y1d2d00*u2d00dk + dBC4Y1d3d10*u3d10dk + dBC4Y1d3d00*u3d00dk + dBC4Y1d4d10*u4d10dk + dBC4Y1d4d00*u4d00dk + dBC4Y1d5d10*u5d10dk + dBC4Y1d5d00*u5d00dk + dBC4Y1d6d10*u6d10dk + dBC4Y1d6d00*u6d00dk + dBC4Y1d7d10*u7d10dk + dBC4Y1d7d00*u7d00dk + dBC4Y1d8d10*u8d10dk + dBC4Y1d8d00*u8d00dk + dBC4Y1d9d10*u9d10dk + dBC4Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC4Y1d0d01*u0d01dk + dBC4Y1d0d00*u0d00dk + dBC4Y1d1d01*u1d01dk + dBC4Y1d1d00*u1d00dk + dBC4Y1d2d01*u2d01dk + dBC4Y1d2d00*u2d00dk + dBC4Y1d3d01*u3d01dk + dBC4Y1d3d00*u3d00dk + dBC4Y1d4d01*u4d01dk + dBC4Y1d4d00*u4d00dk + dBC4Y1d5d01*u5d01dk + dBC4Y1d5d00*u5d00dk + dBC4Y1d6d01*u6d01dk + dBC4Y1d6d00*u6d00dk + dBC4Y1d7d01*u7d01dk + dBC4Y1d7d00*u7d00dk + dBC4Y1d8d01*u8d01dk + dBC4Y1d8d00*u8d00dk + dBC4Y1d9d01*u9d01dk + dBC4Y1d9d00*u9d00dk;
			
			//----- BC5 -----
			BC5Y1       =       BC5Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d0d00 = dBC5Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d0d10 = dBC5Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d0d01 = dBC5Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d1d00 = dBC5Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d1d10 = dBC5Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d1d01 = dBC5Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d2d00 = dBC5Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d2d10 = dBC5Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d2d01 = dBC5Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d3d00 = dBC5Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d3d10 = dBC5Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d3d01 = dBC5Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d4d00 = dBC5Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d4d10 = dBC5Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d4d01 = dBC5Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d5d00 = dBC5Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d5d10 = dBC5Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d5d01 = dBC5Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d6d00 = dBC5Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d6d10 = dBC5Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d6d01 = dBC5Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d7d00 = dBC5Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d7d10 = dBC5Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d7d01 = dBC5Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d8d00 = dBC5Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d8d10 = dBC5Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d8d01 = dBC5Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y1d9d00 = dBC5Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d9d10 = dBC5Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y1d9d01 = dBC5Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 5;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC5Y1d5d00*a5d0[a] + dBC5Y1d5d10*a5d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC5Y1d0d01*b0d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC5Y1);
			Dxvec[k +pn] = dBC5Y1d0d10*u0d10dk + dBC5Y1d0d00*u0d00dk + dBC5Y1d1d10*u1d10dk + dBC5Y1d1d00*u1d00dk + dBC5Y1d2d10*u2d10dk + dBC5Y1d2d00*u2d00dk + dBC5Y1d3d10*u3d10dk + dBC5Y1d3d00*u3d00dk + dBC5Y1d4d10*u4d10dk + dBC5Y1d4d00*u4d00dk + dBC5Y1d5d10*u5d10dk + dBC5Y1d5d00*u5d00dk + dBC5Y1d6d10*u6d10dk + dBC5Y1d6d00*u6d00dk + dBC5Y1d7d10*u7d10dk + dBC5Y1d7d00*u7d00dk + dBC5Y1d8d10*u8d10dk + dBC5Y1d8d00*u8d00dk + dBC5Y1d9d10*u9d10dk + dBC5Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC5Y1d0d01*u0d01dk + dBC5Y1d0d00*u0d00dk + dBC5Y1d1d01*u1d01dk + dBC5Y1d1d00*u1d00dk + dBC5Y1d2d01*u2d01dk + dBC5Y1d2d00*u2d00dk + dBC5Y1d3d01*u3d01dk + dBC5Y1d3d00*u3d00dk + dBC5Y1d4d01*u4d01dk + dBC5Y1d4d00*u4d00dk + dBC5Y1d5d01*u5d01dk + dBC5Y1d5d00*u5d00dk + dBC5Y1d6d01*u6d01dk + dBC5Y1d6d00*u6d00dk + dBC5Y1d7d01*u7d01dk + dBC5Y1d7d00*u7d00dk + dBC5Y1d8d01*u8d01dk + dBC5Y1d8d00*u8d00dk + dBC5Y1d9d01*u9d01dk + dBC5Y1d9d00*u9d00dk;
			
			//----- BC6 -----
			BC6Y1       =       BC6Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d0d00 = dBC6Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d0d10 = dBC6Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d0d01 = dBC6Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d1d00 = dBC6Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d1d10 = dBC6Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d1d01 = dBC6Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d2d00 = dBC6Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d2d10 = dBC6Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d2d01 = dBC6Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d3d00 = dBC6Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d3d10 = dBC6Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d3d01 = dBC6Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d4d00 = dBC6Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d4d10 = dBC6Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d4d01 = dBC6Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d5d00 = dBC6Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d5d10 = dBC6Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d5d01 = dBC6Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d6d00 = dBC6Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d6d10 = dBC6Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d6d01 = dBC6Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d7d00 = dBC6Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d7d10 = dBC6Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d7d01 = dBC6Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d8d00 = dBC6Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d8d10 = dBC6Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d8d01 = dBC6Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y1d9d00 = dBC6Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d9d10 = dBC6Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y1d9d01 = dBC6Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 6;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC6Y1d6d00*a6d0[a] + dBC6Y1d6d10*a6d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +6;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC6Y1d1d01*b1d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC6Y1);
			Dxvec[k +pn] = dBC6Y1d0d10*u0d10dk + dBC6Y1d0d00*u0d00dk + dBC6Y1d1d10*u1d10dk + dBC6Y1d1d00*u1d00dk + dBC6Y1d2d10*u2d10dk + dBC6Y1d2d00*u2d00dk + dBC6Y1d3d10*u3d10dk + dBC6Y1d3d00*u3d00dk + dBC6Y1d4d10*u4d10dk + dBC6Y1d4d00*u4d00dk + dBC6Y1d5d10*u5d10dk + dBC6Y1d5d00*u5d00dk + dBC6Y1d6d10*u6d10dk + dBC6Y1d6d00*u6d00dk + dBC6Y1d7d10*u7d10dk + dBC6Y1d7d00*u7d00dk + dBC6Y1d8d10*u8d10dk + dBC6Y1d8d00*u8d00dk + dBC6Y1d9d10*u9d10dk + dBC6Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC6Y1d0d01*u0d01dk + dBC6Y1d0d00*u0d00dk + dBC6Y1d1d01*u1d01dk + dBC6Y1d1d00*u1d00dk + dBC6Y1d2d01*u2d01dk + dBC6Y1d2d00*u2d00dk + dBC6Y1d3d01*u3d01dk + dBC6Y1d3d00*u3d00dk + dBC6Y1d4d01*u4d01dk + dBC6Y1d4d00*u4d00dk + dBC6Y1d5d01*u5d01dk + dBC6Y1d5d00*u5d00dk + dBC6Y1d6d01*u6d01dk + dBC6Y1d6d00*u6d00dk + dBC6Y1d7d01*u7d01dk + dBC6Y1d7d00*u7d00dk + dBC6Y1d8d01*u8d01dk + dBC6Y1d8d00*u8d00dk + dBC6Y1d9d01*u9d01dk + dBC6Y1d9d00*u9d00dk;
			
			//----- BC7 -----
			BC7Y1       =       BC7Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d0d00 = dBC7Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d0d10 = dBC7Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d0d01 = dBC7Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d1d00 = dBC7Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d1d10 = dBC7Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d1d01 = dBC7Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d2d00 = dBC7Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d2d10 = dBC7Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d2d01 = dBC7Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d3d00 = dBC7Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d3d10 = dBC7Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d3d01 = dBC7Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d4d00 = dBC7Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d4d10 = dBC7Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d4d01 = dBC7Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d5d00 = dBC7Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d5d10 = dBC7Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d5d01 = dBC7Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d6d00 = dBC7Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d6d10 = dBC7Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d6d01 = dBC7Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d7d00 = dBC7Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d7d10 = dBC7Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d7d01 = dBC7Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d8d00 = dBC7Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d8d10 = dBC7Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d8d01 = dBC7Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y1d9d00 = dBC7Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d9d10 = dBC7Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y1d9d01 = dBC7Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 7;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC7Y1d7d00*a7d0[a] + dBC7Y1d7d10*a7d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC7Y1d2d01*b2d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC7Y1);
			Dxvec[k +pn] = dBC7Y1d0d10*u0d10dk + dBC7Y1d0d00*u0d00dk + dBC7Y1d1d10*u1d10dk + dBC7Y1d1d00*u1d00dk + dBC7Y1d2d10*u2d10dk + dBC7Y1d2d00*u2d00dk + dBC7Y1d3d10*u3d10dk + dBC7Y1d3d00*u3d00dk + dBC7Y1d4d10*u4d10dk + dBC7Y1d4d00*u4d00dk + dBC7Y1d5d10*u5d10dk + dBC7Y1d5d00*u5d00dk + dBC7Y1d6d10*u6d10dk + dBC7Y1d6d00*u6d00dk + dBC7Y1d7d10*u7d10dk + dBC7Y1d7d00*u7d00dk + dBC7Y1d8d10*u8d10dk + dBC7Y1d8d00*u8d00dk + dBC7Y1d9d10*u9d10dk + dBC7Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC7Y1d0d01*u0d01dk + dBC7Y1d0d00*u0d00dk + dBC7Y1d1d01*u1d01dk + dBC7Y1d1d00*u1d00dk + dBC7Y1d2d01*u2d01dk + dBC7Y1d2d00*u2d00dk + dBC7Y1d3d01*u3d01dk + dBC7Y1d3d00*u3d00dk + dBC7Y1d4d01*u4d01dk + dBC7Y1d4d00*u4d00dk + dBC7Y1d5d01*u5d01dk + dBC7Y1d5d00*u5d00dk + dBC7Y1d6d01*u6d01dk + dBC7Y1d6d00*u6d00dk + dBC7Y1d7d01*u7d01dk + dBC7Y1d7d00*u7d00dk + dBC7Y1d8d01*u8d01dk + dBC7Y1d8d00*u8d00dk + dBC7Y1d9d01*u9d01dk + dBC7Y1d9d00*u9d00dk;
			
			//----- BC8 -----
			BC8Y1       =       BC8Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d0d00 = dBC8Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d0d10 = dBC8Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d0d01 = dBC8Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d1d00 = dBC8Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d1d10 = dBC8Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d1d01 = dBC8Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d2d00 = dBC8Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d2d10 = dBC8Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d2d01 = dBC8Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d3d00 = dBC8Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d3d10 = dBC8Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d3d01 = dBC8Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d4d00 = dBC8Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d4d10 = dBC8Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d4d01 = dBC8Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d5d00 = dBC8Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d5d10 = dBC8Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d5d01 = dBC8Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d6d00 = dBC8Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d6d10 = dBC8Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d6d01 = dBC8Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d7d00 = dBC8Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d7d10 = dBC8Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d7d01 = dBC8Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d8d00 = dBC8Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d8d10 = dBC8Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d8d01 = dBC8Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y1d9d00 = dBC8Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d9d10 = dBC8Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y1d9d01 = dBC8Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 8;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC8Y1d8d00*a8d0[a] + dBC8Y1d8d10*a8d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC8Y1d3d01*b3d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC8Y1);
			Dxvec[k +pn] = dBC8Y1d0d10*u0d10dk + dBC8Y1d0d00*u0d00dk + dBC8Y1d1d10*u1d10dk + dBC8Y1d1d00*u1d00dk + dBC8Y1d2d10*u2d10dk + dBC8Y1d2d00*u2d00dk + dBC8Y1d3d10*u3d10dk + dBC8Y1d3d00*u3d00dk + dBC8Y1d4d10*u4d10dk + dBC8Y1d4d00*u4d00dk + dBC8Y1d5d10*u5d10dk + dBC8Y1d5d00*u5d00dk + dBC8Y1d6d10*u6d10dk + dBC8Y1d6d00*u6d00dk + dBC8Y1d7d10*u7d10dk + dBC8Y1d7d00*u7d00dk + dBC8Y1d8d10*u8d10dk + dBC8Y1d8d00*u8d00dk + dBC8Y1d9d10*u9d10dk + dBC8Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC8Y1d0d01*u0d01dk + dBC8Y1d0d00*u0d00dk + dBC8Y1d1d01*u1d01dk + dBC8Y1d1d00*u1d00dk + dBC8Y1d2d01*u2d01dk + dBC8Y1d2d00*u2d00dk + dBC8Y1d3d01*u3d01dk + dBC8Y1d3d00*u3d00dk + dBC8Y1d4d01*u4d01dk + dBC8Y1d4d00*u4d00dk + dBC8Y1d5d01*u5d01dk + dBC8Y1d5d00*u5d00dk + dBC8Y1d6d01*u6d01dk + dBC8Y1d6d00*u6d00dk + dBC8Y1d7d01*u7d01dk + dBC8Y1d7d00*u7d00dk + dBC8Y1d8d01*u8d01dk + dBC8Y1d8d00*u8d00dk + dBC8Y1d9d01*u9d01dk + dBC8Y1d9d00*u9d00dk;
			
			//----- BC9 -----
			BC9Y1       =       BC9Y1out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d0d00 = dBC9Y1d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d0d10 = dBC9Y1d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d0d01 = dBC9Y1d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d1d00 = dBC9Y1d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d1d10 = dBC9Y1d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d1d01 = dBC9Y1d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d2d00 = dBC9Y1d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d2d10 = dBC9Y1d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d2d01 = dBC9Y1d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d3d00 = dBC9Y1d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d3d10 = dBC9Y1d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d3d01 = dBC9Y1d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d4d00 = dBC9Y1d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d4d10 = dBC9Y1d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d4d01 = dBC9Y1d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d5d00 = dBC9Y1d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d5d10 = dBC9Y1d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d5d01 = dBC9Y1d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d6d00 = dBC9Y1d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d6d10 = dBC9Y1d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d6d01 = dBC9Y1d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d7d00 = dBC9Y1d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d7d10 = dBC9Y1d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d7d01 = dBC9Y1d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d8d00 = dBC9Y1d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d8d10 = dBC9Y1d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d8d01 = dBC9Y1d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y1d9d00 = dBC9Y1d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d9d10 = dBC9Y1d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y1d9d01 = dBC9Y1d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 9;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC9Y1d9d00*a9d0[a] + dBC9Y1d9d10*a9d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC9Y1d4d01*b4d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC9Y1);
			Dxvec[k +pn] = dBC9Y1d0d10*u0d10dk + dBC9Y1d0d00*u0d00dk + dBC9Y1d1d10*u1d10dk + dBC9Y1d1d00*u1d00dk + dBC9Y1d2d10*u2d10dk + dBC9Y1d2d00*u2d00dk + dBC9Y1d3d10*u3d10dk + dBC9Y1d3d00*u3d00dk + dBC9Y1d4d10*u4d10dk + dBC9Y1d4d00*u4d00dk + dBC9Y1d5d10*u5d10dk + dBC9Y1d5d00*u5d00dk + dBC9Y1d6d10*u6d10dk + dBC9Y1d6d00*u6d00dk + dBC9Y1d7d10*u7d10dk + dBC9Y1d7d00*u7d00dk + dBC9Y1d8d10*u8d10dk + dBC9Y1d8d00*u8d00dk + dBC9Y1d9d10*u9d10dk + dBC9Y1d9d00*u9d00dk;
			Dyvec[k +pn] = dBC9Y1d0d01*u0d01dk + dBC9Y1d0d00*u0d00dk + dBC9Y1d1d01*u1d01dk + dBC9Y1d1d00*u1d00dk + dBC9Y1d2d01*u2d01dk + dBC9Y1d2d00*u2d00dk + dBC9Y1d3d01*u3d01dk + dBC9Y1d3d00*u3d00dk + dBC9Y1d4d01*u4d01dk + dBC9Y1d4d00*u4d00dk + dBC9Y1d5d01*u5d01dk + dBC9Y1d5d00*u5d00dk + dBC9Y1d6d01*u6d01dk + dBC9Y1d6d00*u6d00dk + dBC9Y1d7d01*u7d01dk + dBC9Y1d7d00*u7d00dk + dBC9Y1d8d01*u8d01dk + dBC9Y1d8d00*u8d00dk + dBC9Y1d9d01*u9d01dk + dBC9Y1d9d00*u9d00dk;
			
		}
	}
	
	//:::::::::::::::::::: Y0 ::::::::::::::::::::
	i = m-1;
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
			
			//----- BC0 -----
			BC0Y0       =       BC0Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d0d00 = dBC0Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d0d10 = dBC0Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d0d01 = dBC0Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d1d00 = dBC0Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d1d10 = dBC0Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d1d01 = dBC0Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d2d00 = dBC0Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d2d10 = dBC0Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d2d01 = dBC0Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d3d00 = dBC0Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d3d10 = dBC0Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d3d01 = dBC0Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d4d00 = dBC0Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d4d10 = dBC0Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d4d01 = dBC0Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d5d00 = dBC0Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d5d10 = dBC0Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d5d01 = dBC0Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d6d00 = dBC0Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d6d10 = dBC0Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d6d01 = dBC0Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d7d00 = dBC0Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d7d10 = dBC0Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d7d01 = dBC0Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d8d00 = dBC0Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d8d10 = dBC0Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d8d01 = dBC0Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC0Y0d9d00 = dBC0Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d9d10 = dBC0Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC0Y0d9d01 = dBC0Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 0;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC0Y0d0d01*b0d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC0Y0);
			Dxvec[k +pn] = dBC0Y0d0d10*u0d10dk + dBC0Y0d0d00*u0d00dk + dBC0Y0d1d10*u1d10dk + dBC0Y0d1d00*u1d00dk + dBC0Y0d2d10*u2d10dk + dBC0Y0d2d00*u2d00dk + dBC0Y0d3d10*u3d10dk + dBC0Y0d3d00*u3d00dk + dBC0Y0d4d10*u4d10dk + dBC0Y0d4d00*u4d00dk + dBC0Y0d5d10*u5d10dk + dBC0Y0d5d00*u5d00dk + dBC0Y0d6d10*u6d10dk + dBC0Y0d6d00*u6d00dk + dBC0Y0d7d10*u7d10dk + dBC0Y0d7d00*u7d00dk + dBC0Y0d8d10*u8d10dk + dBC0Y0d8d00*u8d00dk + dBC0Y0d9d10*u9d10dk + dBC0Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC0Y0d0d01*u0d01dk + dBC0Y0d0d00*u0d00dk + dBC0Y0d1d01*u1d01dk + dBC0Y0d1d00*u1d00dk + dBC0Y0d2d01*u2d01dk + dBC0Y0d2d00*u2d00dk + dBC0Y0d3d01*u3d01dk + dBC0Y0d3d00*u3d00dk + dBC0Y0d4d01*u4d01dk + dBC0Y0d4d00*u4d00dk + dBC0Y0d5d01*u5d01dk + dBC0Y0d5d00*u5d00dk + dBC0Y0d6d01*u6d01dk + dBC0Y0d6d00*u6d00dk + dBC0Y0d7d01*u7d01dk + dBC0Y0d7d00*u7d00dk + dBC0Y0d8d01*u8d01dk + dBC0Y0d8d00*u8d00dk + dBC0Y0d9d01*u9d01dk + dBC0Y0d9d00*u9d00dk;
			
			//----- BC1 -----
			BC1Y0       =       BC1Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d0d00 = dBC1Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d0d10 = dBC1Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d0d01 = dBC1Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d1d00 = dBC1Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d1d10 = dBC1Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d1d01 = dBC1Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d2d00 = dBC1Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d2d10 = dBC1Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d2d01 = dBC1Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d3d00 = dBC1Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d3d10 = dBC1Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d3d01 = dBC1Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d4d00 = dBC1Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d4d10 = dBC1Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d4d01 = dBC1Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d5d00 = dBC1Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d5d10 = dBC1Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d5d01 = dBC1Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d6d00 = dBC1Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d6d10 = dBC1Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d6d01 = dBC1Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d7d00 = dBC1Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d7d10 = dBC1Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d7d01 = dBC1Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d8d00 = dBC1Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d8d10 = dBC1Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d8d01 = dBC1Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC1Y0d9d00 = dBC1Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d9d10 = dBC1Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC1Y0d9d01 = dBC1Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 1;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC1Y0d1d01*b1d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC1Y0);
			Dxvec[k +pn] = dBC1Y0d0d10*u0d10dk + dBC1Y0d0d00*u0d00dk + dBC1Y0d1d10*u1d10dk + dBC1Y0d1d00*u1d00dk + dBC1Y0d2d10*u2d10dk + dBC1Y0d2d00*u2d00dk + dBC1Y0d3d10*u3d10dk + dBC1Y0d3d00*u3d00dk + dBC1Y0d4d10*u4d10dk + dBC1Y0d4d00*u4d00dk + dBC1Y0d5d10*u5d10dk + dBC1Y0d5d00*u5d00dk + dBC1Y0d6d10*u6d10dk + dBC1Y0d6d00*u6d00dk + dBC1Y0d7d10*u7d10dk + dBC1Y0d7d00*u7d00dk + dBC1Y0d8d10*u8d10dk + dBC1Y0d8d00*u8d00dk + dBC1Y0d9d10*u9d10dk + dBC1Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC1Y0d0d01*u0d01dk + dBC1Y0d0d00*u0d00dk + dBC1Y0d1d01*u1d01dk + dBC1Y0d1d00*u1d00dk + dBC1Y0d2d01*u2d01dk + dBC1Y0d2d00*u2d00dk + dBC1Y0d3d01*u3d01dk + dBC1Y0d3d00*u3d00dk + dBC1Y0d4d01*u4d01dk + dBC1Y0d4d00*u4d00dk + dBC1Y0d5d01*u5d01dk + dBC1Y0d5d00*u5d00dk + dBC1Y0d6d01*u6d01dk + dBC1Y0d6d00*u6d00dk + dBC1Y0d7d01*u7d01dk + dBC1Y0d7d00*u7d00dk + dBC1Y0d8d01*u8d01dk + dBC1Y0d8d00*u8d00dk + dBC1Y0d9d01*u9d01dk + dBC1Y0d9d00*u9d00dk;
			
			//----- BC2 -----
			BC2Y0       =       BC2Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d0d00 = dBC2Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d0d10 = dBC2Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d0d01 = dBC2Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d1d00 = dBC2Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d1d10 = dBC2Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d1d01 = dBC2Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d2d00 = dBC2Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d2d10 = dBC2Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d2d01 = dBC2Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d3d00 = dBC2Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d3d10 = dBC2Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d3d01 = dBC2Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d4d00 = dBC2Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d4d10 = dBC2Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d4d01 = dBC2Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d5d00 = dBC2Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d5d10 = dBC2Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d5d01 = dBC2Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d6d00 = dBC2Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d6d10 = dBC2Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d6d01 = dBC2Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d7d00 = dBC2Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d7d10 = dBC2Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d7d01 = dBC2Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d8d00 = dBC2Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d8d10 = dBC2Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d8d01 = dBC2Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC2Y0d9d00 = dBC2Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d9d10 = dBC2Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC2Y0d9d01 = dBC2Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 2;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC2Y0d2d01*b2d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC2Y0);
			Dxvec[k +pn] = dBC2Y0d0d10*u0d10dk + dBC2Y0d0d00*u0d00dk + dBC2Y0d1d10*u1d10dk + dBC2Y0d1d00*u1d00dk + dBC2Y0d2d10*u2d10dk + dBC2Y0d2d00*u2d00dk + dBC2Y0d3d10*u3d10dk + dBC2Y0d3d00*u3d00dk + dBC2Y0d4d10*u4d10dk + dBC2Y0d4d00*u4d00dk + dBC2Y0d5d10*u5d10dk + dBC2Y0d5d00*u5d00dk + dBC2Y0d6d10*u6d10dk + dBC2Y0d6d00*u6d00dk + dBC2Y0d7d10*u7d10dk + dBC2Y0d7d00*u7d00dk + dBC2Y0d8d10*u8d10dk + dBC2Y0d8d00*u8d00dk + dBC2Y0d9d10*u9d10dk + dBC2Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC2Y0d0d01*u0d01dk + dBC2Y0d0d00*u0d00dk + dBC2Y0d1d01*u1d01dk + dBC2Y0d1d00*u1d00dk + dBC2Y0d2d01*u2d01dk + dBC2Y0d2d00*u2d00dk + dBC2Y0d3d01*u3d01dk + dBC2Y0d3d00*u3d00dk + dBC2Y0d4d01*u4d01dk + dBC2Y0d4d00*u4d00dk + dBC2Y0d5d01*u5d01dk + dBC2Y0d5d00*u5d00dk + dBC2Y0d6d01*u6d01dk + dBC2Y0d6d00*u6d00dk + dBC2Y0d7d01*u7d01dk + dBC2Y0d7d00*u7d00dk + dBC2Y0d8d01*u8d01dk + dBC2Y0d8d00*u8d00dk + dBC2Y0d9d01*u9d01dk + dBC2Y0d9d00*u9d00dk;
			
			//----- BC3 -----
			BC3Y0       =       BC3Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d0d00 = dBC3Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d0d10 = dBC3Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d0d01 = dBC3Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d1d00 = dBC3Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d1d10 = dBC3Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d1d01 = dBC3Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d2d00 = dBC3Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d2d10 = dBC3Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d2d01 = dBC3Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d3d00 = dBC3Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d3d10 = dBC3Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d3d01 = dBC3Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d4d00 = dBC3Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d4d10 = dBC3Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d4d01 = dBC3Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d5d00 = dBC3Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d5d10 = dBC3Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d5d01 = dBC3Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d6d00 = dBC3Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d6d10 = dBC3Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d6d01 = dBC3Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d7d00 = dBC3Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d7d10 = dBC3Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d7d01 = dBC3Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d8d00 = dBC3Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d8d10 = dBC3Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d8d01 = dBC3Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC3Y0d9d00 = dBC3Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d9d10 = dBC3Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC3Y0d9d01 = dBC3Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 3;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC3Y0d3d01*b3d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC3Y0);
			Dxvec[k +pn] = dBC3Y0d0d10*u0d10dk + dBC3Y0d0d00*u0d00dk + dBC3Y0d1d10*u1d10dk + dBC3Y0d1d00*u1d00dk + dBC3Y0d2d10*u2d10dk + dBC3Y0d2d00*u2d00dk + dBC3Y0d3d10*u3d10dk + dBC3Y0d3d00*u3d00dk + dBC3Y0d4d10*u4d10dk + dBC3Y0d4d00*u4d00dk + dBC3Y0d5d10*u5d10dk + dBC3Y0d5d00*u5d00dk + dBC3Y0d6d10*u6d10dk + dBC3Y0d6d00*u6d00dk + dBC3Y0d7d10*u7d10dk + dBC3Y0d7d00*u7d00dk + dBC3Y0d8d10*u8d10dk + dBC3Y0d8d00*u8d00dk + dBC3Y0d9d10*u9d10dk + dBC3Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC3Y0d0d01*u0d01dk + dBC3Y0d0d00*u0d00dk + dBC3Y0d1d01*u1d01dk + dBC3Y0d1d00*u1d00dk + dBC3Y0d2d01*u2d01dk + dBC3Y0d2d00*u2d00dk + dBC3Y0d3d01*u3d01dk + dBC3Y0d3d00*u3d00dk + dBC3Y0d4d01*u4d01dk + dBC3Y0d4d00*u4d00dk + dBC3Y0d5d01*u5d01dk + dBC3Y0d5d00*u5d00dk + dBC3Y0d6d01*u6d01dk + dBC3Y0d6d00*u6d00dk + dBC3Y0d7d01*u7d01dk + dBC3Y0d7d00*u7d00dk + dBC3Y0d8d01*u8d01dk + dBC3Y0d8d00*u8d00dk + dBC3Y0d9d01*u9d01dk + dBC3Y0d9d00*u9d00dk;
			
			//----- BC4 -----
			BC4Y0       =       BC4Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d0d00 = dBC4Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d0d10 = dBC4Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d0d01 = dBC4Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d1d00 = dBC4Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d1d10 = dBC4Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d1d01 = dBC4Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d2d00 = dBC4Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d2d10 = dBC4Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d2d01 = dBC4Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d3d00 = dBC4Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d3d10 = dBC4Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d3d01 = dBC4Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d4d00 = dBC4Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d4d10 = dBC4Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d4d01 = dBC4Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d5d00 = dBC4Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d5d10 = dBC4Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d5d01 = dBC4Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d6d00 = dBC4Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d6d10 = dBC4Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d6d01 = dBC4Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d7d00 = dBC4Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d7d10 = dBC4Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d7d01 = dBC4Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d8d00 = dBC4Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d8d10 = dBC4Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d8d01 = dBC4Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC4Y0d9d00 = dBC4Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d9d10 = dBC4Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC4Y0d9d01 = dBC4Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 4;
			
			for (a = 0; a <= r+onexside; ++a)
			{
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC4Y0d4d01*b4d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC4Y0);
			Dxvec[k +pn] = dBC4Y0d0d10*u0d10dk + dBC4Y0d0d00*u0d00dk + dBC4Y0d1d10*u1d10dk + dBC4Y0d1d00*u1d00dk + dBC4Y0d2d10*u2d10dk + dBC4Y0d2d00*u2d00dk + dBC4Y0d3d10*u3d10dk + dBC4Y0d3d00*u3d00dk + dBC4Y0d4d10*u4d10dk + dBC4Y0d4d00*u4d00dk + dBC4Y0d5d10*u5d10dk + dBC4Y0d5d00*u5d00dk + dBC4Y0d6d10*u6d10dk + dBC4Y0d6d00*u6d00dk + dBC4Y0d7d10*u7d10dk + dBC4Y0d7d00*u7d00dk + dBC4Y0d8d10*u8d10dk + dBC4Y0d8d00*u8d00dk + dBC4Y0d9d10*u9d10dk + dBC4Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC4Y0d0d01*u0d01dk + dBC4Y0d0d00*u0d00dk + dBC4Y0d1d01*u1d01dk + dBC4Y0d1d00*u1d00dk + dBC4Y0d2d01*u2d01dk + dBC4Y0d2d00*u2d00dk + dBC4Y0d3d01*u3d01dk + dBC4Y0d3d00*u3d00dk + dBC4Y0d4d01*u4d01dk + dBC4Y0d4d00*u4d00dk + dBC4Y0d5d01*u5d01dk + dBC4Y0d5d00*u5d00dk + dBC4Y0d6d01*u6d01dk + dBC4Y0d6d00*u6d00dk + dBC4Y0d7d01*u7d01dk + dBC4Y0d7d00*u7d00dk + dBC4Y0d8d01*u8d01dk + dBC4Y0d8d00*u8d00dk + dBC4Y0d9d01*u9d01dk + dBC4Y0d9d00*u9d00dk;
			
			//----- BC5 -----
			BC5Y0       =       BC5Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d0d00 = dBC5Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d0d10 = dBC5Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d0d01 = dBC5Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d1d00 = dBC5Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d1d10 = dBC5Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d1d01 = dBC5Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d2d00 = dBC5Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d2d10 = dBC5Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d2d01 = dBC5Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d3d00 = dBC5Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d3d10 = dBC5Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d3d01 = dBC5Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d4d00 = dBC5Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d4d10 = dBC5Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d4d01 = dBC5Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d5d00 = dBC5Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d5d10 = dBC5Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d5d01 = dBC5Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d6d00 = dBC5Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d6d10 = dBC5Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d6d01 = dBC5Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d7d00 = dBC5Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d7d10 = dBC5Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d7d01 = dBC5Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d8d00 = dBC5Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d8d10 = dBC5Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d8d01 = dBC5Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC5Y0d9d00 = dBC5Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d9d10 = dBC5Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC5Y0d9d01 = dBC5Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 5;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC5Y0d5d00*a5d0[a] + dBC5Y0d5d10*a5d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +5;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC5Y0d0d01*b0d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +0;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC5Y0);
			Dxvec[k +pn] = dBC5Y0d0d10*u0d10dk + dBC5Y0d0d00*u0d00dk + dBC5Y0d1d10*u1d10dk + dBC5Y0d1d00*u1d00dk + dBC5Y0d2d10*u2d10dk + dBC5Y0d2d00*u2d00dk + dBC5Y0d3d10*u3d10dk + dBC5Y0d3d00*u3d00dk + dBC5Y0d4d10*u4d10dk + dBC5Y0d4d00*u4d00dk + dBC5Y0d5d10*u5d10dk + dBC5Y0d5d00*u5d00dk + dBC5Y0d6d10*u6d10dk + dBC5Y0d6d00*u6d00dk + dBC5Y0d7d10*u7d10dk + dBC5Y0d7d00*u7d00dk + dBC5Y0d8d10*u8d10dk + dBC5Y0d8d00*u8d00dk + dBC5Y0d9d10*u9d10dk + dBC5Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC5Y0d0d01*u0d01dk + dBC5Y0d0d00*u0d00dk + dBC5Y0d1d01*u1d01dk + dBC5Y0d1d00*u1d00dk + dBC5Y0d2d01*u2d01dk + dBC5Y0d2d00*u2d00dk + dBC5Y0d3d01*u3d01dk + dBC5Y0d3d00*u3d00dk + dBC5Y0d4d01*u4d01dk + dBC5Y0d4d00*u4d00dk + dBC5Y0d5d01*u5d01dk + dBC5Y0d5d00*u5d00dk + dBC5Y0d6d01*u6d01dk + dBC5Y0d6d00*u6d00dk + dBC5Y0d7d01*u7d01dk + dBC5Y0d7d00*u7d00dk + dBC5Y0d8d01*u8d01dk + dBC5Y0d8d00*u8d00dk + dBC5Y0d9d01*u9d01dk + dBC5Y0d9d00*u9d00dk;
			
			//----- BC6 -----
			BC6Y0       =       BC6Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d0d00 = dBC6Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d0d10 = dBC6Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d0d01 = dBC6Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d1d00 = dBC6Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d1d10 = dBC6Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d1d01 = dBC6Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d2d00 = dBC6Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d2d10 = dBC6Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d2d01 = dBC6Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d3d00 = dBC6Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d3d10 = dBC6Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d3d01 = dBC6Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d4d00 = dBC6Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d4d10 = dBC6Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d4d01 = dBC6Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d5d00 = dBC6Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d5d10 = dBC6Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d5d01 = dBC6Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d6d00 = dBC6Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d6d10 = dBC6Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d6d01 = dBC6Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d7d00 = dBC6Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d7d10 = dBC6Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d7d01 = dBC6Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d8d00 = dBC6Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d8d10 = dBC6Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d8d01 = dBC6Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC6Y0d9d00 = dBC6Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d9d10 = dBC6Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC6Y0d9d01 = dBC6Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 6;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC6Y0d6d00*a6d0[a] + dBC6Y0d6d10*a6d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +6;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC6Y0d1d01*b1d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +1;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC6Y0);
			Dxvec[k +pn] = dBC6Y0d0d10*u0d10dk + dBC6Y0d0d00*u0d00dk + dBC6Y0d1d10*u1d10dk + dBC6Y0d1d00*u1d00dk + dBC6Y0d2d10*u2d10dk + dBC6Y0d2d00*u2d00dk + dBC6Y0d3d10*u3d10dk + dBC6Y0d3d00*u3d00dk + dBC6Y0d4d10*u4d10dk + dBC6Y0d4d00*u4d00dk + dBC6Y0d5d10*u5d10dk + dBC6Y0d5d00*u5d00dk + dBC6Y0d6d10*u6d10dk + dBC6Y0d6d00*u6d00dk + dBC6Y0d7d10*u7d10dk + dBC6Y0d7d00*u7d00dk + dBC6Y0d8d10*u8d10dk + dBC6Y0d8d00*u8d00dk + dBC6Y0d9d10*u9d10dk + dBC6Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC6Y0d0d01*u0d01dk + dBC6Y0d0d00*u0d00dk + dBC6Y0d1d01*u1d01dk + dBC6Y0d1d00*u1d00dk + dBC6Y0d2d01*u2d01dk + dBC6Y0d2d00*u2d00dk + dBC6Y0d3d01*u3d01dk + dBC6Y0d3d00*u3d00dk + dBC6Y0d4d01*u4d01dk + dBC6Y0d4d00*u4d00dk + dBC6Y0d5d01*u5d01dk + dBC6Y0d5d00*u5d00dk + dBC6Y0d6d01*u6d01dk + dBC6Y0d6d00*u6d00dk + dBC6Y0d7d01*u7d01dk + dBC6Y0d7d00*u7d00dk + dBC6Y0d8d01*u8d01dk + dBC6Y0d8d00*u8d00dk + dBC6Y0d9d01*u9d01dk + dBC6Y0d9d00*u9d00dk;
			
			//----- BC7 -----
			BC7Y0       =       BC7Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d0d00 = dBC7Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d0d10 = dBC7Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d0d01 = dBC7Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d1d00 = dBC7Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d1d10 = dBC7Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d1d01 = dBC7Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d2d00 = dBC7Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d2d10 = dBC7Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d2d01 = dBC7Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d3d00 = dBC7Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d3d10 = dBC7Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d3d01 = dBC7Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d4d00 = dBC7Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d4d10 = dBC7Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d4d01 = dBC7Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d5d00 = dBC7Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d5d10 = dBC7Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d5d01 = dBC7Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d6d00 = dBC7Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d6d10 = dBC7Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d6d01 = dBC7Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d7d00 = dBC7Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d7d10 = dBC7Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d7d01 = dBC7Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d8d00 = dBC7Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d8d10 = dBC7Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d8d01 = dBC7Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC7Y0d9d00 = dBC7Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d9d10 = dBC7Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC7Y0d9d01 = dBC7Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 7;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC7Y0d7d00*a7d0[a] + dBC7Y0d7d10*a7d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +7;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC7Y0d2d01*b2d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +2;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC7Y0);
			Dxvec[k +pn] = dBC7Y0d0d10*u0d10dk + dBC7Y0d0d00*u0d00dk + dBC7Y0d1d10*u1d10dk + dBC7Y0d1d00*u1d00dk + dBC7Y0d2d10*u2d10dk + dBC7Y0d2d00*u2d00dk + dBC7Y0d3d10*u3d10dk + dBC7Y0d3d00*u3d00dk + dBC7Y0d4d10*u4d10dk + dBC7Y0d4d00*u4d00dk + dBC7Y0d5d10*u5d10dk + dBC7Y0d5d00*u5d00dk + dBC7Y0d6d10*u6d10dk + dBC7Y0d6d00*u6d00dk + dBC7Y0d7d10*u7d10dk + dBC7Y0d7d00*u7d00dk + dBC7Y0d8d10*u8d10dk + dBC7Y0d8d00*u8d00dk + dBC7Y0d9d10*u9d10dk + dBC7Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC7Y0d0d01*u0d01dk + dBC7Y0d0d00*u0d00dk + dBC7Y0d1d01*u1d01dk + dBC7Y0d1d00*u1d00dk + dBC7Y0d2d01*u2d01dk + dBC7Y0d2d00*u2d00dk + dBC7Y0d3d01*u3d01dk + dBC7Y0d3d00*u3d00dk + dBC7Y0d4d01*u4d01dk + dBC7Y0d4d00*u4d00dk + dBC7Y0d5d01*u5d01dk + dBC7Y0d5d00*u5d00dk + dBC7Y0d6d01*u6d01dk + dBC7Y0d6d00*u6d00dk + dBC7Y0d7d01*u7d01dk + dBC7Y0d7d00*u7d00dk + dBC7Y0d8d01*u8d01dk + dBC7Y0d8d00*u8d00dk + dBC7Y0d9d01*u9d01dk + dBC7Y0d9d00*u9d00dk;
			
			//----- BC8 -----
			BC8Y0       =       BC8Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d0d00 = dBC8Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d0d10 = dBC8Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d0d01 = dBC8Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d1d00 = dBC8Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d1d10 = dBC8Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d1d01 = dBC8Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d2d00 = dBC8Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d2d10 = dBC8Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d2d01 = dBC8Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d3d00 = dBC8Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d3d10 = dBC8Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d3d01 = dBC8Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d4d00 = dBC8Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d4d10 = dBC8Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d4d01 = dBC8Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d5d00 = dBC8Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d5d10 = dBC8Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d5d01 = dBC8Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d6d00 = dBC8Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d6d10 = dBC8Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d6d01 = dBC8Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d7d00 = dBC8Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d7d10 = dBC8Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d7d01 = dBC8Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d8d00 = dBC8Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d8d10 = dBC8Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d8d01 = dBC8Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC8Y0d9d00 = dBC8Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d9d10 = dBC8Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC8Y0d9d01 = dBC8Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 8;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC8Y0d8d00*a8d0[a] + dBC8Y0d8d10*a8d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +8;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC8Y0d3d01*b3d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +3;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC8Y0);
			Dxvec[k +pn] = dBC8Y0d0d10*u0d10dk + dBC8Y0d0d00*u0d00dk + dBC8Y0d1d10*u1d10dk + dBC8Y0d1d00*u1d00dk + dBC8Y0d2d10*u2d10dk + dBC8Y0d2d00*u2d00dk + dBC8Y0d3d10*u3d10dk + dBC8Y0d3d00*u3d00dk + dBC8Y0d4d10*u4d10dk + dBC8Y0d4d00*u4d00dk + dBC8Y0d5d10*u5d10dk + dBC8Y0d5d00*u5d00dk + dBC8Y0d6d10*u6d10dk + dBC8Y0d6d00*u6d00dk + dBC8Y0d7d10*u7d10dk + dBC8Y0d7d00*u7d00dk + dBC8Y0d8d10*u8d10dk + dBC8Y0d8d00*u8d00dk + dBC8Y0d9d10*u9d10dk + dBC8Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC8Y0d0d01*u0d01dk + dBC8Y0d0d00*u0d00dk + dBC8Y0d1d01*u1d01dk + dBC8Y0d1d00*u1d00dk + dBC8Y0d2d01*u2d01dk + dBC8Y0d2d00*u2d00dk + dBC8Y0d3d01*u3d01dk + dBC8Y0d3d00*u3d00dk + dBC8Y0d4d01*u4d01dk + dBC8Y0d4d00*u4d00dk + dBC8Y0d5d01*u5d01dk + dBC8Y0d5d00*u5d00dk + dBC8Y0d6d01*u6d01dk + dBC8Y0d6d00*u6d00dk + dBC8Y0d7d01*u7d01dk + dBC8Y0d7d00*u7d00dk + dBC8Y0d8d01*u8d01dk + dBC8Y0d8d00*u8d00dk + dBC8Y0d9d01*u9d01dk + dBC8Y0d9d00*u9d00dk;
			
			//----- BC9 -----
			BC9Y0       =       BC9Y0out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d0d00 = dBC9Y0d0d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d0d10 = dBC9Y0d0d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d0d01 = dBC9Y0d0d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d1d00 = dBC9Y0d1d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d1d10 = dBC9Y0d1d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d1d01 = dBC9Y0d1d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d2d00 = dBC9Y0d2d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d2d10 = dBC9Y0d2d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d2d01 = dBC9Y0d2d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d3d00 = dBC9Y0d3d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d3d10 = dBC9Y0d3d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d3d01 = dBC9Y0d3d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d4d00 = dBC9Y0d4d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d4d10 = dBC9Y0d4d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d4d01 = dBC9Y0d4d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d5d00 = dBC9Y0d5d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d5d10 = dBC9Y0d5d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d5d01 = dBC9Y0d5d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d6d00 = dBC9Y0d6d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d6d10 = dBC9Y0d6d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d6d01 = dBC9Y0d6d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d7d00 = dBC9Y0d7d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d7d10 = dBC9Y0d7d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d7d01 = dBC9Y0d7d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d8d00 = dBC9Y0d8d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d8d10 = dBC9Y0d8d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d8d01 = dBC9Y0d8d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			dBC9Y0d9d00 = dBC9Y0d9d00out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d9d10 = dBC9Y0d9d10out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			dBC9Y0d9d01 = dBC9Y0d9d01out(xk, yk, r_H, M, alpha, beta, u0d00k, u0d10k, u0d01k, u1d00k, u1d10k, u1d01k, u2d00k, u2d10k, u2d01k, u3d00k, u3d10k, u3d01k, u4d00k, u4d10k, u4d01k, u5d00k, u5d10k, u5d01k, u6d00k, u6d10k, u6d01k, u7d00k, u7d10k, u7d01k, u8d00k, u8d10k, u8d01k, u9d00k, u9d10k, u9d01k);
			
			//-------------------------------------------------------
			pn = 9;
			
			for (a = 0; a <= r+onexside; ++a)
			{
				JacVal[nnztemp] = dBC9Y0d9d00*a9d0[a] + dBC9Y0d9d10*a9d1[a];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = (I*n +J+sx[a])*p +9;
				nnztemp += 1;
				
			}
			for (b = 0; b <= r+oneyside; ++b)
			{
				JacVal[nnztemp] = dBC9Y0d4d01*b4d1[b];
				JacRow[nnztemp] = (k +pn);
				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +4;
				nnztemp += 1;
			}
			
			bvec[k +pn]  = sgn*(BC9Y0);
			Dxvec[k +pn] = dBC9Y0d0d10*u0d10dk + dBC9Y0d0d00*u0d00dk + dBC9Y0d1d10*u1d10dk + dBC9Y0d1d00*u1d00dk + dBC9Y0d2d10*u2d10dk + dBC9Y0d2d00*u2d00dk + dBC9Y0d3d10*u3d10dk + dBC9Y0d3d00*u3d00dk + dBC9Y0d4d10*u4d10dk + dBC9Y0d4d00*u4d00dk + dBC9Y0d5d10*u5d10dk + dBC9Y0d5d00*u5d00dk + dBC9Y0d6d10*u6d10dk + dBC9Y0d6d00*u6d00dk + dBC9Y0d7d10*u7d10dk + dBC9Y0d7d00*u7d00dk + dBC9Y0d8d10*u8d10dk + dBC9Y0d8d00*u8d00dk + dBC9Y0d9d10*u9d10dk + dBC9Y0d9d00*u9d00dk;
			Dyvec[k +pn] = dBC9Y0d0d01*u0d01dk + dBC9Y0d0d00*u0d00dk + dBC9Y0d1d01*u1d01dk + dBC9Y0d1d00*u1d00dk + dBC9Y0d2d01*u2d01dk + dBC9Y0d2d00*u2d00dk + dBC9Y0d3d01*u3d01dk + dBC9Y0d3d00*u3d00dk + dBC9Y0d4d01*u4d01dk + dBC9Y0d4d00*u4d00dk + dBC9Y0d5d01*u5d01dk + dBC9Y0d5d00*u5d00dk + dBC9Y0d6d01*u6d01dk + dBC9Y0d6d00*u6d00dk + dBC9Y0d7d01*u7d01dk + dBC9Y0d7d00*u7d00dk + dBC9Y0d8d01*u8d01dk + dBC9Y0d8d00*u8d00dk + dBC9Y0d9d01*u9d01dk + dBC9Y0d9d00*u9d00dk;
			
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



int BVP_GENIC(struct param_type *params, struct ICparam_type *ICparams, double x[], double y[], double Psi[], double GRIDdx[], double GRIDdy[])
{
	int i,j,k,l, status;
	const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;
	const double r_H = (*ICparams).r_HICparam;
	const double M = (*ICparams).MICparam;
	const double alpha = (*ICparams).alphaICparam;
	const double beta = (*ICparams).betaICparam;
	const double delta = (*ICparams).deltaICparam;
	int pn;
	double xk, yk;
	
	//Generate IC
	for (i = 0; i < m; ++i)
	{
		y[i] = ((double) i/(m-1)) * (M_PI/2.0); //0 < y < M_PI/2
		yk = y[i];
		for (j = 0; j < n; ++j)
		{
			x[j] = ((double) j/(n-1));
			xk = x[j];
			k = (i*n+j)*p;
			
			pn = 0;
			Psi[k +pn] = IC0out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 1;
			Psi[k +pn] = IC1out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 2;
			Psi[k +pn] = IC2out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 3;
			Psi[k +pn] = IC3out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 4;
			Psi[k +pn] = IC4out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 5;
			Psi[k +pn] = IC5out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 6;
			Psi[k +pn] = IC6out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 7;
			Psi[k +pn] = IC7out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 8;
			Psi[k +pn] = IC8out(xk, yk, r_H, M, alpha, beta, delta);
			
			pn = 9;
			Psi[k +pn] = IC9out(xk, yk, r_H, M, alpha, beta, delta);
			
		}
	}
	
	//Calculate uniform grid x dimension
	for (j = 0; j < n; ++j)
	{
		if (j == 0)
		{
			GRIDdx[j] = (x[j+1] - x[j]);
		}
		else if (j == n-1)
		{
			GRIDdx[j] = (x[j] - x[j-1]);
		}
		else
		{
			GRIDdx[j] = (x[j+1] - x[j-1])/2.0;
		}
	}
	//Calculate uniform grid y dimension
	for (i = 0; i < m; ++i)
	{
		if (i == 0)
		{
			GRIDdy[i] = (y[i+1] - y[i]);
		}
		else if (i == m-1)
		{
			GRIDdy[i] = (y[i] - y[i-1]);
		}
		else
		{
			GRIDdy[i] = (y[i+1] - y[i-1])/2.0;
		}
	}
	
	
	return 0;
}




