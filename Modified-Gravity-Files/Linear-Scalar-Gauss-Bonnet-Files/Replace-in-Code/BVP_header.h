#ifndef HEADER_BVP
#define HEADER_BVP

/*
 *File List and Description
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_permute_vector.h>

//--------------- Initializing struct ---------------
struct param_type
{
	int nparam;
	int mparam;
	int pparam;
	int Nparam;
	int rparam;
	int nnzparam;
	double tolparam;
	double LStolparam;
	int LSimaxparam;
	double r_Hparam;
	double Mparam;
	double chiparam;
	double alphaparam;
	double betaparam;
	double wminparam;
	double wmaxparam;
	double wshrinkparam;
	double wgrowparam;
	double ICalphaparam;
	double ICchiparam;
};



//----- Main iteration loop -----
int BVP_main(const int plotvar, const int linsolvar, const int r, int n, int m, const int p, int N, const double tol, const double chi, const double ICdiff, const double r_H, const double M, const double alpha, const double beta, const double SCHWIC);

//----- Initial and boundary conditions -----
int BVP_GENIC(struct param_type *params, double x[], double y[], double Psi[], double GRIDdx[], double GRIDdy[]);
int BVP_GENBC(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[]);

//----- Output -----
int BVP_out(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], double Dx[], double Dy[], int it);
int BVP_sysout(struct param_type *params, double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], double Ddx[], double Dx[], double Ddy[], double Dy[], int it);

//----- Evaluation of system at current iteration -----
int BVP_GENsys(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[]);

//----- Grid resize -----
int BVP_grid(struct param_type *params, const int gridresizevar, double **x, double **y, double **Psi, double **GRIDdxvec, double **GRIDdyvec, double **Ddxvec, double **Ddyvec, double GRIDtol);
int BVP_realloc(struct param_type *params, double **b, double **Dx, double **Dy, double **dPsi, double **Ddx, double **Ddy, double **JacVal, int **JacRow, int **JacCol);

//----- Calculates new stepsize in 1 dimension -----
int BVP_resizecalc(struct param_type *params, const int gridresizevar, const double GRIDtol, double **x, double **y, double **dxold, double **dyold, double **Ddx, double **Ddy, double **dxnew, double **dynew);

//----- Calculates and interpolates new grid size -----
int BVP_gridsize(struct param_type *params, const int gridresizevar, double **x, double **y, double **Psi, double **dxold, double **dyold, double **dxnew, double **dynew);

//----- Interpolates system at new grid size -----
int BVP_interp(struct param_type *params, double x[], double y[], double Psi[], int nnew, int mnew, double xnew[], double ynew[], double Psinew[]);

//----- Evaluation of system at next iteration -----
int BVP_conv(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[], double bnorm, double *w, double dPsi[]);

//Sparse Matrix Storage Conversion to Numerical Recipes Code
int BVP_COOmodtoCOO(struct param_type *params, int *nnz, double **Val, int **arr1, int **arr2);
int COOinsertionSort(int *nnz, double **Val, int **arr1, int **arr2);
int BVP_COOtoCSR(struct param_type *params, int *nnz, double **Val, int **arr1, int **arr2, double **sa, unsigned long **ija);

//Sparse Matrix Calculation
void linbcg(double sa[], unsigned long ija[], unsigned long n, double b[], double x[], int itol, double tol, int itmax, int *iter, double *err);
void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n);
void dsprstx(double sa[], unsigned long ija[], double x[], double b[], unsigned long n);
void dsprsin(double a[], int n, double thresh, unsigned long nmax, double sa[],
	unsigned long ija[]);

//----- Calculates physical properties of solution -----
//int BVP_phys(void *params, double x[], double y[], double Psi[], double Jac[], double bvec[], double Dxvec[], double Dyvec[]);

//----- Calculates stencil -----
int steninit(int s[], int *oneside, int *twoside, const int k, const int r, const int N);

//----- Linear Solver -----
int BVP_LSsolver(struct param_type *params, int linsolvar, double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], double Ddx[], double Dx[], double Ddy[], double Dy[], double bnorm, double Dxnorm, double Dynorm, int it, int plotvar);

//----- LU DeCoMPosition -----
int NumRec_BCG(struct param_type *params, double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], int nnzJac);
int GSL_LUdcmp(const int N, double JacVal[], int JacRow[], int JacCol[], double x[], double b[], int nnzJac);
int GSL_GMRES(struct param_type *params, double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], int nnzJac);
int BVP_LUdcmp(const int N, double A[], double b[], double x[]);
int BVPGSL_LUdcmp(const int N, double JacVal[], int JacRow[], int JacCol[], double x[], double b[], const int nnzJac);

//----- Count non-zero elements of sparse matrix -----
int BVP_count(struct param_type *params);

//----- Misc functions -----
double norm(const int N, double v[]);

//--------------- FE Declaration ---------------
#include "../Funcs/Decl_CFuncs.c"
#include "../Funcs/Decl_JSFuncs.c"


void nrerror(char error_text[]);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);

#endif
