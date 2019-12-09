#include "BVP_header.h"

//----- List of defined functions -----
/*
BVP_LSsolver
GSL_GMRES
GSL_LUdcmp
NumRec_BCG
linbcg
asolve
atimes
snrm
dsprsin
dsprsax
dsprstx
nrerror
*dvector
free_dvector
*/

/*
linsolvar = :
0 = Combination of GSL library LU decomposition algorithm and backsubstitution algorithm from Numerical Recipes
1 = LU decomposition solely from GSL (Calls fortran wrapped backsubstitution which typecasts matrix size to int which can cause overflow for large grid size)
2 = QR decomposition using GPUs in CUDA. can only be uncommented on GPU hardware.
3 = GMRES from GSL
4 = Lin BCG from Numerical Recpies
5 = LU decomposition solely from Numerical Recipes
See Linear Solvers.pdf for mre details on some of the algorithms
*/


int BVP_LSsolver(struct param_type *params, int linsolvar, double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], double Ddx[], double Dx[], double Ddy[], double Dy[], double bnorm, double Dxnorm, double Dynorm, int it, int plotvar)
{
	int i,j,k,l,status;
	long krow;
	long kcol;
	const int N = (*params).Nparam;
	const long Nl = (long) N;
	int nnzJac = (*params).nnzparam;
	const double tol = (*params).tolparam;
	const double LStol = (*params).LStolparam;
	const int LSimax = (*params).LSimaxparam;

	printf("Begin Solver for nxm*p = N = %i. Lin solver variable = %i.\n", N, linsolvar);
	printf("nnz = %i. N^2 = %lu. rough sparsity = %.4f %%.\n",nnzJac, Nl*Nl, (float) (nnzJac)/(Nl*Nl)*100.0);

	//---------- LU decomposition combination of GSL and BVP ----------
	if (linsolvar == 0)
	{
		status = BVPGSL_LUdcmp(N, JacVal, JacRow, JacCol, dPsi, b, nnzJac);
		if (status != 0)
		{
			printf("ERROR: Chosen linear solver.\n");
			return -1;
		}

		if (Dxnorm < tol)
		{
			printf("|Dx| < tol.\n");
		}
		else
		{
			status = BVPGSL_LUdcmp(N, JacVal, JacRow, JacCol, Ddx, Dx, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dynorm < tol)
		{
			printf("|Dy| < tol.\n");
		}
		else
		{
			status = BVPGSL_LUdcmp(N, JacVal, JacRow, JacCol, Ddy, Dy, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}


		//Cleanup


	}


	//---------- LU decomposition from GSL ----------
	else if (linsolvar == 1)
	{
		status = GSL_LUdcmp(N, JacVal, JacRow, JacCol, dPsi, b, nnzJac);
		if (status != 0)
		{
			printf("ERROR: Chosen linear solver.\n");
			return -1;
		}
		if (Dxnorm < tol)
		{
			printf("|Dx| < tol.\n");
		}
		else
		{
			status = GSL_LUdcmp(N, JacVal, JacRow, JacCol, Ddx, Dx, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}
		if (Dynorm < tol)
		{
			printf("|Dy| < tol.\n");
		}
		else
		{
			status = GSL_LUdcmp(N, JacVal, JacRow, JacCol, Ddy, Dy, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}
	}


	//---------- QR decomposition from CUDA ----------
	else if (linsolvar == 2)
	{
		double *Jac_dense;
	  Jac_dense = calloc(Nl*Nl, sizeof(double));
		if(!Jac_dense){
	        printf("Error! Memory not allocated. Exiting");
			return -1;
	    }
	  for (i = 0; i < nnzJac; ++i)
	  {
	    kcol = (long) JacCol[i] * Nl + (long) JacRow[i];
	    Jac_dense[kcol] = JacVal[i];
	    //printf("Jac[%2lu] = %.2f.\n",kcol,Jac_dense[kcol]);
	  }

		status = 0;
		//status = QRsolve(N, Jac_dense, b, dPsi);
		if (status != 0)
		{
			printf("ERROR: Chosen linear solver.\n");
			return -1;
		}

		if (Dxnorm < tol)
		{
			printf("|Dx| < tol.\n");
		}
		else
		{
			//status = QRsolve(N, Jac_CUdense, Dx, Ddx);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dynorm < tol)
		{
			printf("|Dy| < tol.\n");
		}
		else
		{
			//status = QRsolve(N, Jac_CUdense, Dy, Ddy);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}


		//Cleanup
		free(Jac_dense);

	}


	//---------- GMRES from GSL ----------
	else if (linsolvar == 3)
	{
		if (bnorm*LStol < DBL_EPSILON)
		{
			printf("Linear Solver: ||b||*tol = %e < DBL_EPSILON. Returning zero vector.\n",bnorm*LStol);
			for (i = 0; i < N; ++i)
			{
				dPsi[i] = 0.0;
			}
		}
		else
		{
			status = GSL_GMRES(params, JacVal, JacRow, JacCol, dPsi, b, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dxnorm*LStol < DBL_EPSILON)
		{
			printf("Linear Solver: ||Dx||*tol = %e < DBL_EPSILON. Returning zero vector.\n",Dxnorm*LStol);
			for (i = 0; i < N; ++i)
			{
				Ddx[i] = 0.0;
			}
		}
		else
		{
			status = GSL_GMRES(params, JacVal, JacRow, JacCol, Ddx, Dx, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dynorm*LStol < DBL_EPSILON)
		{
			printf("Linear Solver: ||Dy||*tol = %e < DBL_EPSILON. Returning zero vector.\n",Dynorm*LStol);
			for (i = 0; i < N; ++i)
			{
				Ddy[i] = 0.0;
			}
		}
		else
		{
			status = GSL_GMRES(params, JacVal, JacRow, JacCol, Ddy, Dy, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}
	}


	//---------- Lin BCG from Numerical Recipes ----------
	else if (linsolvar == 4)
	{
		if (bnorm*LStol < DBL_EPSILON)
		{
			printf("Linear Solver: ||b||*tol = %e < DBL_EPSILON. Returning zero vector.\n",bnorm*LStol);
			for (i = 0; i < N; ++i)
			{
				dPsi[i] = 0.0;
			}
		}
		else
		{
			status = NumRec_BCG(params, JacVal, JacRow, JacCol, dPsi, b, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dxnorm*LStol < DBL_EPSILON)
		{
			printf("Linear Solver: ||Dx||*tol = %e < DBL_EPSILON. Returning zero vector.\n",Dxnorm*LStol);
			for (i = 0; i < N; ++i)
			{
				Ddx[i] = 0.0;
			}
		}
		else
		{
			status = NumRec_BCG(params, JacVal, JacRow, JacCol, Ddx, Dx, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dynorm*LStol < DBL_EPSILON)
		{
			printf("Linear Solver: ||Dy||*tol = %e < DBL_EPSILON. Returning zero vector.\n",Dynorm*LStol);
			for (i = 0; i < N; ++i)
			{
				Ddy[i] = 0.0;
			}
		}
		else
		{
			status = NumRec_BCG(params, JacVal, JacRow, JacCol, Ddy, Dy, nnzJac);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}
	}


	//---------- LU decomposition from Numerical Recipes ----------
	else if (linsolvar == 5)
	{
		double *Jac_dense;
	  Jac_dense = calloc(Nl*Nl, sizeof(double));
		if(!Jac_dense){
	        printf("Error! Memory not allocated. Exiting");
			return -1;
	    }
	  for (i = 0; i < nnzJac; ++i)
	  {
	    krow = (long) JacRow[i] * Nl + (long) JacCol[i];
	    Jac_dense[krow] = JacVal[i];
	    //printf("Jac[%2lu] = %.2f.\n",krow,Jac_dense[krow]);
	  }


		status = BVP_LUdcmp(N, Jac_dense, b, dPsi);
		if (status != 0)
		{
			printf("ERROR: Chosen linear solver.\n");
			return -1;
		}

		if (Dxnorm < tol)
		{
			printf("|Dx| < tol.\n");
		}
		else
		{
			status = BVP_LUdcmp(N, Jac_dense, Dx, Ddx);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}

		if (Dynorm < tol)
		{
			printf("|Dy| < tol.\n");
		}
		else
		{
			status = BVP_LUdcmp(N, Jac_dense, Dy, Ddy);
			if (status != 0)
			{
				printf("ERROR: Chosen linear solver.\n");
				return -1;
			}
		}


		//Cleanup
		free(Jac_dense);

	}


	if (plotvar == 1)
	{
		status = BVP_sysout(params, JacVal, JacRow, JacCol, dPsi, b, Ddx, Dx, Ddy, Dy, it); //Output
	}


	//Cleanup


	return 0;
}





//------------------------------------------
//           Function Definitions
//------------------------------------------
int BVPGSL_LUdcmp(const int N, double JacVal[], int JacRow[], int JacCol[], double x[], double b[], const int nnzJac)
{
	int status;
	long il,jl;
	const long Nl = (long) N;
	const long nnzJacl = (long) nnzJac;

	clock_t begin, end;
	double time_spent1, time_spent2;
	int sigperm;

	gsl_permutation *p = gsl_permutation_alloc(Nl);
	gsl_matrix *JacDense = gsl_matrix_calloc(Nl, Nl);
	gsl_spmatrix *JacCOO = gsl_spmatrix_alloc_nzmax(Nl, Nl, nnzJacl, GSL_SPMATRIX_TRIPLET);
	gsl_vector *bGSL = gsl_vector_alloc(Nl);
	gsl_vector *xGSL = gsl_vector_alloc(Nl);

	//Import sparse matrix in COO format
	for (il = 0; il < nnzJacl; ++il)
	{
		gsl_spmatrix_set(JacCOO, JacRow[il], JacCol[il], JacVal[il]);
	}

	//Construct r.h.s. vector
	for (il = 0; il < Nl; ++il)
	{
		gsl_vector_set(bGSL, il, b[il]);
	}

	//set initial guess to zero
	gsl_vector_set_zero(xGSL);

	//Convert COO to dense format
	status = gsl_spmatrix_sp2d(JacDense, JacCOO);
	if (status != 0)
	{
		return -1;
	}


	//Solve linear system using GSL
	begin = clock();
	status = gsl_linalg_LU_decomp(JacDense, p, &sigperm);
	if (status != 0)
	{
		return -2;
	}
	end = clock();
	time_spent1 = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for GSL LU_decomp.\n",time_spent1/60.0,time_spent1);


	//Solve by backsubstitution using Numerical Recipes
	status = gsl_vector_memcpy (xGSL, bGSL);
	if (status != 0)
	{
		return -1;
	}
	status = gsl_permute_vector (p, xGSL);
	if (status != 0)
	{
		return -1;
	}

	for (il = 0; il < Nl; il++)
	{
		//x[il] = b[il];
  	//y[i] = b[i]/al[i*N+i];
		for (jl = 0; jl < il; jl++)
		{
			gsl_vector_set(xGSL,il, gsl_vector_get(xGSL,il) - gsl_matrix_get(JacDense,il,jl)*gsl_vector_get(xGSL,jl));
      //x[il] = x[il] - (JacDense->data)[il*Nl+jl]*x[jl];
      //y[i] = y[i] - al[i*N+j]*y[j]/al[i*N+i];
		}
		//printf("y[%3lu] = %7.2f\n",il,x[il]);
	}
	for (il = Nl-1l; il > -1; il--)
	{
		gsl_vector_set(xGSL,il, gsl_vector_get(xGSL,il)/gsl_matrix_get(JacDense,il,il));
		//x[il] = x[il]/(JacDense->data)[il*Nl+il];
    //x[i] = y[i]/be[i*N+i];
		for (jl = il+1; jl < Nl; jl++)
		{
			gsl_vector_set(xGSL,il, gsl_vector_get(xGSL,il) - gsl_matrix_get(JacDense,il,jl) * gsl_vector_get(xGSL,jl) / gsl_matrix_get(JacDense,il,il));
			//x[il] = x[il] - (JacDense->data)[il*Nl+jl]*x[jl]/(JacDense->data)[il*Nl+il];
      //x[i] = x[i] - be[i*N+j]*x[j]/be[i*N+i];
		}
		//printf("x[%3lu] = %7.2f\n",il,x[il]);
	}
	end = clock();
	time_spent2 = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for BVP LU_solve.\n",(time_spent2-time_spent1)/60.0,(time_spent2-time_spent1));
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for BVP/GSL LU decomposition Solver.\n",time_spent2/60.0,time_spent2);


	//Output solution
	for (il = 0; il < Nl; ++il)
	{
		x[il] = gsl_vector_get(xGSL,il);
		//printf("%11.4e\n", gsl_vector_get(xGSL, il));
	}

	//Cleanup
	gsl_permutation_free(p);
	gsl_matrix_free(JacDense);
	gsl_spmatrix_free(JacCOO);
	gsl_vector_free(bGSL);
	gsl_vector_free(xGSL);


	return 0;
}

int BVP_LUdcmp(const int N, double A[], double b[], double x[])
{
	int status;
  long il,jl,kl;
	const long Nl = (long) N;

	clock_t begin, end;
	double time_spent1, time_spent2;


	begin = clock();
	for (jl = 0; jl < Nl; jl++)
	{
		for (il = 0; il < jl+1; il++)
		{
			//be[i*N+j] = A[i*N+j];
			for (kl = 0; kl < il; kl++)
			{
			  A[il*Nl+jl] = A[il*Nl+jl] - A[il*Nl+kl]*A[kl*Nl+jl];
        //be[i*N+j] = be[i*N+j] - al[i*N+k]*be[k*N+j];
			}
			//printf("be[%3lu,%3lu] = %7.2f\n",il,jl,A[il*Nl+jl]);
		}
		for (il = jl+1; il < Nl; il++)
		{
			A[il*Nl+jl] = A[il*Nl+jl]/A[jl*Nl+jl];
      //al[i*N+j] = A[i*N+j]/be[j*N+j];
			for (kl = 0; kl < jl; kl++)
			{
        A[il*Nl+jl] = A[il*Nl+jl] - A[il*Nl+kl]*A[kl*Nl+jl]/A[jl*Nl+jl];
        //al[i*N+j] = al[i*N+j] - al[i*N+k]*be[k*N+j]/be[j*N+j];
			}
			//printf("al[%3lu,%3lu] = %7.2f\n",il,jl,A[il*Nl+jl]);
		}
	}
	end = clock();
	time_spent1 = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for BVP LU decomp.\n",time_spent1/60.0,time_spent1);

	for (il = 0; il < Nl; ++il)
	{
		x[il] = b[il];
	}

	for (il = 0; il < Nl; il++)
	{
		x[il] = b[il];
    //y[i] = b[i]/al[i*N+i];
		for (jl = 0; jl < il; jl++)
		{
			x[il] = x[il] - A[il*Nl+jl]*x[jl];
      //y[i] = y[i] - al[i*N+j]*y[j]/al[i*N+i];
		}
		//printf("y[%3lu] = %7.2f\n",il,x[il]);
	}
	for (il = Nl-1l; il > -1; il--)
	{
		x[il] = x[il]/A[il*Nl+il];
    //x[i] = y[i]/be[i*N+i];
		for (jl = il+1; jl < Nl; jl++)
		{
			x[il] = x[il] - A[il*Nl+jl]*x[jl]/A[il*Nl+il];
      //x[i] = x[i] - be[i*N+j]*x[j]/be[i*N+i];
		}
		//printf("x[%3lu] = %7.2f\n",il,x[il]);
	}
	end = clock();
	time_spent2 = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for BVP LU decomposition Solver.\n",(time_spent2-time_spent1)/60.0,(time_spent2-time_spent1));
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for BVP LU decomposition Solver.\n",time_spent2/60.0,time_spent2);

  //Print L U
	/*
  for (il = 0; il < Nl; ++il)
  {
    for (jl = 0; jl < Nl; ++jl)
    {
      kl = il*Nl+jl;
      printf("A[%3lu,%3lu] = %7.2f.\n",il,jl,A[kl]);
    }
  }
	*/


	//Cleanup


	return 0;
}

int GSL_LUdcmp(const int N, double JacVal[], int JacRow[], int JacCol[], double x[], double b[], int nnzJac)
{
	int i,j,k,l, status;
	const long int Nl = (long int) N;

	clock_t begin, end;
	double time_spent1, time_spent2;
	int sigperm;

	gsl_permutation *p = gsl_permutation_alloc(Nl);
	gsl_matrix *JacDense = gsl_matrix_calloc(Nl, Nl);
	gsl_spmatrix *JacCOO = gsl_spmatrix_alloc_nzmax(Nl, Nl, nnzJac, GSL_SPMATRIX_TRIPLET);
	gsl_vector *bGSL = gsl_vector_alloc(Nl);
	gsl_vector *xGSL = gsl_vector_alloc(Nl);

	//Import sparse matrix in COO format
	for (i = 0; i < nnzJac; ++i)
	{
		gsl_spmatrix_set(JacCOO, JacRow[i], JacCol[i], JacVal[i]);
	}

	//Construct r.h.s. vector
	for (i = 0; i < (int) Nl; ++i)
	{
		gsl_vector_set(bGSL, i, b[i]);
	}

	//set initial guess to zero
	gsl_vector_set_zero(xGSL);

	//Convert COO to dense format
	status = gsl_spmatrix_sp2d(JacDense, JacCOO);
	if (status != 0)
	{
		return -1;
	}


	//Solve linear system
	begin = clock();
	status = gsl_linalg_LU_decomp(JacDense, p, &sigperm);
	if (status != 0)
	{
		return -2;
	}
	end = clock();
	time_spent1 = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for GSL LU_decomp.\n",time_spent1/60.0,time_spent1);

	status = gsl_linalg_LU_solve(JacDense, p, bGSL, xGSL);
	if (status != 0)
	{
		return -2;
	}
	end = clock();
	time_spent2 = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for GSL LU_solve.\n",(time_spent2-time_spent1)/60.0,(time_spent2-time_spent1));
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for GSL LU decompositio Solver.\n",time_spent2/60.0,time_spent2);


	//Output solution
	for (i = 0; i < (int) Nl; ++i)
	{
		x[i] = gsl_vector_get(xGSL,i);
		//printf("%11.4e\n", gsl_vector_get(xGSL, i));
	}

	//Cleanup
	gsl_permutation_free(p);
	gsl_matrix_free(JacDense);
	gsl_spmatrix_free(JacCOO);
	gsl_vector_free(bGSL);
	gsl_vector_free(xGSL);


	return 0;
}

int GSL_GMRES(struct param_type *params, double JacVal[], int JacRow[], int JacCol[], double x[], double b[], int nnzJac)
{
	int i,j,k,l, status;
	const int N = (*params).Nparam;
	const double LStol = (*params).LStolparam;
	const int LSimax = (*params).LSimaxparam;
	const int nKrylov = N;

	clock_t begin, end;
	double time_spent;
	double residual;

	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, N, nKrylov);
	size_t iter = 0;

	gsl_spmatrix *JacCOO = gsl_spmatrix_alloc_nzmax(N, N, nnzJac, GSL_SPMATRIX_TRIPLET);
	gsl_spmatrix *JacCRS;
	gsl_vector *bGSL = gsl_vector_alloc(N);
	gsl_vector *xGSL = gsl_vector_alloc(N);

	//Import sparse matrix in COO format
	for (i = 0; i < nnzJac; ++i)
	{
		gsl_spmatrix_set(JacCOO, JacRow[i], JacCol[i], JacVal[i]);
	}

	//Construct r.h.s. vector
	for (i = 0; i < N; ++i)
	{
		gsl_vector_set(bGSL, i, b[i]);
	}

	//set initial guess to zero
	gsl_vector_set_zero(xGSL);

	//Convert COO to CRS format
	JacCRS = gsl_spmatrix_crs(JacCOO);

	//Solve linear system
	begin = clock();
	do
	{
		status = gsl_splinalg_itersolve_iterate(JacCRS, bGSL, LStol, xGSL, work);
		//print out residual norm ||A*x - b||
		residual = gsl_splinalg_itersolve_normr(work);
		//fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);
		printf("iter %zu residual = %.12e\n", iter, residual);

		if (status == GSL_SUCCESS)
		{
			//fprintf(stderr, "Converged\n");
			printf("Converged\n");
		}
	}
	while (status == GSL_CONTINUE && ++iter < LSimax);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time: %10.3f [m], %11.4e [s] for GSL GMRES Solver.\n",time_spent/60.0,time_spent);
	if (iter == LSimax)
	{
		printf("\n\nERROR!! Number of iterations is equal to max iterations, %i.\n\n\n",LSimax);
		return -1;
	}

	//output solution
	for (i = 0; i < N; ++i)
	{
		x[i] = gsl_vector_get(xGSL,i);
		//printf("%11.4e\n", gsl_vector_get(xGSL, i));
	}


	//Cleanup
	gsl_splinalg_itersolve_free(work);
	gsl_spmatrix_free(JacCOO);
	gsl_spmatrix_free(JacCRS);
	gsl_vector_free(bGSL);
	gsl_vector_free(xGSL);


	return 0;
}

int NumRec_BCG(struct param_type *params, double JacVal[], int JacRow[], int JacCol[], double x[], double b[], int nnzJac)
{
	int i,j,k,l, status;
	const int N = (*params).Nparam;
	const double LStol = (*params).LStolparam;
	const int LSimax = (*params).LSimaxparam;

	clock_t begin, end;
	double time_spent;
	int nnzCSR;
	int iter, itertot;
	double err = 1.0;

	int ITOL = 1;
	/*
	Tolerance condition:
	1 = |A x - b| / |b| < TOL
	2 = |A^{-1} (A x - b)| / |A^{-1} b| < TOL
	3 = "Uses its own estimate of the error in x, and requires its magnitude, divided by the magnitude of x, to be less than TOL
	4 = same as 3 except that the largest abs value component of the error and largest component of x used instead of vector magnitude. L_\infty norm vs L_2 norm.
	*/

	nnzCSR = nnzJac + 2;
	//---------- Calloc BEGIN ----------
	unsigned long *ija;
	ija = (unsigned long*)calloc((size_t) nnzCSR , sizeof(unsigned long));
	double *sa;
	sa = (double*)calloc((size_t) nnzCSR , sizeof(double));
	double *bSPRS;
	bSPRS = (double*)calloc((size_t) N+1 , sizeof(double));
	double *xSPRS;
	xSPRS = (double*)calloc((size_t) N+1 , sizeof(double));
	//---------- Calloc END ----------

	//status = BVP_COOtoCSR(params, &nnzJac, &JacVal, &JacRow, &JacCol, &sa, &ija);
	gsl_matrix *JacDense = gsl_matrix_alloc(N, N);
	gsl_spmatrix *JacCOO = gsl_spmatrix_alloc_nzmax(N, N, nnzJac, GSL_SPMATRIX_TRIPLET);
	double *JacFull;
	JacFull = (double*)calloc((size_t) N*N , sizeof(double));

	//Import sparse matrix in COO format
	for (i = 0; i < nnzJac; ++i)
	{
		gsl_spmatrix_set(JacCOO, JacRow[i], JacCol[i], JacVal[i]);
	}

	//Convert COO to dense format
	status = gsl_spmatrix_sp2d(JacDense, JacCOO);

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			k = i*N+j;
			JacFull[k] = gsl_matrix_get(JacDense,i,j);
		}
	}
	dsprsin(JacFull, N, DBL_EPSILON, nnzCSR, sa, ija);

	/*
	printf("CSR format.\n");
	printf("index	ija	sa\n");
	for (i = 0; i < N; ++i)//nnzCSR; ++i)
	{
		printf("%2i\t%lu\t%11.4e.\n",i,ija[i],sa[i]);
	}
	*/

	for (i = 0; i < N; ++i)
	{
		bSPRS[i+1] = b[i];
		//printf("b[%i] = %6.3f\n",i+1,bSPRS[i+1]);
	}

	//-------------------- Solve --------------------
	itertot = 0;
	begin = clock();
	while ((err > LStol) && (itertot < LSimax))
	{
		linbcg(sa,ija,N,bSPRS,xSPRS,ITOL,LStol,2*LSimax,&iter,&err);
		printf("%s %6d, %s %15e\n","Iteration",iter+itertot*LSimax,"Estimated error:",err);
		itertot++;
	}
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Elapsed Time for Numerical Recipes BCG Solver = %10.3f [m], %11.4e [s].\n",time_spent/60.0,time_spent);
	if (itertot == LSimax)
	{
		printf("\n\nERROR!! Number of BCG restarts is equal to max iterations, %i.\n\n\n",LSimax);
	}
	printf("%s %15e\n","Estimated error:",err);
	printf("%s %6d %s %6d\n","Number of restarts:",itertot,"of",LSimax);
	printf("%s %6d\n","Total iterations needed:",iter + itertot*LSimax);
	//printf("\nSolution vector:\n");
	//for (i = 1; i < N+1; ++i)
	//{
		//printf("x[%2i] = %7.4f\n",i,xSPRS[i]);
	//}
	//printf("\n");
	for (i = 1; i < N+1; ++i)
	{
		x[i-1] = xSPRS[i];
	}


	//Check
	/*
	double *bcmp;
	bcmp = (double*)calloc((size_t) N+1 , sizeof(double));
	dsprsax(sa,ija,xSPRS,bcmp,N);
	printf("test of solution vector:\n");
	printf("%9s %12s %5s\n","a*x","b","diff");
	for (i = 1; i < N+1; ++i)
	{
		printf("%11.4e\t%11.4e\t%11.4e\n",bcmp[i],bSPRS[i],fabs(bcmp[i]-bSPRS[i]));
	}
	free(bcmp);
	*/

	//Cleanup
	gsl_matrix_free(JacDense);
	gsl_spmatrix_free(JacCOO);
	free(JacFull);
	free(ija);
	free(sa);
	free(bSPRS);
	free(xSPRS);


	return 0;
}

#define EPS 1.0e-14
#define NR_END 1
#define FREE_ARG char*

void linbcg(double sa[], unsigned long ija[], unsigned long n, double b[], double x[], int itol, double tol,
	int itmax, int *iter, double *err)
{
	void asolve(double sa[], unsigned long ija[], unsigned long n, double b[], double x[], int itrnsp);
	void atimes(double sa[], unsigned long ija[], unsigned long n, double x[], double r[], int itrnsp);
	double snrm(unsigned long n, double sx[], int itol);
	unsigned long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p=dvector(1,n);
	pp=dvector(1,n);
	r=dvector(1,n);
	rr=dvector(1,n);
	z=dvector(1,n);
	zz=dvector(1,n);

	*iter=0;
	atimes(sa,ija,n,x,r,0);
	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	atimes(sa,ija,n,r,rr,0); //Uncommented means "minimum residual method"
	if (itol == 1) {
		bnrm=snrm(n,b,itol);
		asolve(sa,ija,n,r,z,0);
	}
	else if (itol == 2) {
		asolve(sa,ija,n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(sa,ija,n,r,z,0);
	}
	else if (itol == 3 || itol == 4) {
		asolve(sa,ija,n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(sa,ija,n,r,z,0);
		znrm=snrm(n,z,itol);
	} else nrerror((char*)"illegal itol in linbcg");
	while (*iter <= itmax) {
		++(*iter);
		asolve(sa,ija,n,rr,zz,1);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(sa,ija,n,p,z,0);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(sa,ija,n,pp,zz,1);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(sa,ija,n,r,z,0);
		if (itol == 1)
			*err=snrm(n,r,itol)/bnrm;
 		else if (itol == 2)
			*err=snrm(n,z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
		//printf("iter=%4d err=%12.6f\n",*iter,*err);
	if (*err <= tol) break;
	}

	free_dvector(p,1,n);
	free_dvector(pp,1,n);
	free_dvector(r,1,n);
	free_dvector(rr,1,n);
	free_dvector(z,1,n);
	free_dvector(zz,1,n);
}

void asolve(double sa[], unsigned long ija[], unsigned long n, double b[], double x[], int itrnsp)
{
	unsigned long i;

	for(i=1;i<=n;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]); //Preconditioner A-tilde is diagonal elements of b
	//for(i=1;i<=n;i++) x[i]=b[i]; //Preconditioner A-tilde is identity matrix. Equivalent to no preconditioner.
}

void atimes(double sa[], unsigned long ija[], unsigned long n, double x[], double r[], int itrnsp)
{
	if (itrnsp) dsprstx(sa,ija,x,r,n);
	else dsprsax(sa,ija,x,r,n);
}

double snrm(unsigned long n, double sx[], int itol)
{
	unsigned long i,isamax;
	double ans;

	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}

void dsprsin(double a[], int n, double thresh, unsigned long nmax, double sa[],
	unsigned long ija[])
{
	void nrerror(char error_text[]);
	int i,j;
	unsigned long k;

	for (j=1;j<=n;j++) sa[j]=a[(j-1)*n+(j-1)];
	ija[1]=n+2;
	k=n+1;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) {
			if (fabs(a[(i-1)*n+(j-1)]) >= thresh && i != j) {
				if (++k > nmax) nrerror((char*)"sprsin: nmax too small");
				sa[k]=a[(i-1)*n+(j-1)];
				ija[k]=j;
	//			printf("k = %ld\n",k);
			}
		}
		ija[i+1]=k+1;
	}
}

void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,k;

	if (ija[1] != n+2) nrerror((char*)"dsprsax: mismatched vector and matrix");
	for (i=1;i<=n;i++)
	{
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++) b[i] += sa[k]*x[ija[k]];
	}
}

void dsprstx(double sa[], unsigned long ija[], double x[], double b[], unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,j,k;
	if (ija[1] != n+2) nrerror((char*)"mismatched vector and matrix in dsprstx");
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i];
	for (i=1;i<=n;i++) {
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror((char*)"allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

#undef EPS
#undef NR_END
#undef FREE_ARG
