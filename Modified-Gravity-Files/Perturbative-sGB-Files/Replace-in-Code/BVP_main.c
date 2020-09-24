/*
Andrew Sullivan
Montana State University
Newton-Rapson method for solving PDE BVP in 2D as coupled equations
Procedure:
Notes:
*/

#include "BVP_header.h"

//----- List of defined functions -----
/*
BVP_main
norm
BVP_count
BVP_COOmodtoCOO
COOinsertionSort
BVP_COOtoCSR
BVP_conv
*/

int BVP_main(const int plotvar, const int linsolvar, const int r, int n, int m, const int p, int N, const double tol, const double chi, const double ICalpha, const double r_H, const double M, const double alpha, const double beta, const double ICchi)
{
	//-------------------- Initialization --------------------
	int i,j,k,l, status; //For loop counters
	int it, gridit, gridresizevar; //Iteration number
	double bnorm, Jnorm, Dxnorm, Dynorm, unorm, dPsinorm, duDxnorm, duDynorm; //Sys norms
	double temp; //Temp variable
	int nnzJac, nnzJacCount, nnzJacTemp;
	clock_t begin, end;
	double time_spent;

	double w = 1.0; //Relaxation start value
	const double wmin = 5.0E-3; //Minimum w
	const double wmax = 1.0; //Max w
	const double wshrink = 0.5; //If not converging shrink w
	const double wgrow = 1.5; //If converging grow w

	const double LStol = 1.0E-7; //Linear solver relative tolerance
	const int bvpimax = 50; //max it
	const int griditmax = 5; //maximum number of grid resize attempts
	int LSimax = N; //Sparse linear solver iteration max

	//---------- Setting struct ----------
	struct param_type my_params; //parameter structure
	my_params.nparam = n;
	my_params.mparam = m;
	my_params.pparam = p;
	my_params.Nparam = N;
	my_params.rparam = r;
	my_params.tolparam = tol;
	my_params.LStolparam = LStol;
	my_params.LSimaxparam = LSimax;
	my_params.chiparam = chi;
	my_params.r_Hparam = r_H;
	my_params.Mparam = M;
	my_params.alphaparam = alpha;
	my_params.betaparam = beta;
	my_params.wminparam = wmin;
	my_params.wmaxparam = wmax;
	my_params.wshrinkparam = wshrink;
	my_params.wgrowparam = wgrow;
	my_params.ICalphaparam = ICalpha;
	my_params.ICchiparam = ICchi;

	nnzJacCount = BVP_count(&my_params); //calculate estimated Jacobian nonzeros
	my_params.nnzparam = nnzJacCount;
	nnzJac = nnzJacCount;
	printf("Estimated modified non-zero count, and initial nnzparam: %i.\n",nnzJacCount);

	//----- Calloc BEGIN -----
	status = 0;
	//Grid
	double *x;
	x = (double*)calloc(n, sizeof(double));
	if (!x)
		status = -1;
	double *y;
	y = (double*)calloc(m, sizeof(double));
	if (!y)
		status = -1;
	double *Psi;
	Psi = (double*)calloc(N, sizeof(double));
	if (!Psi)
		status = -1;
	double *GRIDdx;
	GRIDdx = (double*)calloc(n, sizeof(double));
	if (!GRIDdx)
		status = -1;
	double *GRIDdy;
	GRIDdy = (double*)calloc(m, sizeof(double));
	if (!GRIDdy)
		status = -1;
	//r.h.s. Newton's equation
	double *b;
	b = (double*)calloc(N, sizeof(double));
	if (!b)
		status = -1;
	double *Dx;
	Dx = (double*)calloc(N, sizeof(double));
	if (!Dx)
		status = -1;
	double *Dy;
	Dy = (double*)calloc(N, sizeof(double));
	if (!Dy)
		status = -1;
	//l.h.s of Newton's equation
	double *dPsi;
	dPsi = (double*)calloc(N, sizeof(double));
	if (!dPsi)
		status = -1;
	double *Ddx;
	Ddx = (double*)calloc(N, sizeof(double));
	if (!Ddx)
		status = -1;
	double *Ddy;
	Ddy = (double*)calloc(N, sizeof(double));
	if (!Ddy)
		status = -1;
	//Jacobian
	double *JacVal;
	JacVal = (double*)calloc(nnzJacCount, sizeof(double));
	if (!JacVal)
		status = -1;
	int *JacRow;
	JacRow = (int*)calloc(nnzJacCount, sizeof(int));
	if (!JacRow)
		status = -1;
	int *JacCol;
	JacCol = (int*)calloc(nnzJacCount, sizeof(int));
	if (!JacCol)
		status = -1;
	//----- Calloc END -----
	if (status != 0)
	{
		printf("Calloc failure in BVP_main!!\n");
		return -2;
	}


	//Initial Conditions
	status = BVP_GENIC(&my_params, x, y, Psi, GRIDdx, GRIDdy); //Calculate initial condition
	//status = BVP_out(&my_params, x, y, Psi, JacVal, JacRow, JacCol, dPsi, b, Dx, Dy, -3); //Converged solution in different directory


	for (it = 0; it <= bvpimax; ++it)
	{
		gridit = 0;
		printf("\n----- Iteration: %03i -----\n\n",it);

		do //---------- Grid Size Loop ----------
		{
			//Linear system evaluation
			begin = clock();
			nnzJacTemp = BVP_GENsys(&my_params, x, y, Psi, JacVal, JacRow, JacCol, b, Dx, Dy);
			end = clock();
			time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			printf("Elapsed Time: %10.3f [m], %11.4e [s] for Jacobian evaluation.\n",time_spent/60.0,time_spent);
			if (nnzJacTemp != nnzJacCount)
			{
				printf("ERROR!! the number of nonzeros from BVP_count should equal those from BVP_GENsys.\n");
				break;
			}
			//Convert from modified COO storage to COO storage
			begin = clock();
			status = BVP_COOmodtoCOO(&my_params, &nnzJacTemp, &JacVal, &JacRow, &JacCol);
			end = clock();
			time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			printf("Elapsed Time: %10.3f [m], %11.4e [s] for COO sort.\n",time_spent/60.0,time_spent);
			printf("Sorted COO (row major) format number of non-zeros: %i.\n",nnzJacTemp);
			if (nnzJac != nnzJacTemp)
			{
				nnzJac = nnzJacTemp;
				my_params.nnzparam = nnzJac;
				printf("Resizing number of nonzeros, nnzJac = %i.\n",nnzJac);
			}

			bnorm = norm(N, b);
			Jnorm = norm(nnzJac, JacVal);
			unorm = norm(N,Psi);
			Dxnorm = norm(N, Dx);
			Dynorm = norm(N, Dy);
			printf("||b||    = %11.4e.\n",bnorm);
			printf("||u||    = %11.4e.\n",unorm);
			printf("||J||    = %11.4e.\n",Jnorm);
			printf("||Dx||   = %11.4e.\n",Dxnorm);
			printf("||Dy||   = %11.4e.\n",Dynorm);

			//Triple linear system solver: J dPsi = -b, J Ddx = -Dx, J Ddy = -Dy
			status = BVP_LSsolver(&my_params, linsolvar, JacVal, JacRow, JacCol, dPsi, b, Ddx, Dx, Ddy, Dy, bnorm, Dxnorm, Dynorm, it, plotvar);
			if (status != 0)
			{
				printf("ERROR!!: In BVP_LSsolver.\n");
				break;
			}
			dPsinorm = norm(N,dPsi);
			printf("||dPsi|| = %11.4e.\n",dPsinorm);
			duDxnorm = norm(N,Ddx);
			printf("||dDx||  = %11.4e.\n",duDxnorm);
			duDynorm = norm(N,Ddy);
			printf("||dDy||  = %11.4e.\n",duDynorm);


			printf("Ratio comparison x: Need %11.4e < %11.4e.\n",duDxnorm,tol*unorm);
			printf("Ratio comparison y: Need %11.4e < %11.4e.\n",duDynorm,tol*unorm);

			//---------- Check Solution ----------
			if (bnorm <= tol) //First convergence condition
			{
				break; //Grid does not need to be resized
			}

			//---------- Check for grid resize----------
			if ((duDxnorm > tol*unorm) && (duDynorm > tol*unorm))
			{
				gridresizevar = 3;
				printf("\n----------------------------------------------\n");
				printf("Resizing x and y dimension of %i and %i points.\n",n,m);
			}
			else if (duDxnorm > tol*unorm)
			{
				gridresizevar = 1;
				printf("\n----------------------------------------------\n");
				printf("Resizing x dimension of %i points.\n",n);
			}
			else if (duDynorm > tol*unorm)
			{
				gridresizevar = 2;
				printf("\n----------------------------------------------\n");
				printf("Resizing y dimension of %i points.\n",m);
			}
			else
			{
				printf("Grid size is sufficient.\n");
				gridresizevar = 0;
			}

			//---------- Grid resize ----------
			if (gridresizevar != 0)
			{
				status = BVP_grid(&my_params, gridresizevar, &x, &y, &Psi, &GRIDdx, &GRIDdy, &Ddx, &Ddy, tol*unorm);
				if (status == -1)
				{
					printf("ERROR!!: In BVP_grid.\n");
					break;
				}
				n = my_params.nparam;
		    m = my_params.mparam;
		    N = my_params.Nparam;
				LSimax = my_params.LSimaxparam;
				nnzJacCount = my_params.nnzparam;

				status = BVP_realloc(&my_params, &b, &Dx, &Dy, &dPsi, &Ddx, &Ddy, &JacVal, &JacRow, &JacCol);

			}
			if (gridit++ >= griditmax-1)
			{
				printf("ERROR!!: BVP grid resize is not converging.\n");
				break;
			}


		} while (((duDxnorm > tol*unorm) || (duDynorm > tol*unorm)) && (status == 0));
		if (gridit == griditmax)
		{
			status = -3;
			break;
		}
		if (status != 0)
		{
			break;
		}


		//---------- Output current system ----------
		if (plotvar == 1)
		{
			status = BVP_out(&my_params, x, y, Psi, JacVal, JacRow, JacCol, dPsi, b, Dx, Dy, it); //Current iteration
		}

		//---------- Check Solution ----------
		if (bnorm <= tol) //First convergence condition
		{
			printf("\nSUCCESS!! ||b|| = %11.4e <= tol = %11.4e.\n\n",bnorm, tol);
			//printf("BVP_phys eval\n");
			//status = BVP_phys(&my_params, x, y, Psi, Jac, b, Dx, Dy);
			status = BVP_out(&my_params, x, y, Psi, JacVal, JacRow, JacCol, dPsi, b, Dx, Dy, -3); //Converged solution in different directory
			break;
		}

		//---------- Check convergence & Update Solution ----------
		status = BVP_conv(&my_params, x, y, Psi, JacVal, JacRow, JacCol, b, Dx, Dy, bnorm, &w, dPsi);
		if (status == -1)
		{
			printf("ERROR!!: In BVP_conv.\n");
			break;
		}
		if (status == -2)
		{
			printf("BVP not converging. w < %.4f.\n",wmin);
			break;
		}

	} //::::::::::::::::::::::::: END OF MAIN LOOP :::::::::::::::::::::::::
	if (it == bvpimax+1)
	{
		printf("\n\nBVP FAILURE! Reached max iteration %3i.\n\n\n",it-1);
	}
	if (status == -1)
	{
		printf("\n\nBVP FAILURE! ERROR in evaluation at %3i. Check code.\n\n\n",it);
	}
	if (status == -2)
	{
		printf("\n\nBVP FAILURE! Not converging at %3i.\n\n\n",it);
	}
	if (status == -3)
	{
		printf("\n\nBVP FAILURE! Grid not converging at %3i.\n\n\n",it);
	}


	//Cleanup
	free(x);
	free(y);
	free(Psi);
	free(GRIDdx);
	free(GRIDdy);
	free(b);
	free(Dx);
	free(Dy);
	free(dPsi);
	free(Ddx);
	free(Ddy);

	free(JacVal);
	free(JacRow);
	free(JacCol);


	return 0;
}



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::: Function Definitions :::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double norm(const int N, double v[])
{
	int i;
	double temp = 0.0;
	//Euclidean Norm
	/*
	for (i = 0; i < N; ++i)
	{
		temp += v[i]*v[i];
	}
	temp = sqrt(temp);
	*/
	//Max Norm
	for (i = 0; i < N; ++i)
	{
		if (temp < fabs(v[i]))
		{
			temp = fabs(v[i]);
		}
	}

	return temp;
}

int BVP_count(struct param_type *params)
{
	int i,j,k,l, status;
	const int n = params->nparam;
	const int m = params->mparam;
	const int p = params->pparam;
	const int N = params->Nparam;
	const int r = params->rparam;

	int nnztemp = 0;

	int dFEx  = dFExout();
	int dFEy  = dFEyout();
	int dBCXx = dBCXxout();
	int dBCXy = dBCXyout();
	int dBCYx = dBCYxout();
	int dBCYy = dBCYyout();

	//-------------------- Boundaries --------------------

	nnztemp += m*(dBCXx + dBCXy)*(r+1);
	nnztemp += m*dBCXx; //one-sided x
	nnztemp += r*dBCXy; //one-sided y

	nnztemp += (n-2)*(dBCYx + dBCYy)*(r+1); //The -2 is because X Boundary is applied on corners of grid
	nnztemp += (n-2)*dBCYy; //one-sided x
	nnztemp += (r-2)*dBCYx; //one-sided y

	//-------------------- Residual --------------------

	nnztemp += (n-2)*(m-2)*(dFEx + dFEy)*(r+1); //Diagonals (double counts d/dx and d/dy diagonals for storage)
	nnztemp += (r-2)*(m-2)*dFEx; //one-sided x
	nnztemp += (n-2)*(r-2)*dFEy; //one-sided y


	return nnztemp;
}

int BVP_COOmodtoCOO(struct param_type *params, int *nnz, double **Val, int **arr1, int **arr2)
{
	int i,j,k,l, status;
	const int N = (*params).Nparam;
	int nnztemp = *nnz;

	status = COOinsertionSort(&nnztemp, Val, arr1, arr2);

	/*
	//Realloc sparse Jac
	*Val  = realloc(*Val, (size_t) nnztemp * sizeof(double));
	if (*Val == NULL) {
		printf("ERROR in realloc!!\n");
		return -1;
	}
	*arr1  = realloc(*arr1, (size_t) nnztemp * sizeof(int));
	if (*arr1 == NULL) {
		printf("ERROR in realloc!!\n");
		return -1;
	}
	*arr2  = realloc(*arr2, (size_t) nnztemp * sizeof(int));
	if (*arr2 == NULL) {
		printf("ERROR in realloc!!\n");
		return -1;
	}
	*/
	*nnz = nnztemp;


	return 0;
}

/* Function to sort a sparse matrix in COO format using insertion sort*/
//Modified from https://www.geeksforgeeks.org/insertion-sort/
int COOinsertionSort(int *nnz, double **Val, int **arr1, int **arr2)
{
	int i,j, key1, key2;
	double keyval;
	int tempnnz, tempi;

	//Sort arr1
	for (i = 1; i < *nnz; ++i)
	{
		key1   = (*arr1)[i];
		key2   = (*arr2)[i];
		keyval = (*Val)[i];
		j = i - 1;

		//Move elements of arr1[0..i-1], that are greater than key, to one position ahead of their current position
		while (j >= 0 && (*arr1)[j] > key1)
		{
			(*arr1)[j+1] = (*arr1)[j];
			(*arr2)[j+1] = (*arr2)[j];
			(*Val)[j+1]  = (*Val)[j];
			--j;
		}
		(*arr1)[j+1] = key1;
		(*arr2)[j+1] = key2;
		(*Val)[j+1]  = keyval;
	}

	//Sort arr2
	tempnnz = 0;
	tempi = 0;

	while (tempi < *nnz-1)
	{
		//Find index of next row
		while (((*arr1)[tempi] != ((*arr1)[tempnnz] + 1)) && (tempi < *nnz-1) )
		{
			++tempi;
		}

		//Sort arr2 from tempnnz to tempi-1
		for (i = tempnnz; i < tempi; ++i)
		{
			key1   = (*arr1)[i];
			key2   = (*arr2)[i];
			keyval = (*Val)[i];
			j = i - 1;

			while (j >= tempnnz && (*arr2)[j] > key2)
			{
				(*arr1)[j+1] = (*arr1)[j];
				(*arr2)[j+1] = (*arr2)[j];
				(*Val)[j+1]  = (*Val)[j];
				--j;
			}
			(*arr1)[j+1] = key1;
			(*arr2)[j+1] = key2;
			(*Val)[j+1]  = keyval;
		}

		tempnnz = tempi;
	}

	//Combine duplicate matrix elements
	tempnnz = 0;
	for (i = 0; i < *nnz-1; ++i)
	{
		if ( ((*arr1)[i] == (*arr1)[i+1]) && ((*arr2)[i] == (*arr2)[i+1]) )
		{
			(*Val)[tempnnz] += (*Val)[i+1];
		}
		else
		{
			++tempnnz;
			(*Val)[tempnnz]  = (*Val)[i+1];
			(*arr1)[tempnnz] = (*arr1)[i+1];
			(*arr2)[tempnnz] = (*arr2)[i+1];
		}

	}
	++tempnnz;
	*nnz = tempnnz;


	return 0;
}

int BVP_COOtoCSR(struct param_type *params, int *nnz, double **Val, int **arr1, int **arr2, double **sa, unsigned long **ija)
{
	int i,j,k,l, status;
	unsigned long il,jl;
	const int N = (*params).Nparam;
	int nnztemp = *nnz;
	int nnzCSR = *nnz+2;
	int itemp;

	k = N + 2;
	printf("k = %i.\n",k);
	(*ija)[1] = (unsigned long) k;
	//(*ija)[itemp-1] = nnzCSR;

	for (i = 0; i < nnztemp; ++i)
	{
		if ((*arr1)[i] == (*arr2)[i])
		{
			(*sa)[(*arr1)[i]+1] = (*Val)[i];
		}
	}

	i = 0;
	while (i < nnztemp)
	{
		itemp = i + 1;

		if (itemp < (nnztemp-1))
		{
			while (((*arr1)[itemp] == (*arr1)[itemp-1]) && (itemp < (nnztemp-1)))
			{
				++itemp;
			}
		}
		if (itemp == (nnztemp-1))
		{
			if ((*arr1)[itemp] == (*arr1)[itemp-1])
			{
				++itemp;
			}
		}
		//printf("i = %i. itemp = %i.\n",i,itemp);

		for (j = i; j < itemp; ++j)
		{
			if ((*arr1)[j] != (*arr2)[j])
			{
				(*sa)[k] = (*Val)[j];
				(*ija)[k] = (unsigned long) ((*arr2)[j] + 1);
				++k;
			}
		}
		(*ija)[(*arr1)[i] + 2] = k;
		i = itemp;
	}


	/*
	printf("CSR format.\n");
	for (i = 0; i < nnztemp+2; ++i)
	{
		printf("sa[%2i]  = %11.4e. ija = %lu.\n",i,(*sa)[i],(*ija)[i]);
	}
	*/

	return 0;
}

int BVP_conv(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[], double bnorm, double *w, double dPsi[])
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
	const double wmin = (*params).wminparam;
	const double wmax = (*params).wmaxparam;
	const double wshrink = (*params).wshrinkparam;
	const double wgrow = (*params).wgrowparam;

	int nnzJacTemp;
	double bnormtemp;
	clock_t begin, end;
	double time_spent;

	double *Psitemp;
	Psitemp = (double*)calloc(N, sizeof(double));

	do
	{
		for (i = 0; i < N; ++i)
		{
			Psitemp[i] = Psi[i] + *w*dPsi[i];
		}

		//Linear system evaluation
		begin = clock();
		nnzJacTemp = BVP_GENsys(params, x, y, Psitemp, JacVal, JacRow, JacCol, bvec, Dxvec, Dyvec);
		end = clock();
		time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

		bnormtemp = norm(N, bvec);
		printf("New btemp = %.3e.\n",bnormtemp);
		*w *= wshrink;
		if (*w < wmin)
		{
			//Exiting
			free(Psitemp);
			return -2;
		}
	} while (bnormtemp > bnorm);

	//---------- Update Solution ----------
	*w /= wshrink; //undo decrement of do...while loop
	printf("Converged with ||btemp||    = %11.4e and w = %.3f.\n", bnormtemp, *w);
	for (i = 0; i < N; ++i)
	{
		Psi[i] += *w*dPsi[i];
	}
	*w *= wgrow; //Grow w upon successful convergence
	if (*w > wmax)
	{
		*w = wmax;
	}


	//Cleanup
	free(Psitemp);


	return 0;
}
