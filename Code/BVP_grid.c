#include "BVP_header.h"

//----- List of defined functions -----
/*
BVP_grid
BVP_resizecalc
BVP_gridsize
BVP_realloc
*/

int BVP_grid(struct param_type *params, const int gridresizevar, double **x, double **y, double **Psi, double **GRIDdxvec, double **GRIDdyvec, double **Ddxvec, double **Ddyvec, double GRIDtol)
{
	int i,j,k,l, status;
  const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;

  //----- Calloc BEGIN -----
  double *dxnew;
	dxnew = (double*)calloc(n, sizeof(double));
  double *dynew;
	dynew = (double*)calloc(m, sizeof(double));
  //----- Calloc END -----

  status = BVP_resizecalc(params, gridresizevar, GRIDtol, x, y, GRIDdxvec, GRIDdyvec, Ddxvec, Ddyvec, &dxnew, &dynew);
  status = BVP_gridsize(params, gridresizevar, x, y, Psi, GRIDdxvec, GRIDdyvec, &dxnew, &dynew);

  //Cleanup
  free(dxnew);
  free(dynew);


  return 0;
}

int BVP_resizecalc(struct param_type *params, const int gridresizevar, const double GRIDtol, double **x, double **y, double **dxold, double **dyold, double **Ddx, double **Ddy, double **dxnew, double **dynew)
{
  const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;
	int i,j,l,k;
	double duDtemp, errorratio;

	//---------- x dimension ----------
  if (gridresizevar != 2)
  {
    for (j = 0; j < n; ++j)
  	{
  		duDtemp = 0.0;
  		for (i = 0; i < m; ++i)
  		{
  			for (l = 0; l < p; ++l)
  			{
  				k = (i*n+j)*p;
  				duDtemp = MAX(duDtemp,fabs((*Ddx)[k +l]));
  			}
  		}
  		if (duDtemp <= 0.1*GRIDtol)
  		{
  			errorratio = 1.0;
  		}
  		else
  		{
  			errorratio = (1.0/3.0*GRIDtol)/(duDtemp);
  		}
  		(*dxnew)[j] = pow(errorratio,1.0/((double) r))*(*dxold)[j];
      //(*dxnew)[j] = ((*dxold)[j])/3.0*2.0;
  		//printf("duDtemp = %11.4e. errorratio = %11.4e. dxold[%i] = %11.4e. dxnew[%i] = %11.4e.\n",duDtemp,errorratio,j,(*dxold)[j],j,(*dxnew)[j]);
  	}
  }

	//---------- y dimension ----------
  if (gridresizevar != 1)
  {
  	for (i = 0; i < m; ++i)
  	{
  		duDtemp = 0.0;
  		for (j = 0; j < n; ++j)
  		{
  			for (l = 0; l < p; ++l)
  			{
  				k = (i*n+j)*p;
  				duDtemp = MAX(duDtemp,fabs((*Ddy)[k +l]));
  			}
  		}
  		if (duDtemp <= 0.1*GRIDtol)
  		{
  			errorratio = 1.0;
  		}
  		else
  		{
  			errorratio = (1.0/3.0*GRIDtol)/(duDtemp);
  		}
  		(*dynew)[i] = pow(errorratio,1.0/((double) r))*(*dyold)[i];
      //(*dynew)[i] = ((*dyold)[i])/3.0*2.0;
  		//printf("duDtemp = %11.4e. errorratio = %11.4e. dyold[%i] = %11.4e. dynew[%i] = %11.4e.\n",duDtemp,errorratio,i,(*dyold)[i],i,(*dynew)[i]);
  	}
  }

  if (gridresizevar == 1)
  {
    for (i = 0; i < m; ++i)
    {
      (*dynew)[i] = (*dyold)[i];
    }
  }
  if (gridresizevar == 2)
  {
    for (j = 0; j < n; ++j)
    {
      (*dxnew)[j] = (*dxold)[j];
    }
  }


	return 0;
}

int BVP_gridsize(struct param_type *params, const int gridresizevar, double **x, double **y, double **Psi, double **dxold, double **dyold, double **dxnew, double **dynew)
{
  int i,j,k,l, status;
  int n = (*params).nparam;
	int m = (*params).mparam;
	const int p = (*params).pparam;
	int N = (*params).Nparam;
	const int r = (*params).rparam;
  int nnzJacCount = (*params).nnzparam;
	int pn;
  double xk, yk;
	double dxtemp, dytemp, vi;
	double x0, x1;
	double y0, y1;
	double temp, temp2;
	int nmax, mmax;
  int ntemp, mtemp, Ntemp;
  int nold = n;
  int mold = m;
  int Nold = N;

  //----- Calloc BEGIN -----
	double *xtemp;
	xtemp = (double*)calloc(nold, sizeof(double));//(*x) Right
	double *ytemp;
	ytemp = (double*)calloc(mold, sizeof(double));//(*x) Right
	double *psitemp;
	psitemp = (double*)calloc(Nold, sizeof(double));
  //----- Calloc END -----

  //Temporary storage
	for (j = 0; j < nold; ++j)
	{
		xtemp[j] = (*x)[j];
		//printf("xtemp[%3i] = %7.4f.\n",j,xtemp[j]);
	}
	for (i = 0; i < mold; ++i)
	{
		ytemp[i] = (*y)[i];
		//printf("ytemp[%3i] = %7.4f.\n",i,ytemp[i]);
	}
	for (k = 0; k < Nold; ++k)
	{
		psitemp[k] = (*Psi)[k];
		//printf("psitemp[%3i] = %11.4e.\n",k,psitemp[k]);
	}


  //---------- New grid size estimation ----------
  //x dimension
  if (gridresizevar != 2)
  {
  	nmax = 0;
  	temp = 0.0;
  	temp2 = 0.0;
  	for (i = 0; i < nold; ++i)
  	{
  		if (temp2 < fabs((*dxnew)[i]))
  		{
  			temp2 = fabs((*dxnew)[i]);
  			nmax = i;
  		}
  		//printf("%3i: x = %6.3f, dxold = %6.3f, dxnew = %7.4f.\n",i, (*x)[i], (*dxold)[i], (*dxnew)[i]);
  		if ((i == 0) || (i == nold-1))
  		{
  			temp += (*dxnew)[i]/2.0;
  		}
  		else
  		{
  			temp += (*dxnew)[i];
  		}
  	}
  	if ((nmax == nold-1) || (nmax == 0))
  	{
  		nmax = (int) round(nold/2);
  	}
    ntemp = (int) (((double) nold)/(temp/((*x)[nold-1]))); //Safety factor of 2 to make temp array larger than needed
  	//printf("total step distance = %11.4e. Estimated new points = %i.\n",temp,ntemp);
  	//printf("max stepsize is %7.4f at %i.\n",temp2,nmax);

    //---------- Allocate and calculate new stencil ----------
    //----- Calloc BEGIN -----
    double *xL;
  	xL = (double*)calloc(ntemp, sizeof(double));//x Left
  	double *dxL;
  	dxL = (double*)calloc(ntemp, sizeof(double));//Delta x Left
  	double *xR;
  	xR = (double*)calloc(ntemp, sizeof(double));//x Right
  	double *dxR;
  	dxR = (double*)calloc(ntemp, sizeof(double));//Delta x Right
    //----- Calloc END -----

    //--------------- x grid construction ---------------
  	//Left to middle
  	j = 0;
  	xL[j] = xtemp[0];
  	dxL[j] = (*dxnew)[0];

  	for (i = 0; i < nmax; ++i)
  	{
  		x0 = xL[j] + dxL[j]/2.0;
  		dxtemp = (*dxnew)[i];
  		x1 = xtemp[i] + (*dxold)[i]/2.0;
  		while (x0 + dxtemp/2.0 < x1)
  		{
  			xL[j+1] = x0 + dxtemp/2.0;
  			dxL[j+1] = dxtemp;
  			x0 = xL[j+1] + dxtemp/2.0;
  			++j;
  		}
  	}

  	//Right to middle
  	k = 0;
  	xR[k] = xtemp[nold-1];
  	dxR[k] = (*dxnew)[nold-1];
  	for (i = nold-1; i > nmax; --i)
  	{
  		x1 = xR[k] - dxR[k]/2.0;
  		dxtemp = (*dxnew)[i];
  		x0 = xtemp[i] - (*dxold)[i]/2.0;
  		while (x1 - dxtemp/2.0 > x0)
  		{
  			xR[k+1] = x1 - dxtemp/2.0;
  			dxR[k+1] = dxtemp;
  			x1 = xR[k+1] - dxtemp/2.0;
  			++k;
  		}
  	}

  	//Center
  	dxtemp = (*dxnew)[nmax];
  	x0 = xL[j] + dxL[j]/2.0;
  	x1 = xR[k] - dxR[k]/2.0;
  	vi = (x1 - x0)/dxtemp;
  	while (vi > 2.0)
  	{
  		//Left
  		xL[j+1] = x0 + dxtemp/2.0;
  		dxL[j+1] = dxtemp;
  		x0 = xL[j+1] + dxtemp/2.0;
  		++j;
  		//Right
  		xR[k+1] = x1 - dxtemp/2.0;
  		dxR[k+1] = dxtemp;
  		x1 = xR[k+1] - dxtemp/2.0;
  		++k;
  		vi = vi - 2.0;
  	}
  	if (vi > 1.0)
  	{
  		dxtemp = (x1-x0)/2.0;
  		//Left
  		xL[j+1] = x0 + dxtemp/2.0;
  		dxL[j+1] = dxtemp;
  		x0 = xL[j+1] + dxtemp/2.0;
  		++j;
  		//Right
  		xR[k+1] = x1 - dxtemp/2.0;
  		dxR[k+1] = dxtemp;
  		x1 = xR[k+1] - dxtemp/2.0;
  		++k;
  	}
  	else
  	{
  		dxtemp = (x1 - x0);
  		//Left
  		xL[j+1] = x0 + dxtemp/2.0;
  		dxL[j+1] = dxtemp;
  		x0 = xL[j+1] + dxtemp/2.0;
  		++j;
  	}
  	ntemp = j+k+2;
    //printf("new n = %i.\n",ntemp);

  	*x = (double*)realloc(*x, (size_t) ntemp * sizeof(double));
    if (*x == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
    *dxold = (double*)realloc(*dxold, (size_t) ntemp * sizeof(double));
    if (*dxold == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
  	l = 0;
  	for (i = 0; i <= j; ++i)
  	{
  		//printf("xL  = %7.4f.\n",xL[i]);
      //printf("dxL = %7.4f.\n",dxL[i]);
  		(*x)[l] = xL[i];
      (*dxold)[l] = dxL[i];
  		++l;
  	}
  	for (i = k; i >= 0; --i)
  	{
  		//printf("xR  = %7.4f.\n",xR[i]);
      //printf("dxR = %7.4f.\n",dxR[i]);
  		(*x)[l] = xR[i];
      (*dxold)[l] = dxR[i];
  		++l;
  	}

    //Cleanup
    free(xL);
  	free(dxL);
  	free(xR);
  	free(dxR);
  }//End x dimension


  //y dimension
  if (gridresizevar != 1)
  {
  	mmax = 0;
  	temp = 0.0;
  	temp2 = 0.0;
  	for (j = 0; j < mold; ++j)
  	{
  		if (temp2 < fabs((*dynew)[j]))
  		{
  			temp2 = fabs((*dynew)[j]);
  			mmax = j;
  		}
  		//printf("%3i: y = %6.3f, dyold = %6.3f, dynew = %7.4f.\n",j, (*y)[j], (*dyold)[j], (*dynew)[j]);
  		if (j == 0 || j == mold-1)
  		{
  			temp += (*dynew)[j]/2.0;
  		}
  		else
  		{
  			temp += (*dynew)[j];
  		}
  	}
  	if ((mmax == mold-1) || (mmax == 0))
  	{
  		mmax = (int) round(mold/2);
  	}
  	mtemp = (int) (((double) mold)/(temp/((*y)[mold-1]))); //Safety factor of 2 to make temp array larger than needed
  	//printf("total step distance = %11.4e. Estimated new points = %i.\n",temp,mtemp);
  	//printf("max stepsize is %7.4f at %i.\n",temp2,mmax);

    //---------- Allocate and calculate new stencil ----------
    //----- Calloc BEGIN -----
  	double *yL;
  	yL = (double*)calloc(mtemp, sizeof(double));//y Left
  	double *dyL;
  	dyL = (double*)calloc(mtemp, sizeof(double));//Delta y Left
  	double *yR;
  	yR = (double*)calloc(mtemp, sizeof(double));//y Right
  	double *dyR;
  	dyR = (double*)calloc(mtemp, sizeof(double));//Delta y Right
    //----- Calloc END -----

  	//--------------- y grid construction ---------------
  	//Left to middle
  	j = 0;
  	yL[j] = ytemp[0];
  	dyL[j] = (*dynew)[0];

  	for (i = 0; i < mmax; ++i)
  	{
  		y0 = yL[j] + dyL[j]/2.0;
  		dytemp = (*dynew)[i];
  		y1 = ytemp[i] + (*dyold)[i]/2.0;
  		while (y0 + dytemp/2.0 < y1)
  		{
  			yL[j+1] = y0 + dytemp/2.0;
  			dyL[j+1] = dytemp;
  			y0 = yL[j+1] + dytemp/2.0;
  			++j;
  		}
  	}

  	//Right to middle
  	k = 0;
  	yR[k] = ytemp[mold-1];
  	dyR[k] = (*dynew)[mold-1];
  	for (i = mold-1; i > mmax; --i)
  	{
  		y1 = yR[k] - dyR[k]/2.0;
  		dytemp = (*dynew)[i];
  		y0 = ytemp[i] - (*dyold)[i]/2.0;
  		while (y1 - dytemp/2.0 > y0)
  		{
  			yR[k+1] = y1 - dytemp/2.0;
  			dyR[k+1] = dytemp;
  			y1 = yR[k+1] - dytemp/2.0;
  			++k;
  		}
  	}

  	//Center
  	dytemp = (*dynew)[mmax];
  	y0 = yL[j] + dyL[j]/2.0;
  	y1 = yR[k] - dyR[k]/2.0;
  	vi = (y1 - y0)/dytemp;
  	while (vi > 2.0)
  	{
  		//Left
  		yL[j+1] = y0 + dytemp/2.0;
  		dyL[j+1] = dytemp;
  		y0 = yL[j+1] + dytemp/2.0;
  		++j;
  		//Right
  		yR[k+1] = y1 - dytemp/2.0;
  		dyR[k+1] = dytemp;
  		y1 = yR[k+1] - dytemp/2.0;
  		++k;
  		vi = vi - 2.0;
  	}
  	if (vi > 1.0)
  	{
  		dytemp = (y1-y0)/2.0;
  		//Left
  		yL[j+1] = y0 + dytemp/2.0;
  		dyL[j+1] = dytemp;
  		y0 = yL[j+1] + dytemp/2.0;
  		++j;
  		//Right
  		yR[k+1] = y1 - dytemp/2.0;
  		dyR[k+1] = dytemp;
  		y1 = yR[k+1] - dytemp/2.0;
  		++k;
  	}
  	else
  	{
  		dytemp = (y1 - y0);
  		//Left
  		yL[j+1] = y0 + dytemp/2.0;
  		dyL[j+1] = dytemp;
  		y0 = yL[j+1] + dytemp/2.0;
  		++j;
  	}
  	mtemp = j+k+2;
    //printf("new m = %i.\n",mtemp);

  	*y = (double*)realloc(*y, (size_t) mtemp * sizeof(double));
    if (*y == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
    *dyold = (double*)realloc(*dyold, (size_t) mtemp * sizeof(double));
    if (*dyold == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
  	l = 0;
  	for (i = 0; i <= j; ++i)
  	{
  		//printf("yL  = %7.4f.\n",yL[i]);
      //printf("dyL = %7.4f.\n",dyL[i]);
  		(*y)[l] = yL[i];
      (*dyold)[l] = dyL[i];
  		++l;
  	}
  	for (i = k; i >= 0; --i)
  	{
  		//printf("yR  = %7.4f.\n",yR[i]);
      //printf("dyR = %7.4f.\n",dyR[i]);
  		(*y)[l] = yR[i];
      (*dyold)[l] = dyR[i];
  		++l;
  	}

    //Cleanup
    free(yL);
  	free(dyL);
  	free(yR);
  	free(dyR);
  }//End y dimension


  //---------- Interpolate solution ----------
  if (gridresizevar == 1) //Only x adjust
  {
    mtemp = mold;
    //New x-grid
	  Ntemp = ntemp*mold*p;
    printf("New points: n = %i, m = same. N = %i.\n",ntemp,Ntemp);
    printf("----------------------------------------------\n\n");
    *Psi = (double*)realloc(*Psi, (size_t) Ntemp * sizeof(double));
    if (*Psi == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
		for (k = 0; k < Ntemp; ++k)
		{
			(*Psi)[k] = 0.0;
		}
		//Interpolate each Psi to new x-grid
		status = BVP_interp(params, xtemp, ytemp, psitemp, ntemp, mold, *x, ytemp, *Psi);
    (*params).nparam = ntemp;
  }
  else if (gridresizevar == 2) //Only y adjust
  {
    ntemp = nold;
    //New y-grid
		Ntemp = nold*mtemp*p;
		printf("New points: n = same, m = %i. N = %i.\n",mtemp,Ntemp);
    printf("----------------------------------------------\n\n");
		*Psi = (double*)realloc(*Psi, (size_t) Ntemp * sizeof(double));
    if (*Psi == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
    for (k = 0; k < Ntemp; ++k)
		{
			(*Psi)[k] = 0.0;
		}
		//Interpolate each Psi to new y-grid
    status = BVP_interp(params, xtemp, ytemp, psitemp, nold, mtemp, xtemp, *y, *Psi);
    (*params).mparam = mtemp;
  }
  else if (gridresizevar == 3) //With both dimensions changing need to resize 1 at a time.
  {
    //New x-grid first
    Ntemp = ntemp*mold*p;
    printf("New points: n = %i, m = same. N = %i.\n",ntemp,Ntemp);
		*Psi = (double*)realloc(*Psi, (size_t) Ntemp * sizeof(double));
    if (*Psi == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
    for (k = 0; k < Ntemp; ++k)
		{
			(*Psi)[k] = 0.0;
		}
		//Interpolate each Psi to new x-grid
		status = BVP_interp(params, xtemp, ytemp, psitemp, ntemp, mold, *x, ytemp, *Psi);
    (*params).nparam = ntemp;

		//Set psitemp to new x-grid
		psitemp = (double*)realloc(psitemp, (size_t) Ntemp * sizeof(double));
    if (psitemp == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
    for (k = 0; k < Ntemp; ++k)
		{
      psitemp[k] = (*Psi)[k];
		}

    //New y-grid next
		Ntemp = ntemp*mtemp*p;
		printf("New points: n = same, m = %i. N = %i.\n",mtemp,Ntemp);
    printf("----------------------------------------------\n\n");
		*Psi = (double*)realloc(*Psi, (size_t) Ntemp * sizeof(double));
    if (*Psi == NULL) {
      printf("ERROR in realloc!!\n");
      return -1;
    }
    for (k = 0; k < Ntemp; ++k)
		{
			(*Psi)[k] = 0.0;
		}
		//Interpolate each Psi to new y-grid
    status = BVP_interp(params, *x, ytemp, psitemp, ntemp, mtemp, *x, *y, *Psi);
    (*params).mparam = mtemp;
  }
  else
  {
    printf("ERROR. gridresizevar invalid.\n");
    return -1;
  }
  //Set parameters to new grid
  (*params).Nparam = Ntemp;
  nnzJacCount = BVP_count(params);
  (*params).LSimaxparam = N;
  (*params).nnzparam = nnzJacCount;

  //Cleanup
	free(xtemp);
	free(ytemp);
	free(psitemp);


	return 0;
}

int BVP_realloc(struct param_type *params, double **b, double **Dx, double **Dy, double **dPsi, double **Ddx, double **Ddy, double **JacVal, int **JacRow, int **JacCol)
{
  int i,j,k,l, status;
	const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;
  const int nnzJac = (*params).nnzparam;

  *b  = (double*)realloc(*b, (size_t) N * sizeof(double));
  if (*b == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *Dx  = (double*)realloc(*Dx, (size_t) N * sizeof(double));
  if (*Dx == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *Dy  = (double*)realloc(*Dy, (size_t) N * sizeof(double));
  if (*Dy == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *dPsi  = (double*)realloc(*dPsi, (size_t) N * sizeof(double));
  if (*dPsi == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *Ddx  = (double*)realloc(*Ddx, (size_t) N * sizeof(double));
  if (*Ddx == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *Ddy  = (double*)realloc(*Ddy, (size_t) N * sizeof(double));
  if (*Ddy == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }

  *JacVal  = (double*)realloc(*JacVal, (size_t) nnzJac * sizeof(double));
  if (*JacVal == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *JacRow  = (int*)realloc(*JacRow, (size_t) nnzJac * sizeof(int));
  if (*JacRow == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }
  *JacCol  = (int*)realloc(*JacCol, (size_t) nnzJac * sizeof(int));
  if (*JacCol == NULL) {
    printf("ERROR in realloc!!\n");
    return -1;
  }

  //Initialize memory
  for (k = 0; k < N; ++k)
  {
    (*b)[k] = 0.0;
    (*Dx)[k] = 0.0;
    (*Dy)[k] = 0.0;
    (*dPsi)[k] = 0.0;
    (*Ddx)[k] = 0.0;
    (*Ddy)[k] = 0.0;
  }
  for (k = 0; k < nnzJac; ++k)
  {
    (*JacVal)[k] = 0.0;
    (*JacRow)[k] = 0;
    (*JacCol)[k] = 0;
  }


  return 0;
}
