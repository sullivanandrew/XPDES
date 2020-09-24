#include "BVP_header.h"

int BVP_GENIC(struct param_type *params, double x[], double y[], double Psi[], double GRIDdx[], double GRIDdy[])
{
	int i,j,k,l, status;
	const int n = (*params).nparam;
	const int m = (*params).mparam;
	const int p = (*params).pparam;
	const int N = (*params).Nparam;
	const int r = (*params).rparam;
	const int nnzJac = (*params).nnzparam;
	const double tol = (*params).tolparam;
	const double alpha = (*params).alphaparam;
	const double beta = (*params).betaparam;
	const double r_H = (*params).r_Hparam;
	const double M = (*params).Mparam;
	const double chi = (*params).chiparam;
	const double ICalpha = (*params).ICalphaparam;
	const double ICchi = (*params).ICchiparam;

	int pn;
	double xk, yk;

	FILE *infile;
	double temp0,temp1,temp2,temp3,temp4,temp5;
	double *xGrid;
	xGrid = (double*)calloc(n*m, sizeof(double));
	if (!xGrid)
		status = -1;
  double *yGrid;
  yGrid = (double*)calloc(n*m, sizeof(double));
	if (!yGrid)
		status = -1;


	if (ICalpha < 0.0)
	{
		//Generate IC
		printf("Generating IC with delta = %.4f.\n",ICchi);
		for (i = 0; i < m; ++i)
		{
			y[i] = ((double) i/(m-1)) * (M_PI/2.0); //0 < y < M_PI/2
			yk = y[i];
			for (j = 0; j < n; ++j)
			{
				x[j] = ((double) j/(n-1));
				xk = x[j];
				k = (i*n+j)*p;

				//1 Dimensional equations
/*
				if (i == 0 || i == m-1)
				{
					pn = 0;
					//Psi[k +pn] = 1.0;
					//Psi[k +pn] = pow(xk,2.0);
					//Psi[k +pn] = 16.0/pow(1.0+xk,4.0);
					Psi[k +pn] = 1.0/pow(2.0-xk,2.0);
					//Psi[k +pn] = IC0out(xk, yk, r_H, M, alpha, beta, ICchi);

					pn = 1;
					//Psi[k +pn] = 16.0/pow(1.0+xk,4.0);
					Psi[k +pn] = pow(2.0-xk,2.0);
					//Psi[k +pn] = IC1out(xk, yk, r_H, M, alpha, beta, ICchi);

				}
				else
				{
					pn = 0;
					Psi[k +pn] = IC0out(xk, yk, r_H, M, alpha, beta, ICchi);

					pn = 1;
					Psi[k +pn] = IC1out(xk, yk, r_H, M, alpha, beta, ICchi);

				}
*/

				pn = 0;
				Psi[k +pn] = IC0out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 1;
				Psi[k +pn] = IC1out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 2;
				Psi[k +pn] = IC2out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 3;
				Psi[k +pn] = IC3out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 4;
				Psi[k +pn] = IC4out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 5;
				Psi[k +pn] = IC5out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 6;
				Psi[k +pn] = IC6out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 7;
				Psi[k +pn] = IC7out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 8;
				Psi[k +pn] = IC8out(xk, yk, r_H, M, alpha, beta, ICchi);

				pn = 9;
				Psi[k +pn] = IC9out(xk, yk, r_H, M, alpha, beta, ICchi);

			}
		}
	}
	else
	{
		char InputName[256];
		sprintf(InputName,"./Data/BVPoutLIN/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round((alpha-ICalpha)*1.0E4), (int) round(-log10(tol)), -3, (int) round(r_H*1.0E2), (int) round((chi-ICchi)*1.0E2), n, m);
		printf("Importing from = %s.\n",InputName);
		infile = fopen(InputName,"r");
		if(infile == NULL)
		{//Open file
			perror("Error opening file.\n");
		}//Open file


		//Import solution
	  for (i = 0; i < m; ++i)
		{
			for (j = 0; j < n; ++j)
			{
	      l = (i*n+j);
	      fscanf(infile, "%lf %lf",&temp0,&temp1);
	      xGrid[l] = temp0;
	      yGrid[l] = temp1;
	      k = (i*n+j)*p;
	      for (pn = 0; pn < p; ++pn)
				{
	        fscanf(infile, "%lf %lf %lf %lf",&temp2,&temp3,&temp4,&temp5);
	        Psi[k +pn] = temp2;
	        //printf("x[%4i] = %.2f, y[%4i] = %.2f. Psi_%i = %11.4e.\n",l,xGrid[l],l,yGrid[l],pn,Psi[k+pn]);
				}
			}
		}
		j = 0;
		for (i = 0; i < m; ++i)
		{
			k = (i*n+j);
			y[i] = yGrid[k];
			//printf("y[%2i] = %f.\n",i,y[i]);
		}
		i = 0;
		for (j = 0; j < n; ++j)
		{
			k = (i*n+j);
			x[j] = xGrid[k];
			//printf("x[%3i] = %f.\n",j,x[j]);
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


	//Cleanup
	free(xGrid);
	free(yGrid);


	return 0;
}
