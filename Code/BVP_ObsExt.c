#include "BVP_header.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include "../Funcs/Decl_PHYS.c"
#include "../Funcs/Decl_KERR.c"
//#include "../Funcs/Defn_PHYSout.c"

int main(int argc, char *argv[])
{
	int i,j,k,l, status;
	int pn;

	char InputName[256];
	char SolType[16];

	int SolTypeparam = atoi(argv[1]);
	if (SolTypeparam == 0)
	{
		strcpy(SolType, "KERR");
	}
	else if (SolTypeparam == 1)
	{
		strcpy(SolType, "LIN");
	}
	else if (SolTypeparam == 2)
	{
		strcpy(SolType, "EXP");
	}


	strcpy(InputName, "./Data/BVPout");
	strcat(InputName, SolType);
	strcat(InputName, "/BVPout_sols/");
	//strcat(InputName, "/BVPout_iterations/");
  strcat(InputName, argv[2]);
  strcat(InputName, ".dat");
  printf("Input File name: %s\n",InputName);

  const int p = atoi(argv[3]); //Number of field equations. Will always be twice the number of fields due to mixed derivatives
	const int n = atoi(argv[4]); //x grid size
	const int m = atoi(argv[5]); //theta grid size
	const int r = atoi(argv[6]); //Newton polynomial order
	const double tol = pow(10.0,atof(argv[7])); //absolute tolerance of solver
	const double alpha = atof(argv[8]); //coupling parameter alpha
	const double beta = 1.0; //coupling beta set to 1
	const double r_H = atof(argv[9]); //rH value. Note: this is technically \bar{r}_H = r_H/4 <- Schwarzschild
	const double chi = atof(argv[10]); //dimensionless spin parameter
	const double M = 2.0*r_H/sqrt(1.0-pow(chi,2.0)); //bare mass calculated from spin and rH
	const double Omega_H = sqrt(pow(M,2.0)-4.0*pow(r_H,2.0))/(2.0*M*(M+2.0*r_H));
	const double J0 = M*sqrt(pow(M,2.0) - 4.0*pow(r_H,2.0));
	double M0;

	const int N = n*m*p;

	printf("\nRunning BVP with tol = %8.1e,\n",tol);
	printf("Coupling value alpha = %.5f, beta = %.2f.\n",alpha,beta);
	printf("Black hole parameters are r_bar_H = %5.2f [Msun].\n",r_H);
	printf("M = %5.2f [Msun]. J = %5.2f [Msun^2].\n",M,J0);
	printf("Chi = a/M = %5.2f []. a = %5.2f [Msun].\n",chi,chi*M);
	printf("Derivative order r = %2i.\n",r);
	printf("Grid size nxm = %3i x %3i. Total pts = %i.\n",n,m,N);


  FILE *infile;
  double temp0,temp1,temp2,temp3,temp4,temp5;


  infile = fopen(InputName,"r");
	if(infile == NULL)
	{//Open file
		perror("Error opening file.\n");
	}//Open file


  status = 0;
  double *Psi;
	Psi = (double*)calloc(N, sizeof(double));
	if (!Psi)
		status = -1;
	double *Psibar;
	Psibar = (double*)calloc(N, sizeof(double));
	if (!Psibar)
		status = -1;
  double *xGrid;
	xGrid = (double*)calloc(n*m, sizeof(double));
	if (!xGrid)
		status = -1;
  double *yGrid;
  yGrid = (double*)calloc(n*m, sizeof(double));
	if (!yGrid)
		status = -1;
	double *x;
	x = (double*)calloc(n, sizeof(double));
	if (!x)
		status = -1;
	double *y;
	y = (double*)calloc(m, sizeof(double));
	if (!y)
		status = -1;


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
        Psibar[k +pn] = temp2;
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


	//Convert bar quantities to unbar quantities. f = x^2*fbar, m = x^2*mbar, l = x^2*lbar
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			k = (i*n+j)*p;
			pn = 0;
			Psi[k +pn] = pow(x[j],2.0)*Psibar[k +pn];
			pn = 1;
			Psi[k +pn] = pow(x[j],2.0)*Psibar[k +pn];
			pn = 2;
			Psi[k +pn] = pow(x[j],2.0)*Psibar[k +pn];
			pn = 3;
			Psi[k +pn] = Psibar[k +pn];
			pn = 4;
			Psi[k +pn] = Psibar[k +pn];

			pn = 5;
			Psi[k +pn] = pow(x[j],2.0)*Psibar[k +pn];
			pn = 6;
			Psi[k +pn] = pow(x[j],2.0)*Psibar[k +pn];
			pn = 7;
			Psi[k +pn] = pow(x[j],2.0)*Psibar[k +pn];
			pn = 8;
			Psi[k +pn] = Psibar[k +pn];
			pn = 9;
			Psi[k +pn] = Psibar[k +pn];
		}
	}


	const size_t slen = r+4;
	int a, b; //a for x, b for y(theta)
	int I, onexside, twoxside;
	int J, oneyside, twoyside;
	double xk, yk;
	const int dimx = 1;
	const int dimy = n;
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

	//New variables
	double rhoT;
	double xT;
	double gtt,gtph,gpsi;
	double Minf;
	double Dinf;
	double C1inf, Mtemp;
	double Jinf, ainf, chiinf;
	double zeta, dM, dMcoeff, dD, dJ, dJcoeff;
	double xISCO, xISCOSCHW, xISCOKERR;
	double xLR, xLRSCHW, xLRKERR;
	double rhoISCO, rhoISCOSCHW, rhoISCOKERR;
	double rISCOKERR, rLRKERR, rLRSCHW, rISCOSCHW;
	double rhoLR, rhoLRSCHW, rhoLRKERR;
	double rLRKERRsol, rhoLRKERRsol, xLRKERRsol;
	double rISCOKERRsol, rhoISCOKERRsol, xISCOKERRsol;
	double dxISCOcoeff, dxLRcoeff;
	double OmegaLR, OmegaISCO;

	double xk1,yk1;
	double chifKERR = 0.0;
	double chimKERR = 0.0;
	double chilKERR = 0.0;
	double chioKERR = 0.0;

	double chifSCHW = 0.0;
	double chimSCHW = 0.0;
	double chilSCHW = 0.0;
	double chioSCHW = 0.0;


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	//-------------------- Compare solution to Kerr --------------------

	for (i = 0; i < m-1; ++i)
  {
    for (j = 0; j < n-2; ++j) //(Data-expected)^2 at x=1 will be exactly zero by boundary conditions.
    {
			k = (i*n+j)*p;
      yk = y[i];
      yk1 = y[i+1];
      xk = x[j];
      xk1 = x[j+1];

			pn = 0;
			chifSCHW += (xk*pow(Psibar[(i*n+j)*p +pn] - IC0out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC0out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC0out(xk,yk1,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC0out(xk,yk1,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC0out(xk1,yk,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC0out(xk1,yk,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC0out(xk1,yk1,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC0out(xk1,yk1,r_H,2.0*r_H,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

			chifKERR += (xk*pow(Psibar[(i*n+j)*p +pn] - IC0out(xk,yk,r_H,M,0.0,1.0,0.0),2.0)/IC0out(xk,yk,r_H,M,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC0out(xk,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC0out(xk,yk1,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC0out(xk1,yk,r_H,M,0.0,1.0,0.0),2.0)/IC0out(xk1,yk,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC0out(xk1,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC0out(xk1,yk1,r_H,M,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

			pn = 1;
			chimSCHW += (xk*pow(Psibar[(i*n+j)*p +pn] - IC1out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC1out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC1out(xk,yk1,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC1out(xk,yk1,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC1out(xk1,yk,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC1out(xk1,yk,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC1out(xk1,yk1,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC1out(xk1,yk1,r_H,2.0*r_H,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

			chimKERR += (xk*pow(Psibar[(i*n+j)*p +pn] - IC1out(xk,yk,r_H,M,0.0,1.0,0.0),2.0)/IC1out(xk,yk,r_H,M,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC1out(xk,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC1out(xk,yk1,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC1out(xk1,yk,r_H,M,0.0,1.0,0.0),2.0)/IC1out(xk1,yk,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC1out(xk1,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC1out(xk1,yk1,r_H,M,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

			pn = 2;
			chilSCHW += (xk*pow(Psibar[(i*n+j)*p +pn] - IC2out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC2out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC2out(xk,yk1,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC2out(xk,yk1,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC2out(xk1,yk,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC2out(xk1,yk,r_H,2.0*r_H,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC2out(xk1,yk1,r_H,2.0*r_H,0.0,1.0,0.0),2.0)/IC2out(xk1,yk1,r_H,2.0*r_H,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

			chilKERR += (xk*pow(Psibar[(i*n+j)*p +pn] - IC2out(xk,yk,r_H,M,0.0,1.0,0.0),2.0)/IC2out(xk,yk,r_H,M,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC2out(xk,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC2out(xk,yk1,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC2out(xk1,yk,r_H,M,0.0,1.0,0.0),2.0)/IC2out(xk1,yk,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC2out(xk1,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC2out(xk1,yk1,r_H,M,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

			pn = 3;
			chioKERR += (xk*pow(Psibar[(i*n+j)*p +pn] - IC3out(xk,yk,r_H,M,0.0,1.0,0.0),2.0)/IC3out(xk,yk,r_H,M,0.0,1.0,0.0) +
			xk*pow(Psibar[((i+1)*n+j)*p +pn] - IC3out(xk,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC3out(xk,yk1,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[(i*n+(j+1))*p +pn] - IC3out(xk1,yk,r_H,M,0.0,1.0,0.0),2.0)/IC3out(xk1,yk,r_H,M,0.0,1.0,0.0) +
			xk1*pow(Psibar[((i+1)*n+(j+1))*p +pn] - IC3out(xk1,yk1,r_H,M,0.0,1.0,0.0),2.0)/IC3out(xk1,yk1,r_H,M,0.0,1.0,0.0))*(xk1-xk)*(yk1-yk)/4.0;

    }
  }
	chifSCHW /= M_PI/2.0;
	chifKERR /= M_PI/2.0;

	chimSCHW /= M_PI/2.0;
	chimKERR /= M_PI/2.0;

	chilSCHW /= M_PI/2.0;
	chilKERR /= M_PI/2.0;

	chioKERR /= M_PI/2.0;

	printf("chifSCHW = %11.4e. chifKERR = %11.4e.\n",chifSCHW,chifKERR);
	printf("chimSCHW = %11.4e. chimKERR = %11.4e.\n",chimSCHW,chimKERR);
	printf("chilSCHW = %11.4e. chilKERR = %11.4e.\n",chilSCHW,chilKERR);
	printf("			chioKERR = %11.4e.\n",chioKERR);



	//-------------------- Calculate mass and charge at sufficient infinity --------------------
	//rhoT = r_H/tol;
	rhoT = sqrt(r_H/tol);
	//rhoT = cbrt(r_H/tol);
	//rhoT = sqrt(r_H/tol) - alpha*cbrt(r_H/tol); //Sufficient infinity depends on coupling

	xT = (1.0-r_H/rhoT);

	printf("rho such that O(1/rho) term in g_{tt} = tol: rho =  %11.4e.\n",rhoT);
	printf("xT = %14.7e.\n",xT);
	xk = xT;
	for (j = 0; j < n; ++j)
	{
		if (x[j] >= xk)
		{
			J = j;
			break;
		}
	}

	//for (i = 1; i < m-1; ++i)
	i = m-1; //Evaluate on equator
	{
		{
			k = (i*n+j)*p;
			I = i, J = j;
			yk = y[I];

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

			//Mass and Charge calc
			gtt  =	gmetric0out (xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
			gtph =	gmetric1out (xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
			gpsi =	gmetric2out (xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);

			printf("f = %11.4e. m = %11.4e. l = %11.4e. omega = %11.4e. psi = %11.4e.\n",u0d00k, u1d00k,u2d00k,u3d00k,u4d00k);
			printf("gtt = %11.4e. gtph = %11.4e. gpsi = %11.4e.\n",gtt,gtph,gpsi);

			//O(1/r)
			//Minf = rhoT/2.0*(1.0+gtt);
			//Jinf = -rhoT/2.0/pow(sin(yk),2.0)*gtph;
			//Dinf = rhoT*gpsi;

			//O(1/r^2)
			Minf = rhoT/2.0*(1.0-sqrt(-1.0-2.0*gtt));
			//Jinf = gtph/(-2.0*pow(sin(yk),2.0)/rhoT)/(1.0 - Minf/rhoT);
			//ainf = Jinf/Minf;
			Dinf = rhoT*gpsi;

			//O(1/r^3)
			/*
			C1inf = pow(rhoT,2.0)*(u2d00k-1.0);
			printf("C1inf = %14.7e.\n",C1inf);
			Dinf = u4d00k/(-1.0/rhoT+C1inf/pow(rhoT,3.0));
			printf("Dinf  = %14.7e.\n",Dinf);
			printf("check = %14.7e.\n",-pow(Dinf,2.0)/8.0-2.0*pow(r_H,2.0));
			//Minf = sqrt(pow(rhoT,2.0)*(1.0-u1d00k)-C1inf-pow(Dinf,2.0)/4.0);
			Mtemp = (-6.0*rhoT*pow(Dinf,2.0)+208.0*pow(rhoT,3.0)+288.0*pow(r_H,2.0)*rhoT-432.0*u0d00k*pow(rhoT,3.0)+sqrt((186624.0*pow(u0d00k,2.0)-179712.0*u0d00k+76032.0)*pow(rhoT,6.0)+(5184.0*(pow(Dinf,2.0)-48.0*pow(r_H,2.0)))*(u0d00k+1.0/9.0)*pow(rhoT,4.0)+132.0*pow(pow(Dinf,2.0)-48.0*pow(r_H,2.0),2.0)*pow(rhoT,2.0)+pow(pow(Dinf,2.0)-48.0*pow(r_H,2.0),3.0)));
			Minf = (-pow(Dinf,2.0)-32.0*pow(rhoT,2.0)+4.0*rhoT*cbrt(Mtemp) + 48.0*pow(r_H,2.0) + cbrt(pow(Mtemp,2.0)))/(12.0*cbrt(Mtemp));
			printf("Minf  = %14.7e.\n",Minf);
			ainf = -u3d00k/(-2.0*Minf/pow(rhoT,2.0)+6.0*pow(Minf,2.0)/pow(rhoT,3.0));
			printf("ainf  = %14.7e.\n",ainf);
			Jinf = ainf*Minf;
			printf("Jinf  = %14.7e.\n",Jinf);
			*/

			zeta = pow(alpha,2.0)/(beta*pow(r_H,4.0));
			dMcoeff = 1117.0/18480.0;
			dM = dMcoeff*zeta;
			dD = alpha/beta/r_H;

			printf("With an error of O(r_H/rho)^2 = %11.4e:\n",pow(r_H/rhoT,2.0));
			printf("M   = %14.7e.\n",Minf);
			printf("correction with M = %14.7e: %14.7e.\n",M, Minf/M-1.0);
			printf("D   = %14.7e.\n",Dinf);
			printf("dD  = %14.7e.\n",dD);
			printf("dM  = %14.7e.\n",dM);

		}
	}



	//-------------------- Calculate angular momentum at sufficient infinity --------------------
	//rhoT = r_H/tol;
	//rhoT = sqrt(r_H/tol);
	//rhoT = 6.0*cbrt(r_H/tol);
	rhoT = sqrt(r_H/tol) - 11.5*alpha*cbrt(r_H/tol); //Sufficient infinity depends on coupling

	xT = (1.0-r_H/rhoT);

	printf("rho such that O(1/rho) term in g_{tt} = tol: rho =  %11.4e.\n",rhoT);
	printf("xT = %14.7e.\n",xT);
	xk = xT;
	for (j = 0; j < n; ++j)
	{
		if (x[j] >= xk)
		{
			J = j;
			break;
		}
	}

	//for (i = 1; i < m-1; ++i)
	i = m-1; //Evaluate on equator
	{
		{
			k = (i*n+j)*p;
			I = i, J = j;
			yk = y[I];

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

			//Mass and Charge calc
			gtt  =	gmetric0out (xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
			gtph =	gmetric1out (xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
			gpsi =	gmetric2out (xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);

			printf("f = %11.4e. m = %11.4e. l = %11.4e. omega = %11.4e. psi = %11.4e.\n",u0d00k, u1d00k,u2d00k,u3d00k,u4d00k);
			printf("gtt = %11.4e. gtph = %11.4e. gpsi = %11.4e.\n",gtt,gtph,gpsi);

			//O(1/r)
			//Minf = rhoT/2.0*(1.0+gtt);
			//Jinf = -rhoT/2.0/pow(sin(yk),2.0)*gtph;
			//Dinf = rhoT*gpsi;

			//O(1/r^2)
			//Minf = rhoT/2.0*(1.0-sqrt(-1.0-2.0*gtt));
			Jinf = gtph/(-2.0*pow(sin(yk),2.0)/rhoT)/(1.0 - Minf/rhoT);
			ainf = Jinf/Minf;
			//Dinf = rhoT*gpsi;

			//O(1/r^3)
			/*
			C1inf = pow(rhoT,2.0)*(u2d00k-1.0);
			printf("C1inf = %14.7e.\n",C1inf);
			Dinf = u4d00k/(-1.0/rhoT+C1inf/pow(rhoT,3.0));
			printf("Dinf  = %14.7e.\n",Dinf);
			printf("check = %14.7e.\n",-pow(Dinf,2.0)/8.0-2.0*pow(r_H,2.0));
			//Minf = sqrt(pow(rhoT,2.0)*(1.0-u1d00k)-C1inf-pow(Dinf,2.0)/4.0);
			Mtemp = (-6.0*rhoT*pow(Dinf,2.0)+208.0*pow(rhoT,3.0)+288.0*pow(r_H,2.0)*rhoT-432.0*u0d00k*pow(rhoT,3.0)+sqrt((186624.0*pow(u0d00k,2.0)-179712.0*u0d00k+76032.0)*pow(rhoT,6.0)+(5184.0*(pow(Dinf,2.0)-48.0*pow(r_H,2.0)))*(u0d00k+1.0/9.0)*pow(rhoT,4.0)+132.0*pow(pow(Dinf,2.0)-48.0*pow(r_H,2.0),2.0)*pow(rhoT,2.0)+pow(pow(Dinf,2.0)-48.0*pow(r_H,2.0),3.0)));
			Minf = (-pow(Dinf,2.0)-32.0*pow(rhoT,2.0)+4.0*rhoT*cbrt(Mtemp) + 48.0*pow(r_H,2.0) + cbrt(pow(Mtemp,2.0)))/(12.0*cbrt(Mtemp));
			printf("Minf  = %14.7e.\n",Minf);
			ainf = -u3d00k/(-2.0*Minf/pow(rhoT,2.0)+6.0*pow(Minf,2.0)/pow(rhoT,3.0));
			printf("ainf  = %14.7e.\n",ainf);
			Jinf = ainf*Minf;
			printf("Jinf  = %14.7e.\n",Jinf);
			*/

			chiinf = ainf/Minf;
			dJcoeff = 0.0;
			dJ = dJcoeff*zeta;

			printf("With an error of O(r_H/rho)^3 = %11.4e:\n",pow(r_H/rhoT,3.0));
			printf("J   = %14.7e.\n",Jinf);
			printf("a   = %14.7e.\n",ainf);
			printf("Omega_H = %14.7e.\n",Omega_H);
			printf("chi = %14.7e <= 1.\n",chiinf);
			printf("M0 using ^ = %14.7e.\n",2.0*r_H/sqrt(1.0-pow(chiinf,2.0)));
			M0 = 2.0*r_H/sqrt(1.0-pow(chiinf,2.0));
			printf("J    correction; J = J0(1+dJ), dJ = %12.5e.\n",Jinf/J0 - 1.0);
			printf("compared to LIN solution dJ = M*sqrt(M^2-4*r_H^2) = %12.5e.\n", dJ);
			printf("Mass correction; M = M0(1+dM), dM = %12.5e.\n",Minf/M0 - 1.0);
			printf("compared to LIN solution dM = zeta*1117/18480 = %12.5e.\n", dM);

		}
	}



	//-------------------- ISCO & LR variables --------------------
	xISCOSCHW = -4.0+2.0*sqrt(6.0);
	rhoISCOSCHW = r_H/(1.0-xISCOSCHW);
	rISCOSCHW = rhoISCOSCHW + 2.0*r_H + pow(r_H,2.0)/rhoISCOSCHW;
	xLRSCHW = sqrt(3.0)-1.0;
	rhoLRSCHW = r_H/(1.0-xLRSCHW);
	rLRSCHW = rhoLRSCHW + 2.0*r_H + pow(r_H,2.0)/rhoLRSCHW;

	dxISCOcoeff = (213817.0/17962560.0)*(1.0-sqrt(6.0)/3.0);
	dxLRcoeff = -11833.0/4490640.0*(1.0-sqrt(3.0)/3.0);

	double FLR,dFLRdx;
	double KERRLR,dKERRLRdx;
	double KERRISCO,dKERRISCOdx;
	double FISCO,dFISCOdx;
	double FISCOdx;
	double deltax = 1.0/1.0E5;


	//LR KERR calc
	rLRKERRsol = 2.0*M*(1.0 + cos(2.0/3.0*acos(-chi)));
	rhoLRKERRsol = -M/2.0 + rLRKERRsol/2.0 + sqrt(pow(rLRKERRsol,2.0) - 2.0*M*rLRKERRsol + pow(chi*M,2.0))/2.0;
	xLRKERRsol = (1.0 - r_H/rhoLRKERRsol);
	//xk = xLRSCHW;
	xk = xLRKERRsol;
	do
	{
		KERRLR = KERRLRout(xk,r_H,M);
		dKERRLRdx = dKERRLRdxout(xk,r_H,M);

		xLR = xk - KERRLR/dKERRLRdx;
		printf("KERRLR = %12.5e. dKERRLRdx = %12.5e. xLR = %.3f. xk = %.5f.\n",KERRLR, dKERRLRdx,xLR,xk);
		xk = xLR;
	} while(fabs(KERRLR) > tol);
	xLRKERR = xLR;
	rhoLRKERR = r_H/(1.0-xLRKERR);
	rLRKERR = rhoLRKERR + M + pow(r_H,2.0)/rhoLRKERR;
	printf("rLRSCHW   = %6.2f [M]. rhoLRSCHW   = %6.2f [M]. xLRSCHW   = %6.2f.\n",rLRSCHW,rhoLRSCHW,xLRSCHW);
	printf("rLRKERRsol= %6.2f [M]. rhoLRKERRsol= %6.2f [M]. xLRKERRsol= %6.2f.\n",rLRKERRsol,rhoLRKERRsol,xLRKERRsol);
	printf("rLRKERR   = %6.2f [M]. rhoLRKERR   = %6.2f [M]. xLRKERR   = %6.2f.\n",rLRKERR,rhoLRKERR,xLRKERR);
	printf("delta xLRKERR = %12.5e\n",xLRKERR/xLRKERRsol-1.0);


	//-------------------- Light Ring --------------------
	xk = xLRSCHW*(1.0+dxLRcoeff*zeta);
	//xk = xLRKERRsol*(1.0+dxLRcoeff*zeta);
	l = 0;
	do
	{
		for (j = 0; j < n; ++j)
		{
			if (x[j] >= xk)
			{
				J = j;
				break;
			}
		}

		i = m-1; //Evaluate on equator
		{
			{
				k = (i*n+j)*p;
				I = i, J = j;
				yk = y[I];

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

				FLR = FLRout(xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
				dFLRdx = dFLRdxout(xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);

				xLR = xk - FLR/dFLRdx;
				xk = xLR;
				++l;
				printf("FLR = %e. dFLRdx = %e. xLR = %.3f. xk = %.7f.\n",FLR, dFLRdx,xLR,xk);

			}
		}

	} while(fabs(FLR) > pow(tol,2.0));
	OmegaLR = OrbLRout(xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
	printf("Omega_CLR = %14.7e.\n",OmegaLR);

	printf("Light Ring after %i iterations:\n",l);
	printf("xLRKERR   = %14.7e. rhoLRKERR   = %14.7e.\n", xLRKERRsol, rhoLRKERRsol);
	printf("xLR       = %14.7e. rhoLR       = %14.7e.\n", xLR, r_H/(1.0-xLR));
	printf("Light ring correction; x_LR = x_GR(1+dx_LR), dx_LR = %12.5e.\n",xLR/xLRKERRsol - 1.0);
	printf("compared to LIN solution dLR = %12.5e.\n", dxLRcoeff*zeta);



	//ISCO KERR calc
	double Z1, Z2;
	Z1 = 1.0 + cbrt(1.0-pow(chi,2.0))*(cbrt(1.0+chi) + cbrt(1.0-chi));
	Z2 = sqrt(3.0*pow(chi,2.0) + pow(Z1,2.0));
	printf("Z1 = %f. Z2 = %f.\n",Z1,Z2);
	rISCOKERRsol = M*(3.0 + Z2 - sqrt((3.0-Z1)*(3.0 + Z1 + 2.0*Z2)));
	rhoISCOKERRsol = -M/2.0 + rISCOKERRsol/2.0 + sqrt(pow(rISCOKERRsol,2.0) - 2.0*M*rISCOKERRsol + pow(chi*M,2.0))/2.0;
	xISCOKERRsol = (1.0 - r_H/rhoISCOKERRsol);
	/*
	//xk = xISCOSCHW;
	xk = xISCOKERRsol;
	do
	{
		KERRISCO = KERRISCOout(xk,r_H,M);
		dKERRISCOdx = KERRISCOout(xk+deltax,r_H,M);
		dKERRISCOdx = (dKERRISCOdx-KERRISCO)/deltax;
		//dKERRISCOdx = dKERRISCOdxout(xk,r_H,M);

		xISCO = xk - KERRISCO/dKERRISCOdx;
		printf("KERRISCO = %12.5e. dKERRISCOdx = %12.5e. xISCO = %.3f. xk = %.3f.\n",KERRISCO, dKERRISCOdx,xISCO,xk);
		xk = xISCO;

	} while(fabs(KERRISCO) > tol);
	xISCOKERR = xISCO;
	rhoISCOKERR = r_H/(1.0-xISCOKERR);
	rISCOKERR = rhoISCOKERR + M + pow(r_H,2.0)/rhoISCOKERR;
	*/
	printf("rISCOSCHW    = %6.2f [M]. rhoISCOSCHW    = %6.2f [M]. xISCOSCHW    = %6.2f.\n",rISCOSCHW,rhoISCOSCHW,xISCOSCHW);
	printf("rISCOKERRsol = %6.2f [M]. rhoISCOKERRsol = %6.2f [M]. xISCOKERRsol = %6.2f.\n",rISCOKERRsol,rhoISCOKERRsol,xISCOKERRsol);
	//printf("rISCOKERR    = %6.2f [M]. rhoISCOKERR    = %6.2f [M]. xISCOKERR    = %6.2f.\n",rISCOKERR,rhoISCOKERR,xISCOKERR);
	//printf("delta xISCOKERR = %12.5e\n",xISCOKERR/xISCOKERRsol-1.0);



	//-------------------- ISCO --------------------
	//xk = xISCOSCHW;
	xk = xISCOKERRsol*(1.0+dxISCOcoeff*zeta);
	l = 0;
	do
	{
		for (j = 0; j < n; ++j)
		{
			if (x[j] >= xk)
			{
				J = j;
				break;
			}
		}

		i = m-1; //Evaluate on equator
		{
			{
				k = (i*n+j)*p;
				I = i, J = j;
				yk = y[I];

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

				FISCO = FISCOout(xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
				FISCOdx = FISCOout(xk+deltax,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
				dFISCOdx = (FISCOdx-FISCO)/deltax;

				xISCO = xk - FISCO/dFISCOdx;
				printf("FISCO = %e. dFISCOdx = %e. xISCO = %.3f. xk = %.7f.\n",FISCO, dFISCOdx,xISCO,xk);
				xk = xISCO;
				++l;

			}
		}

	} while(fabs(FISCO) > pow(tol,2.0));
	OmegaISCO = OrbISCOout(xk,yk,r_H,M,alpha,beta, u0d00k, u0d10k, u0d20k, u1d00k, u1d10k, u1d20k, u2d00k, u2d10k, u2d20k, u3d00k, u3d10k, u3d20k, u4d00k, u4d10k, u4d20k);
	printf("Omega_CISCO = %14.7e.\n",OmegaISCO);

	printf("ISCO after %i iterations:\n",l);
	printf("xISCOGR = %14.7e. rhoISCOGR = %14.7e.\n", xISCOKERRsol, rhoISCOKERRsol);
	printf("xISCO   = %14.7e. rhoSICO   = %14.7e.\n",xISCO,r_H/(1.0-xISCO));
	printf("ISCO correction; x_ISCO = x_GR(1+dx_ISCO), dx_ISCO = %12.5e.\n",xISCO/xISCOKERRsol - 1.0);
	printf("compared to LIN solution dISCO = %12.5e.\n", dxISCOcoeff*zeta);



	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	//-------------------- Observables Output --------------------

	char foutname[256];
	FILE *foutfile;

	sprintf(foutname,"./Data/PlotData/PropertiesData/");
	strcat(foutname, SolType);
	char buffer[256];
	sprintf(buffer,"prop_chi%0i_tol%03i_rH%03i.dat",(int) round(chi*1.0E1), (int) round(-log10(tol)), (int) round(r_H*1.0E2));
	strcat(foutname, buffer);
	printf("Observables out file name: %s\n",foutname);
	foutfile = fopen(foutname,"a");
	if(foutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file

	//alpha, r_H, chi, Omega_H, M0, J0
	//M, dM, dMpert,
	//J, dJ, dJpert,
	//D, dDpert,
	//xISCO, dxISCO, dxISCOpert,
	//xLR, dxLR, dxLRpert;
	//OmegaLR, OmegaISCO;

	fprintf(foutfile,"%.3f %.2f %14.7e %14.7e %14.7e %14.7e ", alpha/pow(r_H,2.0), r_H, chiinf, Omega_H, M0, J0);
	fprintf(foutfile,"%14.7e %14.7e %14.7e ",	Minf, Minf/M0-1.0, dMcoeff*zeta);
	fprintf(foutfile,"%14.7e %14.7e %14.7e ",	Jinf, Jinf/J0-1.0, dJcoeff*zeta);
	fprintf(foutfile,"%14.7e %14.7e ", Dinf, dD);
	fprintf(foutfile,"%14.7e %14.7e %14.7e ", xISCO, xISCO/xISCOKERRsol-1.0, dxISCOcoeff*zeta);
	fprintf(foutfile,"%14.7e %14.7e %14.7e ", xLR, xLR/xLRKERRsol-1.0, dxLRcoeff*zeta);
	fprintf(foutfile,"%14.7e %14.7e\n", OmegaLR, OmegaISCO);

	fclose(foutfile);


	//Make directory for comparison angle plots
/*
	int check;
	char dirname[256];
	sprintf(dirname,"./Data/PlotData/angleplotdata_al%06i_chi%03i", (int) round(alpha*1.0E4), (int) round(chi*1.0E2));

	DIR* dir = opendir(dirname);
	if (dir) {
	    //Directory exists
	    closedir(dir);
			printf("Directory %s exists.\n",dirname);

	} else if (ENOENT == errno) {
	    //Directory does not exist.
			check = mkdir(dirname, 0755);
			if (!check)
				 printf("Directory created\n");
			else {
				 printf("Unable to create directory\n");
				 exit(1);
			}
	} else {
	    //opendir() failed for some other reason.
	}


	//-------------------- EDGB from kerr metric plot o6utput --------------------
	int Iyplot = 3;
	int yplot[3] = {0,(m-1)/2,(m-1)};

	char pkoutname[256];
	FILE *pkoutfile;

	strcpy(pkoutname, dirname);
	strcat(pkoutname, "/");
	strcat(pkoutname, SolType);
	strcat(pkoutname, "fromKERR.dat");
	printf("from GR out file name: %s\n",pkoutname);

	pkoutfile = fopen(pkoutname,"w");
	if(pkoutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (j = 0; j < n; ++j)
	{
		xk = x[j];
		fprintf(pkoutfile,"%22.15e",xk);

		for (i = 0; i < Iyplot; ++i)
		{
			yk = y[yplot[i]];
			k = (yplot[i]*n+j)*p;
			fprintf(pkoutfile," %22.15e %22.15e %22.15e %22.15e %22.15e", Psi[k +0] - KERRu0Xout (xk,yk,r_H,M),
			Psi[k +1] - KERRu1Xout (xk,yk,r_H,M),
			Psi[k +2] - KERRu2Xout (xk,yk,r_H,M),
			Psi[k +3] - KERRu3Xout (xk,yk,r_H,M),
			Psi[k +4] - KERRu4Xout (xk,yk,r_H,M));
		}
		fprintf(pkoutfile,"\n");
	}
	fclose(pkoutfile);


	//-------------------- EDGB diff metric plot output --------------------

	char psoutname[256];
	FILE *psoutfile;

	strcpy(psoutname, dirname);
	strcat(psoutname, "/");
	strcat(psoutname, SolType);
	strcat(psoutname, "fromSPHsGB.dat");
	printf("from SPH Pert out file name: %s\n",psoutname);

	psoutfile = fopen(psoutname,"w");
	if(psoutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (j = 0; j < n; ++j)
	{
		xk = x[j];
		fprintf(psoutfile,"%22.15e",xk);

		for (i = 0; i < Iyplot; ++i)
		{
			yk = y[yplot[i]];
			k = (yplot[i]*n+j)*p;
			fprintf(psoutfile," %22.15e %22.15e %22.15e %22.15e %22.15e", Psi[k +0] - pow(xk,2.0)*IC0out(xk,yk,r_H,M,alpha,1.0,0.0),
			Psi[k +1] - pow(xk,2.0)*IC1out(xk,yk,r_H,M,alpha,1.0,0.0),
			Psi[k +2] - pow(xk,2.0)*IC2out(xk,yk,r_H,M,alpha,1.0,0.0),
			Psi[k +3] - IC3out(xk,yk,r_H,M,alpha,1.0,0.0),
			Psi[k +4] - IC4out(xk,yk,r_H,M,alpha,1.0,0.0));
		}
		fprintf(psoutfile,"\n");
	}
	fclose(psoutfile);
*/

	//-------------------- Kerr metric plot o6utput --------------------
/*
	int yplot[3] = {0,(m-1)/2,(m-1)};
	char poutname[256];
	FILE *poutfile;

	sprintf(poutname,"./Data/PlotData/KERRplotdata/KERRerror_angleplotdata_chi%03i_n%03i_m%02i.dat", (int) round(chi*1.0E2), n, m);
	poutfile = fopen(poutname,"w");
	if(poutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (j = 0; j < n; ++j)
	{
		xk = x[j];
		fprintf(poutfile,"%22.15e",xk);
// - pow(xk,2.0)*IC0out(xk,yk,r_H,2.0*r_H,0.0,1.0,0.0)
		for (i = 0; i < 3; ++i)
		{
			yk = y[yplot[i]];
			k = (yplot[i]*n+j)*p;
			fprintf(poutfile," %22.15e %22.15e %22.15e %22.15e %22.15e", Psi[k +0] - KERRu0Xout (xk,yk,r_H,M),
			Psi[k +1] - KERRu1Xout (xk,yk,r_H,M),
			Psi[k +2] - KERRu2Xout (xk,yk,r_H,M),
			Psi[k +3] - KERRu3Xout (xk,yk,r_H,M),
			Psi[k +4] - KERRu4Xout (xk,yk,r_H,M));
		}
		fprintf(poutfile,"\n");
	}
	fclose(poutfile);
*/


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


	//Cleanup
	fclose(infile);
  free(Psi);
  free(xGrid);
  free(yGrid);
	free(x);
	free(y);


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
