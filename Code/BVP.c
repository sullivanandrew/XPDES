/*
 *Andrew Sullivan
 *Montana State University
 *Newton-Rapson method for solving PDE BVP in 2D as coupled system
 *Procedure:
 *Notes:
 *precision - double
 */

#include "BVP_header.h"

int main(int argc, char *argv[])
{
	int i,j,k,l, status;
	//External Input
	const int p = atoi(argv[1]); //Number of field equations. Will always be twice the number of fields due to mixed derivatives
	int n = atoi(argv[2]); //x grid size
	int m = atoi(argv[3]); //theta grid size
	const int r = atoi(argv[4]); //Newton polynomial order
	const double tol = pow(10.0,atof(argv[5])); //absolute tolerance of solver
	const int plotvar = atoi(argv[6]); //1 to output each iteration. 0 for only solution
	const int linsolvar = atoi(argv[7]); //Linear solver type. See BVP_LSsolver.c
	const double alpha = atof(argv[8]); //coupling parameter alpha
	const double beta = 1.0; //coupling beta set to 1
	const double ICdiff = atof(argv[9]); //delta for initial condition
	const double r_H = atof(argv[10]); //rH value. Note: this is technically \bar{r}_H = r_H/4 <- Schwarzschild
	const double chi = atof(argv[11]); //dimensionless spin parameter
	const double M = 2.0*r_H/sqrt(1.0-pow(chi,2.0)); //bare mass calculated from spin and rH
	const double a = chi*M; //spin
	int N = n*m*p; //total grid size

	//Print statements
	printf("\nRunning BVP with tol = %8.1e,\n",tol);
	printf("Coupling value alpha = %.5f, beta = %.2f.\n",alpha,beta);
	printf("Black hole parameters are r_bar_H = %5.2f [Msun]. M = %5.2f [Msun].\n",r_H,M);
	printf("This corresponds to a GR Schw BH with horizon r_H = %5.2f [Msun] = 2*M.\n",4.0*r_H);
	printf("Chi = a/M = %5.2f []. a = %5.2f [Msun].\n",chi,a);
	printf("Derivative order r = %2i.\n",r);
	printf("Grid size nxm = %3i x %3i. Total pts = %i.\n",n,m,N);

	//Executable
	status = BVP_main(plotvar, linsolvar, r, n, m, p, N, tol, chi, ICdiff, r_H, M, alpha, beta);
	if (status != 0)
	{
		printf("ERROR in BVP. status = %i.\n",status);
	}
	printf("End of BVP.\n");

	return 0;
}
