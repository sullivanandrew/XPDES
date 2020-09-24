#include "BVP_header.h"

int BVP_out(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], double Dx[], double Dy[], int it)
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
	const double r_H = (*params).r_Hparam;
	const double M = (*params).Mparam;
	const double chi = (*params).chiparam;
	const double ICalpha = (*params).ICalphaparam;
	const double ICchi = (*params).ICchiparam;

	int pn;
	double xk, yk;


	char poutname[256];
	FILE *poutfile;

	if (it == -3) //Save solution in separate directory
	{
		sprintf(poutname,"./Data/BVPout/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
		//sprintf(poutname,"./Data/BVPoutLIN/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
		//sprintf(poutname,"./Data/BVPoutEXP/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
		//sprintf(poutname,"./Data/BVPoutKERR/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i_ICdel%03i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m, (int) round(ICchi*1.0E3));
		//sprintf(poutname,"./Data/BVPoutPERT/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	}
	else
	{
		sprintf(poutname,"./Data/BVPout/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
		//sprintf(poutname,"./Data/BVPoutLIN/BVPout_iterations/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
		//sprintf(poutname,"./Data/BVPoutEXP/BVPout_sols/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
		//sprintf(poutname,"./Data/BVPoutKERR/BVPout_iterations/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i_ICdel%03i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m, (int) round(ICchi*1.0E3));
		//sprintf(poutname,"./Data/BVPoutPERT/BVPout_iterations/sol_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	}
	poutfile = fopen(poutname,"w");
	if(poutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (i = 0; i < m; ++i)
	{
		yk = y[i];
		for (j = 0; j < n; ++j)
		{
			xk = x[j];
			fprintf(poutfile,"%22.15e %22.15e",xk,yk);
			for (pn = 0; pn < p; ++pn)
			{
				k = (i*n+j)*p;
				fprintf(poutfile," %22.15e %22.15e %22.15e %22.15e",Psi[k +pn], b[k +pn], Dx[k +pn], Dy[k +pn]);
			}
			fprintf(poutfile,"\n");

		}
	}
	fclose(poutfile);


	return 0;
}


int BVP_sysout(struct param_type *params, double JacVal[], int JacRow[], int JacCol[], double dPsi[], double b[], double Ddx[], double Dx[], double Ddy[], double Dy[], int it)
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
	const double r_H = (*params).r_Hparam;
	const double M = (*params).Mparam;
	const double chi = (*params).chiparam;
	const double ICalpha = (*params).ICalphaparam;
	const double ICchi = (*params).ICchiparam;
	int pn;
	double xk, yk;


	char poutname[256];
	FILE *poutfile;
	char boutname[256];
	FILE *boutfile;

	sprintf(poutname,"./Data/BVPout/BVPout_sys/sysJAC_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	//sprintf(poutname,"./Data/BVPoutLIN/BVPout_sys/sysJAC_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	//sprintf(poutname,"./Data/BVPoutEXP/BVPout_sys/sysJAC_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	//sprintf(poutname,"./Data/BVPoutKERR/BVPout_sys/sysJAC_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i_ICdel%03i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m, (int) round(ICchi*1.0E3));
	//sprintf(poutname,"./Data/BVPoutPERT/BVPout_sys/sysJAC_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	poutfile = fopen(poutname,"w");
	if(poutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (i = 0; i < nnzJac; ++i)
	{
		fprintf(poutfile,"%22.15e %i %i\n", JacVal[i], JacRow[i], JacCol[i]);
	}
	fclose(poutfile);


	sprintf(boutname,"./Data/BVPout/BVPout_sys/sysB_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	//sprintf(boutname,"./Data/BVPoutLIN/BVPout_sys/sysB_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	//sprintf(boutname,"./Data/BVPoutEXP/BVPout_sys/sysB_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	//sprintf(boutname,"./Data/BVPoutKERR/BVPout_sys/sysB_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%03i_m%02i_ICdel%03i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m, (int) round(ICchi*1.0E3));
	//sprintf(boutname,"./Data/BVPoutPERT/BVPout_sys/sysB_al%06i_tol%03i_it%02i_rH%03i_chi%03i_n%02i_m%02i.dat",(int) round(alpha*1.0E4), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2), (int) round(chi*1.0E2), n, m);
	boutfile = fopen(boutname,"w");
	if(boutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (i = 0; i < N; ++i)
	{
			fprintf(boutfile,"%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", dPsi[i], b[i], Ddx[i], Dx[i], Ddy[i], Dy[i]);
	}
	fclose(boutfile);


	return 0;
}
