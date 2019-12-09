#include "BVP_header.h"

//Temporary placeholder

int BVP_physics(void *params, double x[], double Psi[], double bvec[], double Dvec[], double RGB[], double R2[], double Rab2[], double Rabcd2[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int p = my_ptr->pparam;
	const int N = my_ptr->Nparam;
	const int r = my_ptr->rparam;
	const double r_H = my_ptr->r_Hparam;
	const double tol = my_ptr->tolparam;
	const double alpha = my_ptr->alphaparam;
	const double beta = my_ptr->betaparam;
	const double gamma = my_ptr->gammaparam;
	const double kappa = my_ptr->kappaparam;
	int i,j,k,l,status;
	int pn;
	double xk;
	int a;
	const double sgn = -1.0;
	int I, oneside, twoside;

	double *a1_a;
	double *a1x_a;
	double *a1xx_a;
	a1_a = calloc(r+4, sizeof(double));
	a1x_a = calloc(r+4, sizeof(double));
	a1xx_a = calloc(r+4, sizeof(double));
	if(!a1_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a1x_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a1xx_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	double *a2_a;
	double *a2x_a;
	double *a2xx_a;
	a2_a = calloc(r+4, sizeof(double));
	a2x_a = calloc(r+4, sizeof(double));
	a2xx_a = calloc(r+4, sizeof(double));
	if(!a2_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a2x_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a2xx_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	double *a3_a;
	double *a3x_a;
	double *a3xx_a;
	a3_a = calloc(r+4, sizeof(double));
	a3x_a = calloc(r+4, sizeof(double));
	a3xx_a = calloc(r+4, sizeof(double));
	if(!a3_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a3x_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	if(!a3xx_a){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	int *s;
	s = (int*)calloc(r+3, sizeof(int));
	if(!s){
        printf("Error! Memory not allocated. Exiting");
		exit(0);
    }
	double ps1k,dps1k,d2ps1k;
	double ps1dk,dps1dk,d2ps1dk;
	double ps2k,dps2k,d2ps2k;
	double ps2dk,dps2dk,d2ps2dk;
	double ps3k,dps3k,d2ps3k;
	double ps3dk,dps3dk,d2ps3dk;

	double Tcoupl, Kcoupl;
	double Gpsicoupl, GBcoupl;

	double Gtt, Ttt, Ktt;
	double dGttdf2,dGttdf1,dGttdf0;
	double dGttdm2,dGttdm1,dGttdm0;
	double dGttdp2,dGttdp1,dGttdp0;
	double dTttdf2,dTttdf1,dTttdf0;
	double dTttdm2,dTttdm1,dTttdm0;
	double dTttdp2,dTttdp1,dTttdp0;
	double dKttdf2,dKttdf1,dKttdf0;
	double dKttdm2,dKttdm1,dKttdm0;
	double dKttdp2,dKttdp1,dKttdp0;

	double Grr, Trr, Krr;
	double dGrrdf2,dGrrdf1,dGrrdf0;
	double dGrrdm2,dGrrdm1,dGrrdm0;
	double dGrrdp2,dGrrdp1,dGrrdp0;
	double dTrrdf2,dTrrdf1,dTrrdf0;
	double dTrrdm2,dTrrdm1,dTrrdm0;
	double dTrrdp2,dTrrdp1,dTrrdp0;
	double dKrrdf2,dKrrdf1,dKrrdf0;
	double dKrrdm2,dKrrdm1,dKrrdm0;
	double dKrrdp2,dKrrdp1,dKrrdp0;

	double Gpsi, GB, GBcalc, R2calc, Rab2calc, Rabcd2calc;
	double dGpsidf2, dGpsidf1, dGpsidf0;
	double dGpsidm2, dGpsidm1, dGpsidm0;
	double dGpsidp2, dGpsidp1, dGpsidp0;
	double dGBdf2, dGBdf1, dGBdf0;
	double dGBdm2, dGBdm1, dGBdm0;
	double dGBdp2, dGBdp1, dGBdp0;

	double FISCO, FLR;
	double dFISCOdx, dFLRdx;

	//THESE are correct field equations. \beta = \kappa for kinetic term
	Tcoupl = 1.0/2.0*beta/kappa;
	Kcoupl = alpha/kappa;
	Gpsicoupl = beta/kappa;
	GBcoupl = alpha/kappa;


	//New variables
	double rhoT;
	double xT;
	double M0;
	double M;
	double Mf;
	double dMcoeff;
	double D, dD;
	double zeta;
	double xISCO, xISCOGR;
	double xLR, xLRGR;
	double rhoISCO, rhoISCOGR;
	double rhoLR, rhoLRGR;
	double dxISCOcoeff, dxLRcoeff;


	//-------------------- Calculate mass at and charge at sufficient infinity --------------------
	rhoT = r_H/tol;
	xT = (1.0-r_H/4.0/rhoT)/(1.0+r_H/4.0/rhoT);

	printf("rho such that O(1/rho) term in g_{tt} = tol: rho =  %11.4e.\n",rhoT);
	printf("xT = %23.16e.\n",xT);

	j = n-1; //Final point to be over written
	k = j*p;
	xk = xT;
	I = 0;
	for (i = 0; i < n; ++i)
	{
		if (x[i] >= xk)
		{
			I = i;
			break;
		}
	}

	ps1k = 0.0, dps1k = 0.0, d2ps1k = 0.0;
	ps1dk = 0.0, dps1dk = 0.0, d2ps1dk = 0.0;
	ps2k = 0.0, dps2k = 0.0, d2ps2k = 0.0;
	ps2dk = 0.0, dps2dk = 0.0, d2ps2dk = 0.0;
	ps3k = 0.0, dps3k = 0.0, d2ps3k = 0.0;
	ps3dk = 0.0, dps3dk = 0.0, d2ps3dk = 0.0;
	oneside = 0;
	twoside = 0;

	//Compute stencil
	status = steninit(s, &oneside, &twoside, I, r, n);
	//Compute Newton polynomial representation for each function
	status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 0, a1_a, a1x_a, a1xx_a, &ps1dk, &dps1dk, &d2ps1dk, xk);
	status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 1, a2_a, a2x_a, a2xx_a, &ps2dk, &dps2dk, &d2ps2dk, xk);
	status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 2, a3_a, a3x_a, a3xx_a, &ps3dk, &dps3dk, &d2ps3dk, xk);

	for (a = 0; a <= r+oneside; ++a)
	{
		ps1k   += Psi[(I+s[a])*p +0]*a1_a[a];
		dps1k  += Psi[(I+s[a])*p +0]*a1x_a[a];
		d2ps1k += Psi[(I+s[a])*p +0]*a1xx_a[a];
		ps2k   += Psi[(I+s[a])*p +1]*a2_a[a];
		dps2k  += Psi[(I+s[a])*p +1]*a2x_a[a];
		d2ps2k += Psi[(I+s[a])*p +1]*a2xx_a[a];
		ps3k   += Psi[(I+s[a])*p +2]*a3_a[a];
		dps3k  += Psi[(I+s[a])*p +2]*a3x_a[a];
		d2ps3k += Psi[(I+s[a])*p +2]*a3xx_a[a];
	}

	//Mass and Charge calc
	M0 = r_H/2.0;
	zeta = pow(alpha,2.0)/(beta*kappa*pow(r_H,4.0));
	dMcoeff = 49.0/5.0;
	Mf = rhoT/2.0*(1.0 - ps1k);
	D = rhoT*ps3k;
	dD = 4.0*alpha/beta/r_H;

	printf("Mf = %23.16e, MO2 = %23.16e, with an error of O(1/rho) = %11.4e.\n",Mf,M0*(1.0+dMcoeff*zeta),1.0/rhoT);
	printf("Mass correction; M = M0(1+dM), dM = %12.5e. compared to LIN solution dM = zeta*49/5 = %12.5e.\n",Mf/M0 - 1.0, dMcoeff*zeta);
	printf("Charge correction; D = %12.5e. compared to LIN solution dD = %12.5e.\n",D, dD);


	//-------------------- ISCO --------------------
	xISCOGR = sqrt(6.0)/3.0;
	xLRGR = sqrt(3.0)/3.0;
	rhoISCOGR = r_H/4.0*(1.0+xISCOGR)/(1.0-xISCOGR);
	rhoLRGR = r_H/4.0*(1.0+xISCOGR)/(1.0-xISCOGR);
	dxISCOcoeff = 427634.0/841995.0 + 2383.0/13860.0*sqrt(6.0);
	dxLRcoeff = -189328.0/841995.0 + 2383.0/3465.0*sqrt(3.0);

	xk = 0.6;
	do
	{
		I = 0;
		for (i = 0; i < n; ++i)
		{
			if (x[i] >= xk)
			{
				I = i;
				break;
			}
		}

		ps1k = 0.0, dps1k = 0.0, d2ps1k = 0.0;
		ps1dk = 0.0, dps1dk = 0.0, d2ps1dk = 0.0;
		ps2k = 0.0, dps2k = 0.0, d2ps2k = 0.0;
		ps2dk = 0.0, dps2dk = 0.0, d2ps2dk = 0.0;
		ps3k = 0.0, dps3k = 0.0, d2ps3k = 0.0;
		ps3dk = 0.0, dps3dk = 0.0, d2ps3dk = 0.0;
		oneside = 0;
		twoside = 0;

		//Compute stencil
		status = steninit(s, &oneside, &twoside, I, r, n);
		//Compute Newton polynomial representation for each function
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 0, a1_a, a1x_a, a1xx_a, &ps1dk, &dps1dk, &d2ps1dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 1, a2_a, a2x_a, a2xx_a, &ps2dk, &dps2dk, &d2ps2dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 2, a3_a, a3x_a, a3xx_a, &ps3dk, &dps3dk, &d2ps3dk, xk);

		for (a = 0; a <= r+oneside; ++a)
		{
			ps1k   += Psi[(I+s[a])*p +0]*a1_a[a];
			dps1k  += Psi[(I+s[a])*p +0]*a1x_a[a];
			d2ps1k += Psi[(I+s[a])*p +0]*a1xx_a[a];
			ps2k   += Psi[(I+s[a])*p +1]*a2_a[a];
			dps2k  += Psi[(I+s[a])*p +1]*a2x_a[a];
			d2ps2k += Psi[(I+s[a])*p +1]*a2xx_a[a];
			ps3k   += Psi[(I+s[a])*p +2]*a3_a[a];
			dps3k  += Psi[(I+s[a])*p +2]*a3x_a[a];
			d2ps3k += Psi[(I+s[a])*p +2]*a3xx_a[a];
		}

		FISCO = FISCOout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dFISCOdx = dFISCOdxout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		xISCO = xk - FISCO/dFISCOdx;
		xk = xISCO;

	} while(fabs(FISCO) > tol);

	printf("xISCOGR = %23.16e. rhoISCOGR = %23.16e.\n", xISCOGR, rhoISCOGR);
	printf("xISCO   = %23.16e. rhoSICO   = %23.16e.\n",xISCO,r_H/4.0*(1.0+xISCO)/(1.0-xISCOGR));
	printf("ISCO correction; x_ISCO = x_GR(1+dx_ISCO), dx_ISCO = %12.5e. compared to LIN solution dISCO = %12.5e.\n",xISCO/xISCOGR - 1.0, dxISCOcoeff*zeta);


	//-------------------- Light Ring --------------------
	xk = 0.6;
	do
	{
		I = 0;
		for (i = 0; i < n; ++i)
		{
			if (x[i] >= xk)
			{
				I = i;
				break;
			}
		}

		ps1k = 0.0, dps1k = 0.0, d2ps1k = 0.0;
		ps1dk = 0.0, dps1dk = 0.0, d2ps1dk = 0.0;
		ps2k = 0.0, dps2k = 0.0, d2ps2k = 0.0;
		ps2dk = 0.0, dps2dk = 0.0, d2ps2dk = 0.0;
		ps3k = 0.0, dps3k = 0.0, d2ps3k = 0.0;
		ps3dk = 0.0, dps3dk = 0.0, d2ps3dk = 0.0;
		oneside = 0;
		twoside = 0;

		//Compute stencil
		status = steninit(s, &oneside, &twoside, I, r, n);
		//Compute Newton polynomial representation for each function
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 0, a1_a, a1x_a, a1xx_a, &ps1dk, &dps1dk, &d2ps1dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 1, a2_a, a2x_a, a2xx_a, &ps2dk, &dps2dk, &d2ps2dk, xk);
		status = BVP_NPcalc(params, r+oneside, r+oneside+twoside+2, I, s, x, Psi, 2, a3_a, a3x_a, a3xx_a, &ps3dk, &dps3dk, &d2ps3dk, xk);

		for (a = 0; a <= r+oneside; ++a)
		{
			ps1k   += Psi[(I+s[a])*p +0]*a1_a[a];
			dps1k  += Psi[(I+s[a])*p +0]*a1x_a[a];
			d2ps1k += Psi[(I+s[a])*p +0]*a1xx_a[a];
			ps2k   += Psi[(I+s[a])*p +1]*a2_a[a];
			dps2k  += Psi[(I+s[a])*p +1]*a2x_a[a];
			d2ps2k += Psi[(I+s[a])*p +1]*a2xx_a[a];
			ps3k   += Psi[(I+s[a])*p +2]*a3_a[a];
			dps3k  += Psi[(I+s[a])*p +2]*a3x_a[a];
			d2ps3k += Psi[(I+s[a])*p +2]*a3xx_a[a];
		}

		FLR = FLRout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);
		dFLRdx = dFLRdxout(ps1k,ps2k,ps3k,r_H,xk,dps1k,dps2k,dps3k,d2ps1k,d2ps2k,d2ps3k);

		xLR = xk - FLR/dFLRdx;
		xk = xLR;

	} while(fabs(FLR) > tol);

	printf("xLRGR   = %23.16e. rhoLRGR   = %23.16e.\n", xLRGR, rhoLRGR);
	printf("xLR     = %23.16e. rhoLR     = %23.16e.\n", xLR, r_H/4.0*(1.0+xLR)/(1.0-xLR));
	printf("Light ring correction; x_LR = x_GR(1+dx_LR), dx_LR = %12.5e. compared to LIN solution dLR = %12.5e.\n",xLR/xLRGR - 1.0, dxLRcoeff*zeta);

	//-------------------- Output --------------------
	char foutname[256];
	FILE *foutfile;
	sprintf(foutname,"../../Data/BVPout_sols/Sol_Props/prop_tol%03i_rH%03i.dat",(int) round(-log10(tol)), (int) round(r_H*1.0E2));
	foutfile = fopen(foutname,"a");
	if(foutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file

	fprintf(foutfile,"%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n", alpha, r_H, rhoT, Mf, Mf/M0 -1.0, dMcoeff*zeta, D, dD, xISCO, xISCO/xISCOGR -1.0, dxISCOcoeff*zeta, xLR, xLR/xLRGR -1.0, dxLRcoeff*zeta);
	//alpha, r_H, rhoT, M, dM, dMpert, D, dDpert, xISCO, dxISCO, dxISCOpert, xLR, dxLR, dxLRpert;

	fclose(foutfile);


	//Cleanup
	free(a1_a);
	free(a1x_a);
	free(a1xx_a);
	free(a2_a);
	free(a2x_a);
	free(a2xx_a);
	free(a3_a);
	free(a3x_a);
	free(a3xx_a);
	free(s);

	return 0;
}
