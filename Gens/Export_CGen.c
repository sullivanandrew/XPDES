#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int CodeGen_JacInit(FILE *outfile, const int FEn, const int maxderiv);
int CodeGen_NewtCalc(FILE *outfile, int FEn, int maxderiv);
int CodeGen_JacCleanup(FILE *outfile, int FEn, int maxderiv);
int CodeGen_BCDomain(FILE *outfile, int bc);
#include "../Funcs/CodeGen_GENNNZ.c"

int main(int argc, char *argv[])
{
	int i,j,k,l,bc, status;
	char joutname[4096];
	FILE *joutfile;
	char iboutname[4096];
	FILE *iboutfile;
	char var[4096];
	char varBC[4096];
	char vartemp[4096];

	const int FEn = atoi(argv[1]);
	const int maxderiv = atoi(argv[2]);
	const char BCdim[5] = "XXYY";
	const int BCloc[4] = {0,1,1,0};
	const int BClength = 4;
	size_t FEdlength;

	int *FEdx;
	int *FEdy;

	if (maxderiv == 3)
	{
		FEdlength = 7;
		FEdx = (int*)calloc(FEdlength, sizeof(int));
		FEdy = (int*)calloc(FEdlength, sizeof(int));
		FEdx[0] = 0;
		FEdx[1] = 1;
		FEdx[2] = 0;
		FEdx[3] = 2;
		FEdx[4] = 0;
		FEdx[5] = 3;
		FEdx[6] = 0;
		FEdy[0] = 0;
		FEdy[1] = 0;
		FEdy[2] = 1;
		FEdy[3] = 0;
		FEdy[4] = 2;
		FEdy[5] = 0;
		FEdy[6] = 3;
	}
	else if (maxderiv == 2)
	{
		FEdlength = 5;
		FEdx = (int*)calloc(FEdlength, sizeof(int));
		FEdy = (int*)calloc(FEdlength, sizeof(int));
		FEdx[0] = 0;
		FEdx[1] = 1;
		FEdx[2] = 0;
		FEdx[3] = 2;
		FEdx[4] = 0;
		FEdy[0] = 0;
		FEdy[1] = 0;
		FEdy[2] = 1;
		FEdy[3] = 0;
		FEdy[4] = 2;
	}
	else if (maxderiv == 1)
	{
		FEdlength = 3;
		FEdx = (int*)calloc(FEdlength, sizeof(int));
		FEdy = (int*)calloc(FEdlength, sizeof(int));
		FEdx[0] = 0;
		FEdx[1] = 1;
		FEdx[2] = 0;
		FEdy[0] = 0;
		FEdy[1] = 0;
		FEdy[2] = 1;
	}
	else
	{
		printf("ERROR!!! Invalid Maximum Derivative value = %i.\n",maxderiv);
		return -1;
	}
	const int BCdx[3] = {0,1,0};
	const int BCdy[3] = {0,0,1};
	size_t BCdlength = 3;

	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- SYS Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sprintf(joutname,"./Code/BVP_GENsys.c");
	joutfile = fopen(joutname,"w");
	if(joutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file


	fprintf(joutfile,"#include \"BVP_header.h\"\n\n");
	fprintf(joutfile,"int BVP_GENsys(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[])\n");
	fprintf(joutfile,"{\n");


	//-------------------- Jac Init --------------------
	status = CodeGen_JacInit(joutfile, FEn, maxderiv);
	//--------------------------------------------------


	//==================== Jacobian Partials and Structure Initialization ====================
	fprintf(joutfile,"	//----- Jacobian partials -----\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(joutfile,"	double FE%i;\n",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(joutfile,"	double"); //fix. Fix what???
			for (k = 0; k < FEdlength; ++k)
			{
				fprintf(joutfile," dFE%id%id%i%i",i,j,FEdx[k],FEdy[k]); //fix. fix whatt????
				if ((k+1) != FEdlength)
				{
					fprintf(joutfile,",");
				}
			}
			fprintf(joutfile,";\n");
		}
		fprintf(joutfile,"	\n");
	}

	fprintf(joutfile,"	//----- Rezero system -----\n");
	fprintf(joutfile,"	nnztemp = 0;\n");
	fprintf(joutfile,"	for (i = 0; i < N; ++i)\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		bvec[i]  = 0.0;\n");
	fprintf(joutfile,"		Dxvec[i] = 0.0;\n");
	fprintf(joutfile,"		Dyvec[i] = 0.0;\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	for (k = 0; k < nnzJac; ++k)\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		JacVal[k] = 0.0;\n");
	fprintf(joutfile,"		JacRow[k] = 0;\n");
	fprintf(joutfile,"		JacCol[k] = 0;\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	nnztemp = BVP_GENBC(params, x, y, Psi, JacVal, JacRow, JacCol, bvec, Dxvec, Dyvec);\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	\n");

	strcpy(var,"(xk, yk, r_H, M, alpha, beta,");
	for (j = 0 ; j < FEn; ++j)
	{
		for (k = 0; k < FEdlength; ++k)
		{
			sprintf(vartemp," u%id%i%ik",j,FEdx[k],FEdy[k]);
			if ((k+1) != FEdlength)
			{
				strcat(vartemp,",");
			}
			strcat(var, vartemp);
		}
		if ((j+1) != FEn)
		{
			strcat(var,",");
		}
	}
	strcat(var,")");


	//:::::::::::::::::::: Main Loop ::::::::::::::::::::
	fprintf(joutfile,"	//:::::::::::::::::::: Main Loop ::::::::::::::::::::\n");
	fprintf(joutfile,"	for (i = 1; i < m-1; ++i)\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		for (j = 1; j < n-1; ++j)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			k = (i*n+j)*p;\n");
	fprintf(joutfile,"			I = i, J = j;\n");
	fprintf(joutfile,"			xk = x[j];\n");
	fprintf(joutfile,"			yk = y[i];\n");
	fprintf(joutfile,"			\n");


	//-------------------- Newt Calc --------------------
	status = CodeGen_NewtCalc(joutfile, FEn, maxderiv);
	//---------------------------------------------------


	fprintf(joutfile,"			for (a = 0; a <= r+onexside; ++a)\n");
	fprintf(joutfile,"			{\n");

	for (j = 0; j < FEn; ++j)
	{
		fprintf(joutfile,"				u%id%i0k += Psi[(I*n +J+sx[a])*p +%i]*a%id%i[a];\n",j,0,j,j,0);
	}
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			\n");

	//==================== Jacobian Partials and Structure Calculation ====================
	fprintf(joutfile,"			//---------- Evaluate partial derivatives and structure for Jacobian ----------\n");
	fprintf(joutfile,"			\n");

	for (i = 0; i < FEn; ++i)
	{
		fprintf(joutfile,"			//FE%i\n",i);
		fprintf(joutfile,"			FE%i       =       FE%iout%s;\n",i,i,var);
		fprintf(joutfile,"			\n");

		for (j = 0; j < FEn; ++j)
		{
			for (k = 0; k < FEdlength; ++k)
			{
				fprintf(joutfile,"			dFE%id%id%i%i = dFE%id%id%i%iout%s;\n",i,j,FEdx[k],FEdy[k],i,j,FEdx[k],FEdy[k],var);
			}
			fprintf(joutfile,"			\n");
		}


		//-------------------- Jacobian Non-Zeros  --------------------
		status = CodeGen_JacNNZ(joutfile, FEn, maxderiv, i);
		//---------------------------------------------------


		fprintf(joutfile,"			bvec[k +pn]  = sgn*(FE%i);\n",i);

		fprintf(joutfile,"			Dxvec[k +pn] = ");
		for (j = 0; j < FEn; ++j)
		{
			for (k = 0; k <= maxderiv; ++k)
			{
				fprintf(joutfile,"dFE%id%id%i0*u%id%i0dk",i,j,k,j,k);
				if (k != maxderiv)
				{
					fprintf(joutfile," + ");
				}
			}
			if ((j+1) != FEn)
			{
				fprintf(joutfile," + ");
			}
		}
		fprintf(joutfile,";\n");

		fprintf(joutfile,"			Dyvec[k +pn] = ");
		for (j = 0; j < FEn; ++j)
		{
			for (k = 0; k <= maxderiv; ++k)
			{
				fprintf(joutfile,"dFE%id%id0%i*u%id0%idk",i,j,k,j,k);
				if (k != maxderiv)
				{
					fprintf(joutfile," + ");
				}
			}
			if ((j+1) != FEn)
			{
				fprintf(joutfile," + ");
			}
		}
		fprintf(joutfile,";\n");
		fprintf(joutfile,"			\n");

	}
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");


	//-------------------- Jac Cleanup --------------------
	status = CodeGen_JacCleanup(joutfile, FEn, maxderiv);
	//-----------------------------------------------------


	fprintf(joutfile,"	return nnztemp;\n");
	fprintf(joutfile,"}\n");
	fprintf(joutfile,"\n");


	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- INTERP Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	fprintf(joutfile,"int BVP_interp(struct param_type *params, double x[], double y[], double Psi[], int nnew, int mnew, double xnew[], double ynew[], double Psinew[])\n");
	fprintf(joutfile,"{\n");


	//-------------------- Jac Init --------------------
	status = CodeGen_JacInit(joutfile, FEn, maxderiv);
	//--------------------------------------------------


	//:::::::::::::::::::: Main Loop ::::::::::::::::::::
	fprintf(joutfile,"	//:::::::::::::::::::: Main Loop ::::::::::::::::::::\n");
	fprintf(joutfile,"	for (i = 0; i < mnew; ++i)\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		for (j = 0; j < nnew; ++j)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			k = (i*nnew+j)*p;\n");
	fprintf(joutfile,"			xk = xnew[j];\n");
	fprintf(joutfile,"			yk = ynew[i];\n");
	fprintf(joutfile,"			\n");
	fprintf(joutfile,"			I = 0; //Determine I location in old y grid\n");
	fprintf(joutfile,"			for (l = 0; l < m; ++l)\n");
	fprintf(joutfile,"			{\n");
	fprintf(joutfile,"				if (y[l] >= yk)\n");
	fprintf(joutfile,"				{\n");
	fprintf(joutfile,"					I = l;\n");
	fprintf(joutfile,"					break;\n");
	fprintf(joutfile,"				}\n");
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			J = 0; //Determine J location in old x grid\n");
	fprintf(joutfile,"			for (l = 0; l < n; ++l)\n");
	fprintf(joutfile,"			{\n");
	fprintf(joutfile,"				if (x[l] >= xk)\n");
	fprintf(joutfile,"				{\n");
	fprintf(joutfile,"					J = l;\n");
	fprintf(joutfile,"					break;\n");
	fprintf(joutfile,"				}\n");
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			\n");


	//-------------------- Newt Calc --------------------
	status = CodeGen_NewtCalc(joutfile, FEn, maxderiv);
	//---------------------------------------------------

	fprintf(joutfile,"			if (nnew == n)\n");
	fprintf(joutfile,"			{\n");
	fprintf(joutfile,"				for (b = 0; b <= r+oneyside; ++b)\n");
	fprintf(joutfile,"				{\n");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(joutfile,"					u%id0%ik += Psi[((I+sy[b])*n +J)*p +%i]*b%id%i[b];\n",j,0,j,j,0);
	}
	fprintf(joutfile,"				}\n");
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			else if (mnew == m)\n");
	fprintf(joutfile,"			{\n");
	fprintf(joutfile,"				for (a = 0; a <= r+onexside; ++a)\n");
	fprintf(joutfile,"				{\n");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(joutfile,"					u%id%i0k += Psi[(I*n +J+sx[a])*p +%i]*a%id%i[a];\n",j,0,j,j,0);
	}
	fprintf(joutfile,"				}\n");
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			else\n");
	fprintf(joutfile,"			{\n");
	fprintf(joutfile,"				printf(\"ERROR: Interp cannot interpolate x and y simultaneously.\\n\");\n");
	fprintf(joutfile,"				return -1;\n");
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			\n");

	//==================== Field Calculation ====================
	for (i = 0; i < FEn; ++i)
	{
		fprintf(joutfile,"			//FE%i\n",i);
		fprintf(joutfile,"			pn = %i;\n",i);
		fprintf(joutfile,"			Psinew[k +pn] = u%id00k;\n",i);
		fprintf(joutfile,"			\n");

	}
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");


	//-------------------- Jac Cleanup --------------------
	status = CodeGen_JacCleanup(joutfile, FEn, maxderiv);
	//-----------------------------------------------------


	fprintf(joutfile,"	return 0;\n");
	fprintf(joutfile,"}\n");
	fprintf(joutfile,"\n");
	fprintf(joutfile,"\n");
	fprintf(joutfile,"\n");


	//==================== Function Definitions ====================
	fprintf(joutfile,"//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
	fprintf(joutfile,"//::::::::::::::::::::::::::::::::::::::: Function Definitions :::::::::::::::::::::::::::::::::::::::\n");
	fprintf(joutfile,"//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
	fprintf(joutfile,"\n");


	//---------- NP Calc Definition ----------
	fprintf(joutfile,"int BVP_NPcalc(const int p, const int r, const int rn, const int I, int s[], double x[], const int dim, const int IJ, double y[], const int pn,");
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(joutfile," double Q%i_a[],",k);
	}
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(joutfile," double *d%iydk,",k);
	}
	fprintf(joutfile," const double xk)\n");
	fprintf(joutfile,"{\n");
	fprintf(joutfile,"	int k,t,q,l,a,b,c;\n");
	fprintf(joutfile,"	double s_ba;\n");
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(joutfile,"	double P%i_b;\n",k);
	}
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	double temp;\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	for (a = 0; a <= rn; ++a) //[2] loop\n");
	fprintf(joutfile,"	{\n");
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(joutfile,"		Q%i_a[a] = 0.0;\n",k);
	}
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"		for (b = a; b <= rn; ++b)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			s_ba = 1.0;\n");
	fprintf(joutfile,"			for (c = 0; c <= b; ++c) //[3] loop\n");
	fprintf(joutfile,"			{\n");
	fprintf(joutfile,"				if (c != a)\n");
	fprintf(joutfile,"				{\n");
	fprintf(joutfile,"					s_ba /= (x[I+s[a]]-x[I+s[c]]); //[3]\n");
	fprintf(joutfile,"				}\n");
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			\n");
	//----- P -----
	if (maxderiv >= 0)
	{
		fprintf(joutfile,"			//P\n");
		fprintf(joutfile,"			P0_b  = 0.0;\n");
		fprintf(joutfile,"			temp = 1.0;\n");
		fprintf(joutfile,"			for (c = 0; c <= b-1; ++c) //[4] loop\n");
		fprintf(joutfile,"			{\n");
		fprintf(joutfile,"				temp *= (xk - x[I+s[c]]); //[4]\n");
		fprintf(joutfile,"			}\n");
		fprintf(joutfile,"			P0_b += temp;\n");
		fprintf(joutfile,"			\n");
	}
	//----- P' -----
	if (maxderiv >= 1)
	{
		fprintf(joutfile,"			//P'\n");
		fprintf(joutfile,"			P1_b = 0.0;\n");
		fprintf(joutfile,"			for (q = 0; q <= b-1; ++q) //[6] loop\n");
		fprintf(joutfile,"			{\n");
		fprintf(joutfile,"				temp = 1.0;\n");
		fprintf(joutfile,"				for (c = 0; c <= b-1; ++c)\n");
		fprintf(joutfile,"				{\n");
		fprintf(joutfile,"					if (c != q)\n");
		fprintf(joutfile,"					{\n");
		fprintf(joutfile,"						temp *= (xk - x[I+s[c]]); //[6]\n");
		fprintf(joutfile,"					}\n");
		fprintf(joutfile,"				}\n");
		fprintf(joutfile,"				P1_b += temp; //[5]\n");
		fprintf(joutfile,"			}\n");
		fprintf(joutfile,"			\n");
	}
	//----- P'' -----
	if (maxderiv >= 2)
	{
		fprintf(joutfile,"			//P''\n");
		fprintf(joutfile,"			P2_b = 0.0;\n");
		fprintf(joutfile,"			for (t = 0; t <= b-1; ++t) //[7] loop\n");
		fprintf(joutfile,"			{\n");
		fprintf(joutfile,"				for (q = 0; q <= b-1; ++q)\n");
		fprintf(joutfile,"				{\n");
		fprintf(joutfile,"					if (q != t)\n");
		fprintf(joutfile,"					{\n");
		fprintf(joutfile,"						temp = 1.0;\n");
		fprintf(joutfile,"						for (c = 0; c <= b-1; ++c)\n");
		fprintf(joutfile,"						{\n");
		fprintf(joutfile,"							if ((c != q) && (c != t))\n");
		fprintf(joutfile,"							{\n");
		fprintf(joutfile,"								temp *= (xk - x[I+s[c]]); //[7]\n");
		fprintf(joutfile,"							}\n");
		fprintf(joutfile,"						}\n");
		fprintf(joutfile,"						P2_b += temp; //[8]\n");
		fprintf(joutfile,"					}\n");
		fprintf(joutfile,"				}\n");
		fprintf(joutfile,"			}\n");
		fprintf(joutfile,"			\n");
	}
	//----- P''' -----
	if (maxderiv >= 3)
	{
		fprintf(joutfile,"			//P'''\n");
		fprintf(joutfile,"			P3_b = 0.0;\n");
		fprintf(joutfile,"			for (l = 0; l <= b-1; ++l)\n");
		fprintf(joutfile,"			{\n");
		fprintf(joutfile,"				for (t = 0; t <= b-1; ++t)\n");
		fprintf(joutfile,"				{\n");
		fprintf(joutfile,"					if (t != l)\n");
		fprintf(joutfile,"					{\n");
		fprintf(joutfile,"						for (q = 0; q <= b-1; ++q)\n");
		fprintf(joutfile,"						{\n");
		fprintf(joutfile,"							if ((q != t) && (q != l))\n");
		fprintf(joutfile,"							{\n");
		fprintf(joutfile,"								temp = 1.0;\n");
		fprintf(joutfile,"								for (c = 0; c <= b-1; ++c)\n");
		fprintf(joutfile,"								{\n");
		fprintf(joutfile,"									if ((c != q) && (c != t) && (c != l))\n");
		fprintf(joutfile,"									{\n");
		fprintf(joutfile,"										temp *= (xk - x[I+s[c]]);\n");
		fprintf(joutfile,"									}\n");
		fprintf(joutfile,"								}\n");
		fprintf(joutfile,"								P3_b += temp;\n");
		fprintf(joutfile,"							}\n");
		fprintf(joutfile,"						}\n");
		fprintf(joutfile,"					}\n");
		fprintf(joutfile,"				}\n");
		fprintf(joutfile,"			}\n");
		fprintf(joutfile,"			\n");
	}
	fprintf(joutfile,"			if (b > r)\n");
	fprintf(joutfile,"			{\n");
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(joutfile,"				*d%iydk += y[IJ + s[a]*p*dim]*s_ba*P%i_b; //[1]\n",k,k);
	}
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"			else\n");
	fprintf(joutfile,"			{\n");
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(joutfile,"				Q%i_a[a] += s_ba*P%i_b; //[2]\n",k,k);
	}
	fprintf(joutfile,"			}\n");
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	return 0;\n");
	fprintf(joutfile,"}\n");
	fprintf(joutfile,"\n");


	//----- steninit Definition -----
	fprintf(joutfile,"int steninit(int s[], int *oneside, int *twoside, const int k, const int r, const int N)\n");
	fprintf(joutfile,"{\n");
	fprintf(joutfile,"	int l;\n");
	fprintf(joutfile,"	s[0] = 0; //[1]\n");
	fprintf(joutfile,"	s[(r+2)+1] = 0; //[1]\n");
	fprintf(joutfile,"	for (l = 0; l < (r+2)/2; ++l) //[2]\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		s[2*l+1] = l+1;\n");
	fprintf(joutfile,"		s[2*l+2] = -(l+1);\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	if (k < r/2) //[3]\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		*twoside = 0;\n");
	fprintf(joutfile,"		*oneside = 1;\n");
	fprintf(joutfile,"		for (l = 2*k+2; l <= (r+2)+1; ++l)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			s[l] = s[l-1] + 1;\n");
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	else if (k < (r+2)/2) //[4]\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		*twoside = 1;\n");
	fprintf(joutfile,"		*oneside = 0;\n");
	fprintf(joutfile,"		for (l = 2*k+2; l <= (r+2)+1; ++l)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			s[l] = s[l-1] + 1;\n");
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	else if ((N-1) - k < r/2) //[5]\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		*twoside = 0;\n");
	fprintf(joutfile,"		*oneside = 1;\n");
	fprintf(joutfile,"		for (l = 2*((N-1)-k)+1; l <= (r+2)+1; ++l)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			s[l] = s[l-1] - 1;\n");
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	else if ((N-1) - k < (r+2)/2) //[6]\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		*twoside = 1;\n");
	fprintf(joutfile,"		*oneside = 0;\n");
	fprintf(joutfile,"		for (l = 2*((N-1)-k)+1; l <= (r+2)+1; ++l)\n");
	fprintf(joutfile,"		{\n");
	fprintf(joutfile,"			s[l] = s[l-1] - 1;\n");
	fprintf(joutfile,"		}\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	/*\n");
	fprintf(joutfile,"	printf(\"k = %%03i. oneside = %%i. twoside = %%i. s[] = \",k,*oneside,*twoside); //[7]\n");
	fprintf(joutfile,"	for (l = 0; l < (r+2)+1; ++l)\n");
	fprintf(joutfile,"	{\n");
	fprintf(joutfile,"		printf(\" %%2i,\",s[l]);\n");
	fprintf(joutfile,"	}\n");
	fprintf(joutfile,"	printf(\" %%2i.\\n\",s[(r+2)+1]);\n");
	fprintf(joutfile,"	*/\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	\n");
	fprintf(joutfile,"	return 0;\n");
	fprintf(joutfile,"}\n");
	fprintf(joutfile,"\n");
	fprintf(joutfile,"\n");
	fprintf(joutfile,"\n");
	fprintf(joutfile,"\n");


	//Cleanup
	fclose(joutfile);



	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- ICBC Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sprintf(iboutname,"./Code/BVP_GENICBC.c");
	iboutfile = fopen(iboutname,"w");
	if(iboutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file

	fprintf(iboutfile,"#include \"BVP_header.h\"\n\n");

	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- BC Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------

	fprintf(iboutfile,"int BVP_GENBC(struct param_type *params, double x[], double y[], double Psi[], double JacVal[], int JacRow[], int JacCol[], double bvec[], double Dxvec[], double Dyvec[])\n");
	fprintf(iboutfile,"{\n");


	//-------------------- Jac Init --------------------
	status = CodeGen_JacInit(iboutfile, FEn, maxderiv);
	//--------------------------------------------------


	fprintf(iboutfile,"	//----- Boundary partials -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(iboutfile,"	//BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			//----- X0 -----
			fprintf(iboutfile,"	double BC%i%c%i;\n",i,BCdim[bc],BCloc[bc]);
			for (j = 0; j < FEn; ++j)
			{
				fprintf(iboutfile,"	double dBC%i%c%id%id00, dBC%i%c%id%id10, dBC%i%c%id%id01;\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j);
			}
			fprintf(iboutfile,"	\n");
		}
	}

	fprintf(iboutfile,"	nnztemp = 0;\n");
	fprintf(iboutfile,"	\n");
	fprintf(iboutfile,"	\n");

	//----- String variable "varBC" -----
	strcpy(varBC,"(xk, yk, r_H, M, alpha, beta,");
	for (j = 0 ; j < FEn; ++j)
	{
		sprintf(vartemp," u%id00k, u%id10k, u%id01k",j,j,j);
		strcat(varBC, vartemp);

		if ((j+1) != FEn)
		{
			strcat(varBC,",");
		}
	}
	strcat(varBC,")");


	//:::::::::::::::::::: 4 Boundaries ::::::::::::::::::::
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(iboutfile,"	//:::::::::::::::::::: %c%i ::::::::::::::::::::\n",BCdim[bc],BCloc[bc]);


		//-------------------- BC Domain --------------------
		status =  CodeGen_BCDomain(iboutfile, bc);
		//---------------------------------------------------


		fprintf(iboutfile,"			k = (i*n+j)*p;\n");
		fprintf(iboutfile,"			I = i, J = j;\n");
		fprintf(iboutfile,"			xk = x[j];\n");
		fprintf(iboutfile,"			yk = y[i];\n");
		fprintf(iboutfile,"			\n");


		//-------------------- Newt Calc --------------------
		status = CodeGen_NewtCalc(iboutfile, FEn, maxderiv);
		//---------------------------------------------------


		fprintf(iboutfile,"			for (a = 0; a <= r+onexside; ++a)\n");
		fprintf(iboutfile,"			{\n");

		for (j = 0; j < FEn; ++j)
		{
			fprintf(iboutfile,"				u%id%i0k += Psi[(I*n +J+sx[a])*p +%i]*a%id%i[a];\n",j,0,j,j,0);
		}
		fprintf(iboutfile,"			}\n");
		fprintf(iboutfile,"			\n");

		fprintf(iboutfile,"			//---------- Evaluate partial derivatives and structure for Jacobian ----------\n");
		fprintf(iboutfile,"			\n");

		for (i = 0; i < FEn; ++i)
		{
			fprintf(iboutfile,"			//----- BC%i -----\n",i);
			fprintf(iboutfile,"			BC%i%c%i       =       BC%i%c%iout%s;\n",i,BCdim[bc],BCloc[bc],i,BCdim[bc],BCloc[bc],varBC);
			fprintf(iboutfile,"			\n");

			for (j = 0; j < FEn; ++j)
			{
				fprintf(iboutfile,"			dBC%i%c%id%id00 = dBC%i%c%id%id00out%s;\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j,varBC);
				fprintf(iboutfile,"			dBC%i%c%id%id10 = dBC%i%c%id%id10out%s;\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j,varBC);
				fprintf(iboutfile,"			dBC%i%c%id%id01 = dBC%i%c%id%id01out%s;\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j,varBC);
				fprintf(iboutfile,"			\n");
			}


			//-------------------- Boundary Non-Zeros  --------------------
			status = CodeGen_BCNNZ(iboutfile, FEn, maxderiv, i, bc);
			//---------------------------------------------------


			fprintf(iboutfile,"			bvec[k +pn]  = sgn*(BC%i%c%i);\n",i,BCdim[bc],BCloc[bc]);
			fprintf(iboutfile,"			Dxvec[k +pn] = ");
			for (j = 0; j < FEn; ++j)
			{
				fprintf(iboutfile,"dBC%i%c%id%id10*u%id10dk + dBC%i%c%id%id00*u%id00dk",i,BCdim[bc],BCloc[bc],j,j,i,BCdim[bc],BCloc[bc],j,j);
				if ((j+1) != FEn)
				{
					fprintf(iboutfile," + ");
				}
			}
			fprintf(iboutfile,";\n");
			fprintf(iboutfile,"			Dyvec[k +pn] = ");
			for (j = 0; j < FEn; ++j)
			{
				fprintf(iboutfile,"dBC%i%c%id%id01*u%id01dk + dBC%i%c%id%id00*u%id00dk",i,BCdim[bc],BCloc[bc],j,j,i,BCdim[bc],BCloc[bc],j,j);
				if ((j+1) != FEn)
				{
					fprintf(iboutfile," + ");
				}
			}
			fprintf(iboutfile,";\n");
			fprintf(iboutfile,"			\n");
		}
		fprintf(iboutfile,"		}\n");
		fprintf(iboutfile,"	}\n");
		fprintf(iboutfile,"	\n");
	}
	fprintf(iboutfile,"	\n");


	//-------------------- Jac Cleanup --------------------
	status = CodeGen_JacCleanup(iboutfile, FEn, maxderiv);
	//-----------------------------------------------------


	fprintf(iboutfile,"	return nnztemp;\n");
	fprintf(iboutfile,"}\n");
	fprintf(iboutfile,"\n");
	fprintf(iboutfile,"\n");
	fprintf(iboutfile,"\n");



	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- IC Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
/*
	fprintf(iboutfile,"int BVP_GENIC(struct param_type *params, struct ICparam_type *ICparams, double x[], double y[], double Psi[], double GRIDdx[], double GRIDdy[])\n");
	fprintf(iboutfile,"{\n");
	fprintf(iboutfile,"	int i,j,k,l, status;\n");
	fprintf(iboutfile,"	const int n = (*params).nparam;\n");
	fprintf(iboutfile,"	const int m = (*params).mparam;\n");
	fprintf(iboutfile,"	const int p = (*params).pparam;\n");
	fprintf(iboutfile,"	const int N = (*params).Nparam;\n");
	fprintf(iboutfile,"	const int r = (*params).rparam;\n");
	fprintf(iboutfile,"	const double r_H = (*ICparams).r_HICparam;\n");
	fprintf(iboutfile,"	const double M = (*ICparams).MICparam;\n");
	fprintf(iboutfile,"	const double alpha = (*ICparams).alphaICparam;\n");
	fprintf(iboutfile,"	const double beta = (*ICparams).betaICparam;\n");
	fprintf(iboutfile,"	const double delta = (*ICparams).deltaICparam;\n");
	fprintf(iboutfile,"	int pn;\n");
	fprintf(iboutfile,"	double xk, yk;\n");
	fprintf(iboutfile,"	\n");

	fprintf(iboutfile,"	//Generate IC\n");
	fprintf(iboutfile,"	for (i = 0; i < m; ++i)\n");
	fprintf(iboutfile,"	{\n");
	fprintf(iboutfile,"		y[i] = ((double) i/(m-1)) * (M_PI/2.0); //0 < y < M_PI/2\n");
	//fprintf(iboutfile,"		y[i] = (double) i/(m-1); //0 < y < 1\n");
	fprintf(iboutfile,"		yk = y[i];\n");
	fprintf(iboutfile,"		for (j = 0; j < n; ++j)\n");
	fprintf(iboutfile,"		{\n");
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
//Choose what power of x to make grid size
	//fprintf(iboutfile,"			x[j] = pow((double) j/(n-1),2.0);\n");
	fprintf(iboutfile,"			x[j] = ((double) j/(n-1));\n");
//
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	fprintf(iboutfile,"			xk = x[j];\n");
	fprintf(iboutfile,"			k = (i*n+j)*p;\n");
	fprintf(iboutfile,"			\n");

	for (j = 0; j < FEn; ++j)
	{
		fprintf(iboutfile,"			pn = %i;\n",j);
		fprintf(iboutfile,"			Psi[k +pn] = IC%iout(xk, yk, r_H, M, alpha, beta, delta);\n",j);
		fprintf(iboutfile,"			\n");
	}

	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"	}\n");
	fprintf(iboutfile,"	\n");

	fprintf(iboutfile,"	//Calculate uniform grid x dimension\n");
	fprintf(iboutfile,"	for (j = 0; j < n; ++j)\n");
	fprintf(iboutfile,"	{\n");
	fprintf(iboutfile,"		if (j == 0)\n");
	fprintf(iboutfile,"		{\n");
	fprintf(iboutfile,"			GRIDdx[j] = (x[j+1] - x[j]);\n");
	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"		else if (j == n-1)\n");
	fprintf(iboutfile,"		{\n");
	fprintf(iboutfile,"			GRIDdx[j] = (x[j] - x[j-1]);\n");
	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"		else\n");
	fprintf(iboutfile,"		{\n");
	fprintf(iboutfile,"			GRIDdx[j] = (x[j+1] - x[j-1])/2.0;\n");
	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"	}\n");

	fprintf(iboutfile,"	//Calculate uniform grid y dimension\n");
	fprintf(iboutfile,"	for (i = 0; i < m; ++i)\n");
	fprintf(iboutfile,"	{\n");
	fprintf(iboutfile,"		if (i == 0)\n");
	fprintf(iboutfile,"		{\n");
	fprintf(iboutfile,"			GRIDdy[i] = (y[i+1] - y[i]);\n");
	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"		else if (i == m-1)\n");
	fprintf(iboutfile,"		{\n");
	fprintf(iboutfile,"			GRIDdy[i] = (y[i] - y[i-1]);\n");
	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"		else\n");
	fprintf(iboutfile,"		{\n");
	fprintf(iboutfile,"			GRIDdy[i] = (y[i+1] - y[i-1])/2.0;\n");
	fprintf(iboutfile,"		}\n");
	fprintf(iboutfile,"	}\n");

	fprintf(iboutfile,"	\n");
	fprintf(iboutfile,"	\n");

	fprintf(iboutfile,"	return 0;\n");
	fprintf(iboutfile,"}\n");
	fprintf(iboutfile,"\n");
	fprintf(iboutfile,"\n");
	fprintf(iboutfile,"\n");
	fprintf(iboutfile,"\n");
*/

	//Cleanup
	fclose(iboutfile);
	free(FEdx);
	free(FEdy);


	return 0;
}



//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::: Function Definitions :::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int CodeGen_JacInit(FILE *outfile, const int FEn, const int maxderiv)
{
	int i,j,k,l, status;
	fprintf(outfile,"	int i,j,k,l, status;\n");
	fprintf(outfile,"	const int n = (*params).nparam;\n");
	fprintf(outfile,"	const int m = (*params).mparam;\n");
	fprintf(outfile,"	const int p = (*params).pparam;\n");
	fprintf(outfile,"	const int N = (*params).Nparam;\n");
	fprintf(outfile,"	const int r = (*params).rparam;\n");
	fprintf(outfile,"	const int nnzJac = (*params).nnzparam;\n");
	fprintf(outfile,"	const double r_H = (*params).r_Hparam;\n");
	fprintf(outfile,"	const double M = (*params).Mparam;\n");
	fprintf(outfile,"	const double alpha = (*params).alphaparam;\n");
	fprintf(outfile,"	const double beta = (*params).betaparam;\n");
	fprintf(outfile,"	const size_t slen = r+4;\n");
	fprintf(outfile,"	int pn;\n");
	fprintf(outfile,"	const double sgn = -1.0;\n");
	fprintf(outfile,"	int a, b; //a for x, b for y(theta)\n");
	fprintf(outfile,"	int I, onexside, twoxside;\n");
	fprintf(outfile,"	int J, oneyside, twoyside;\n");
	fprintf(outfile,"	double xk, yk;\n");
	fprintf(outfile,"	const int dimx = 1;\n");
	fprintf(outfile,"	const int dimy = n;\n");
	fprintf(outfile,"	int nnztemp;\n");
	fprintf(outfile,"	\n");
	fprintf(outfile,"	int *sx;\n");
	fprintf(outfile,"	sx = (int*)calloc(slen, sizeof(int));\n");
	fprintf(outfile,"	int *sy;\n");
	fprintf(outfile,"	sy = (int*)calloc(slen, sizeof(int));\n");
	fprintf(outfile,"	\n");

	fprintf(outfile,"	//----- Finite difference stencil -----\n");
	for (j = 0 ; j < FEn; ++j)
	{
		fprintf(outfile,"	double");
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile," *a%id%i",j,k);
			if (k != maxderiv)
			{
				fprintf(outfile,",");
			}
		}
		fprintf(outfile,";\n");
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile,"	a%id%i = (double*)calloc(slen, sizeof(double));\n",j,k);
		}

		fprintf(outfile,"	double");
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile," *b%id%i",j,k);
			if (k != maxderiv)
			{
				fprintf(outfile,",");
			}
		}
		fprintf(outfile,";\n");
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile,"	b%id%i = (double*)calloc(slen, sizeof(double));\n",j,k);
		}
		fprintf(outfile,"	\n");
	}

	fprintf(outfile,"	//----- Field derivatives initialization -----\n");
	for (j = 0 ; j < FEn; ++j)
	{
		fprintf(outfile,"	double u%id00k, ",j);
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id%i0k, ",j,k);
		}
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id0%ik",j,k);
			if (k != maxderiv)
			{
				fprintf(outfile,", ");
			}
		}
		fprintf(outfile,";\n");

		fprintf(outfile,"	double u%id00dk,",j);
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id%i0dk,",j,k);
		}
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id0%idk",j,k);
			if (k != maxderiv)
			{
				fprintf(outfile,",");
			}
		}
		fprintf(outfile,";\n");
		fprintf(outfile,"	\n");
	}


	return 0;
}

int CodeGen_NewtCalc(FILE *outfile, int FEn, int maxderiv)
{
	int i,j,k,l, status;
	fprintf(outfile,"			//Zero derivatives of each field\n");
	fprintf(outfile,"			onexside = 0; twoxside = 0;\n");
	fprintf(outfile,"			oneyside = 0; twoyside = 0;\n");
	for (j = 0 ; j < FEn; ++j)
	{
		fprintf(outfile,"			u%id00k  = 0.0,",j);
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id%i0k  = 0.0,",j,k);
		}
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id0%ik  = 0.0",j,k);
			if (k != maxderiv)
			{
				fprintf(outfile,",");
			}
		}
		fprintf(outfile,";\n");

		fprintf(outfile,"			u%id00dk = 0.0,",j);
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id%i0dk = 0.0,",j,k);
		}
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile," u%id0%idk = 0.0",j,k);
			if (k != maxderiv)
			{
				fprintf(outfile,",");
			}
		}
		fprintf(outfile,";\n");
		fprintf(outfile,"			\n");
	}
	fprintf(outfile,"			//Compute stencil\n");
	fprintf(outfile,"			status = steninit(sy, &oneyside, &twoyside, I, r, m);\n");
	fprintf(outfile,"			status = steninit(sx, &onexside, &twoxside, J, r, n);\n");
	fprintf(outfile,"			\n");

	fprintf(outfile,"			//Compute Newton polynomial representation for each field\n");
	for (j = 0 ; j < FEn; ++j)
	{
		//NPcalc calculates to max derivative
		fprintf(outfile,"			status = BVP_NPcalc(p, r+oneyside, r+oneyside+twoyside+2, I, sy, y, dimy, (I*n+J)*p+%i, Psi, %i,",j,j);
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile," b%id%i,",j,k);
		}
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile," &u%id0%idk,",j,k);
		}
		fprintf(outfile," yk);\n");

		fprintf(outfile,"			status = BVP_NPcalc(p, r+onexside, r+onexside+twoxside+2, J, sx, x, dimx, (I*n+J)*p+%i, Psi, %i,",j,j);
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile," a%id%i,",j,k);
		}
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile," &u%id%i0dk,",j,k);
		}
		fprintf(outfile," xk);\n");

		fprintf(outfile,"			\n");
	}

	fprintf(outfile,"			//Calculate derivative of each field\n");
	fprintf(outfile,"			for (b = 0; b <= r+oneyside; ++b)\n");
	fprintf(outfile,"			{\n");

	for (j = 0; j < FEn; ++j)
	{
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile,"				u%id0%ik += Psi[((I+sy[b])*n +J)*p +%i]*b%id%i[b];\n",j,k,j,j,k);
		}
		fprintf(outfile,"				\n");
	}
	fprintf(outfile,"			}\n");

	fprintf(outfile,"			for (a = 0; a <= r+onexside; ++a)\n");
	fprintf(outfile,"			{\n");

	for (j = 0; j < FEn; ++j)
	{
		for (k = 1; k <= maxderiv; ++k)
		{
			fprintf(outfile,"				u%id%i0k += Psi[(I*n +J+sx[a])*p +%i]*a%id%i[a];\n",j,k,j,j,k);
		}
		fprintf(outfile,"				\n");
	}
	fprintf(outfile,"			}\n");
	fprintf(outfile,"			\n");


	return 0;
}

int CodeGen_JacCleanup(FILE *outfile, int FEn, int maxderiv)
{
	int i,j,k,l, status;
	fprintf(outfile,"	\n");
	fprintf(outfile,"	//Cleanup\n");
	fprintf(outfile,"	free(sx);\n");
	fprintf(outfile,"	free(sy);\n");
	fprintf(outfile,"	\n");

	for (j = 0 ; j < FEn; ++j)
	{
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile,"	free(a%id%i);\n",j,k);
		}
		for (k = 0; k <= maxderiv; ++k)
		{
			fprintf(outfile,"	free(b%id%i);\n",j,k);
		}
		fprintf(outfile,"	\n");
	}
	fprintf(outfile,"	\n");


	return 0;
}

int CodeGen_BCDomain(FILE *outfile, int bc)
{
	if (bc == 0)
	{
		//X0
		fprintf(outfile,"	for (i = 0; i < m; ++i)\n");
		fprintf(outfile,"	{\n");
		fprintf(outfile,"		j = 0;\n");
		fprintf(outfile,"		{\n");
	}
	else if (bc == 1)
	{
		//X1
		fprintf(outfile,"	for (i = 0; i < m; ++i)\n");
		fprintf(outfile,"	{\n");
		fprintf(outfile,"		j = n-1;\n");
		fprintf(outfile,"		{\n");
	}
	else if (bc == 2)
	{
		//Y1
		fprintf(outfile,"	i = 0;\n");
		fprintf(outfile,"	{\n");
		fprintf(outfile,"		for (j = 1; j < n-1; ++j)\n"); //Don't double count corners
		fprintf(outfile,"		{\n");
	}
	else if (bc == 3)
	{
		//Y0
		fprintf(outfile,"	i = m-1;\n");
		fprintf(outfile,"	{\n");
		fprintf(outfile,"		for (j = 1; j < n-1; ++j)\n");
		fprintf(outfile,"		{\n");
	}


	return 0;
}
