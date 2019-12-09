#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int i,j,k,l,bc, status;
	char nnzoutname[4096];
	FILE *nnzoutfile;

	const int FEn = atoi(argv[1]);
	const int maxderiv = atoi(argv[2]);
	const char BCdim[5] = "XXYY";
	const int BCloc[4] = {0,1,1,0};
	const int BClength = 4;
	size_t FEdlength;


	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- JacNNZ Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------

	sprintf(nnzoutname,"./Funcs/CodeGen_GENNNZ.c");
	nnzoutfile = fopen(nnzoutname,"w");
	if(nnzoutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file


	fprintf(nnzoutfile,"#include \"Decl_JSFuncs.c\"\n\n");
	fprintf(nnzoutfile,"int CodeGen_JacNNZ(FILE *outfile, int FEn, int maxderiv, int i)\n");
	fprintf(nnzoutfile,"{\n");
	fprintf(nnzoutfile,"	int j,k, status;\n");
	fprintf(nnzoutfile,"	\n");

	fprintf(nnzoutfile,"	//----- Jacobian structure -----\n");
	for (i = 0; i < FEn; ++i)
	{
		for (j = 0; j < FEn; ++j)
		{
			fprintf(nnzoutfile,"	int dFE%id%idX0, dFE%id%id0X;\n",i,j,i,j);
		}
		fprintf(nnzoutfile,"	\n");
	}

	for (i = 0; i < FEn; ++i)
	{
		fprintf(nnzoutfile,"	//-------------------- i = %i --------------------\n",i);
		if (i == 0)
		{
			fprintf(nnzoutfile,"	if (i == %i)\n",i);
		}
		else
		{
			fprintf(nnzoutfile,"	else if (i == %i)\n",i);
		}
		fprintf(nnzoutfile,"	{\n");

		for (j = 0; j < FEn; ++j)
		{
			fprintf(nnzoutfile,"		dFE%id%idX0 = dFE%id%idX0out();\n",i,j,i,j);
			fprintf(nnzoutfile,"		dFE%id%id0X = dFE%id%id0Xout();\n",i,j,i,j);
			fprintf(nnzoutfile,"		\n");
		}

		fprintf(nnzoutfile,"		fprintf(outfile,\"			//-------------------------------------------------------\\n\");\n");
		fprintf(nnzoutfile,"		fprintf(outfile,\"			pn = %%i;\\n\",%i);\n",i);
		fprintf(nnzoutfile,"		fprintf(outfile,\"			\\n\");\n");
		fprintf(nnzoutfile,"		\n");

		fprintf(nnzoutfile,"		fprintf(outfile,\"			for (a = 0; a <= r+onexside; ++a)\\n\");\n");
		fprintf(nnzoutfile,"		fprintf(outfile,\"			{\\n\");\n");
		for (j = 0; j < FEn; ++j)
		{
			fprintf(nnzoutfile,"		if (dFE%id%idX0 != 0)\n",i,j);
			fprintf(nnzoutfile,"		{\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				JacVal[nnztemp] = \");\n");
			fprintf(nnzoutfile,"			for (k = 0; k <= maxderiv; ++k)\n");
			fprintf(nnzoutfile,"			{\n");
			fprintf(nnzoutfile,"				fprintf(outfile,\"dFE%%id%%id%%i0*a%%id%%i[a]\",%i,%i,k,%i,k);\n",i,j,j);
			fprintf(nnzoutfile,"				if (k != maxderiv)\n");
			fprintf(nnzoutfile,"				{\n");
			fprintf(nnzoutfile,"					fprintf(outfile,\" + \");\n");
			fprintf(nnzoutfile,"				}\n");
			fprintf(nnzoutfile,"			}\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\";\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				JacRow[nnztemp] = (k +pn);\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				JacCol[nnztemp] = (I*n +J+sx[a])*p +%%i;\\n\",%i);\n",j);
			fprintf(nnzoutfile,"			fprintf(outfile,\"				nnztemp += 1;\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				\\n\");\n");
			fprintf(nnzoutfile,"		}\n");
		}
		fprintf(nnzoutfile,"		fprintf(outfile,\"			}\\n\");\n");
		fprintf(nnzoutfile,"		\n");

		fprintf(nnzoutfile,"		fprintf(outfile,\"			for (b = 0; b <= r+oneyside; ++b)\\n\");\n");
		fprintf(nnzoutfile,"		fprintf(outfile,\"			{\\n\");\n");
		for (j = 0; j < FEn; ++j)
		{
			fprintf(nnzoutfile,"		if (dFE%id%id0X != 0)\n",i,j);
			fprintf(nnzoutfile,"		{\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				JacVal[nnztemp] = \");\n");
			fprintf(nnzoutfile,"			for (k = 1; k <= maxderiv; ++k)\n");
			fprintf(nnzoutfile,"			{\n");
			fprintf(nnzoutfile,"				fprintf(outfile,\"dFE%%id%%id0%%i*b%%id%%i[b]\",%i,%i,k,%i,k);\n",i,j,j);
			fprintf(nnzoutfile,"				if (k != maxderiv)\n");
			fprintf(nnzoutfile,"				{\n");
			fprintf(nnzoutfile,"					fprintf(outfile,\" + \");\n");
			fprintf(nnzoutfile,"				}\n");
			fprintf(nnzoutfile,"			}\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\";\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				JacRow[nnztemp] = (k +pn);\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +%%i;\\n\",%i);\n",j);
			fprintf(nnzoutfile,"			fprintf(outfile,\"				nnztemp += 1;\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"				\\n\");\n");
			fprintf(nnzoutfile,"		}\n");
		}
		fprintf(nnzoutfile,"		fprintf(outfile,\"			}\\n\");\n");
		fprintf(nnzoutfile,"		fprintf(outfile,\"			\\n\");\n");
		fprintf(nnzoutfile,"		\n");
		fprintf(nnzoutfile,"	}\n");

	}

	fprintf(nnzoutfile,"	else\n");
	fprintf(nnzoutfile,"	{\n");
	fprintf(nnzoutfile,"		printf(\"ERROR!! In CodeGen_NNZ.c, %%i does not match number of field equations.\\n\",i);\n");
	fprintf(nnzoutfile,"		return -1;\n");
	fprintf(nnzoutfile,"	}\n");
	fprintf(nnzoutfile,"	\n");
	fprintf(nnzoutfile,"	\n");
	fprintf(nnzoutfile,"	return 0;\n");
	fprintf(nnzoutfile,"}\n");
	fprintf(nnzoutfile,"\n");



	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- Boundary NNZ Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------


	fprintf(nnzoutfile,"int CodeGen_BCNNZ(FILE *outfile, int FEn, int maxderiv, int i, int bc)\n");
	fprintf(nnzoutfile,"{\n");
	fprintf(nnzoutfile,"	int j,k, status;\n");
	fprintf(nnzoutfile,"	const char BCdim[5] = \"XXYY\";\n");
	fprintf(nnzoutfile,"	const int BCloc[4] = {0,1,1,0};\n");
	fprintf(nnzoutfile,"	const int BClength = 4;\n");
	fprintf(nnzoutfile,"	\n");

	fprintf(nnzoutfile,"	//----- Boundary structure -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(nnzoutfile,"	//BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			for (j = 0; j < FEn; ++j)
			{
				fprintf(nnzoutfile,"	int dBC%i%c%id%idX0, dBC%i%c%id%id0X;\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j);
			}
			fprintf(nnzoutfile,"	\n");
		}
	}

	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(nnzoutfile,"	//:::::::::::::::::::: %c%i ::::::::::::::::::::\n",BCdim[bc],BCloc[bc]);
		if (bc == 0)
		{
			fprintf(nnzoutfile,"	if (bc == %i)\n",bc);
		}
		else
		{
			fprintf(nnzoutfile,"	else if (bc == %i)\n",bc);
		}
		fprintf(nnzoutfile,"	{\n");

		for (i = 0; i < FEn; ++i)
		{
			fprintf(nnzoutfile,"		//-------------------- i = %i --------------------\n",i);
			if (i == 0)
			{
				fprintf(nnzoutfile,"		if (i == %i)\n",i);
			}
			else
			{
				fprintf(nnzoutfile,"		else if (i == %i)\n",i);
			}
			fprintf(nnzoutfile,"		{\n");

			for (j = 0; j < FEn; ++j)
			{
				fprintf(nnzoutfile,"			dBC%i%c%id%idX0 = dBC%i%c%id%idX0out();\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j);
				fprintf(nnzoutfile,"			dBC%i%c%id%id0X = dBC%i%c%id%id0Xout();\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j);
				fprintf(nnzoutfile,"			\n");
			}

			fprintf(nnzoutfile,"			fprintf(outfile,\"			//-------------------------------------------------------\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"			pn = %%i;\\n\",%i);\n",i);
			fprintf(nnzoutfile,"			fprintf(outfile,\"			\\n\");\n");
			fprintf(nnzoutfile,"			\n");

			fprintf(nnzoutfile,"			fprintf(outfile,\"			for (a = 0; a <= r+onexside; ++a)\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"			{\\n\");\n");
			for (j = 0; j < FEn; ++j)
			{
				fprintf(nnzoutfile,"			if (dBC%i%c%id%idX0 != 0)\n",i,BCdim[bc],BCloc[bc],j);
				fprintf(nnzoutfile,"			{\n");
				fprintf(nnzoutfile,"				fprintf(outfile,\"				JacVal[nnztemp] = dBC%%i%%c%%id%%id00*a%%id0[a] + dBC%%i%%c%%id%%id10*a%%id1[a];\\n\",%i,\'%c\',%i,%i,%i,%i,\'%c\',%i,%i,%i);\n",i,BCdim[bc],BCloc[bc],j,j,i,BCdim[bc],BCloc[bc],j,j);
				fprintf(nnzoutfile,"				fprintf(outfile,\"				JacRow[nnztemp] = (k +pn);\\n\");\n");
				fprintf(nnzoutfile,"				fprintf(outfile,\"				JacCol[nnztemp] = (I*n +J+sx[a])*p +%%i;\\n\",%i);\n",j);
				fprintf(nnzoutfile,"				fprintf(outfile,\"				nnztemp += 1;\\n\");\n");
				fprintf(nnzoutfile,"				fprintf(outfile,\"				\\n\");\n");
				fprintf(nnzoutfile,"			}\n");
			}
			fprintf(nnzoutfile,"			fprintf(outfile,\"			}\\n\");\n");
			fprintf(nnzoutfile,"			\n");

			fprintf(nnzoutfile,"			fprintf(outfile,\"			for (b = 0; b <= r+oneyside; ++b)\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"			{\\n\");\n");
			for (j = 0; j < FEn; ++j)
			{
				fprintf(nnzoutfile,"			if (dBC%i%c%id%id0X != 0)\n",i,BCdim[bc],BCloc[bc],j);
				fprintf(nnzoutfile,"			{\n");
				fprintf(nnzoutfile,"				fprintf(outfile,\"				JacVal[nnztemp] = dBC%%i%%c%%id%%id01*b%%id1[b];\\n\",%i,\'%c\',%i,%i,%i);\n",i,BCdim[bc],BCloc[bc],j,j);
				fprintf(nnzoutfile,"				fprintf(outfile,\"				JacRow[nnztemp] = (k +pn);\\n\");\n");
				fprintf(nnzoutfile,"				fprintf(outfile,\"				JacCol[nnztemp] = ((I+sy[b])*n +J)*p +%%i;\\n\",%i);\n",j);
				fprintf(nnzoutfile,"				fprintf(outfile,\"				nnztemp += 1;\\n\");\n");
				fprintf(nnzoutfile,"			}\n");
			}
			fprintf(nnzoutfile,"			fprintf(outfile,\"			}\\n\");\n");
			fprintf(nnzoutfile,"			fprintf(outfile,\"			\\n\");\n");
			fprintf(nnzoutfile,"			\n");
			fprintf(nnzoutfile,"		}\n");

		}

		fprintf(nnzoutfile,"		else\n");
		fprintf(nnzoutfile,"		{\n");
		fprintf(nnzoutfile,"			printf(\"ERROR!! In CodeGen_NNZ.c, %%i does not match number of field equations.\\n\",i);\n");
		fprintf(nnzoutfile,"			return -1;\n");
		fprintf(nnzoutfile,"		}\n");
		fprintf(nnzoutfile,"		\n");
		fprintf(nnzoutfile,"	}\n");

	}
	fprintf(nnzoutfile,"	else\n");
	fprintf(nnzoutfile,"	{\n");
	fprintf(nnzoutfile,"		printf(\"ERROR!! In CodeGen_NNZ.c, %%i does not match number of boundary conditions.\\n\",bc);\n");
	fprintf(nnzoutfile,"		return -1;\n");
	fprintf(nnzoutfile,"	}\n");
	fprintf(nnzoutfile,"	\n");
	fprintf(nnzoutfile,"	\n");
	fprintf(nnzoutfile,"	return 0;\n");
	fprintf(nnzoutfile,"}\n");
	fprintf(nnzoutfile,"\n");
	fprintf(nnzoutfile,"\n");
	fprintf(nnzoutfile,"\n");
	fprintf(nnzoutfile,"\n");


	//Cleanup
	fclose(nnzoutfile);


	return 0;
}
