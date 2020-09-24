#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int i,j,k,l,bc, status;
	char poutname[256];
	FILE *poutfile;
	char qoutname[256];
	FILE *qoutfile;
	char coutname[256];
	FILE *coutfile;
	char JSoutname[256];
	FILE *JSoutfile;
	char vardecl[4096];
	char vardecltemp[4096];

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
	//-------------------- MapleREAD FE Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sprintf(poutname,"./Gens/MapleREAD_FEexport.txt");
	poutfile = fopen(poutname,"w");
	if(poutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file


	//-------------------- Export Variables -------------------------
	fprintf(poutfile,"\n#:::::::::::::::::::: Export Variables ::::::::::::::::::::\n");
	fprintf(poutfile,"FEsub3export := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(poutfile,"diff(u%i(x,y),x,x,x) = u%id30, diff(u%i(x,y),y,y,y) = u%id03",j,j,j,j);
		if ((j+1) != FEn)
		{
			fprintf(poutfile,", ");
		}
	}
	fprintf(poutfile,"};\n");
	fprintf(poutfile,"FEsub2export := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(poutfile,"diff(u%i(x,y),x,x) = u%id20, diff(u%i(x,y),y,y) = u%id02",j,j,j,j);
		if ((j+1) != FEn)
		{
			fprintf(poutfile,", ");
		}
	}
	fprintf(poutfile,"};\n");
	fprintf(poutfile,"FEsub1export := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(poutfile,"diff(u%i(x,y),x) = u%id10, diff(u%i(x,y),y) = u%id01",j,j,j,j);
		if ((j+1) != FEn)
		{
			fprintf(poutfile,", ");
		}
	}
	fprintf(poutfile,"};\n");
	fprintf(poutfile,"FEsub0export := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(poutfile,"u%i(x,y) = u%id00",j,j);
		if ((j+1) != FEn)
		{
			fprintf(poutfile,", ");
		}
	}
	fprintf(poutfile,"};\n");
	fprintf(poutfile,"collectexport := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(poutfile,"u%id30, u%id03, u%id20, u%id02, u%id10, u%id01, u%id00",j,j,j,j,j,j,j);
		if ((j+1) != FEn)
		{
			fprintf(poutfile,", ");
		}
	}
	fprintf(poutfile,"};\n");
	fprintf(poutfile,"booltonum := {true=1,false=0}:\n\n");

	fprintf(poutfile,"FEvar := [x::numeric,y::numeric,r_H::numeric,M::numeric,alpha::numeric,beta::numeric,");
	for (j = 0 ; j < FEn; ++j)
	{
		for (k = 0; k < FEdlength; ++k)
		{
			fprintf(poutfile,"u%id%i%i::numeric",j,FEdx[k],FEdy[k]);
			if ((k+1) != FEdlength)
			{
				fprintf(poutfile,",");
			}
		}
		if ((j+1) != FEn)
		{
			fprintf(poutfile,",");
		}
	}
	fprintf(poutfile,"]:\n\n");


	//-------------------- Export Equations --------------------
	fprintf(poutfile,"#:::::::::::::::::::: Export Equations ::::::::::::::::::::\n");
	fprintf(poutfile,"#----- FE -----\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"FE%ifull := collect(subs(FEsub3export, FEsub2export, FEsub1export, FEsub0export, FE%i),collectexport):\n",i,i);
	}
	fprintf(poutfile,"\n");

	//-------------------- FE deriv calc -------------------------
	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"#----- dFE%i -----\n",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(poutfile,"#FE%i du%i\n",i,j);
			for (k = 0; k < FEdlength; ++k)
			{
				fprintf(poutfile,"dFE%id%id%i%ifull := collect(simplify(diff(FE%ifull,u%id%i%i)),collectexport):\n",i,j,FEdx[k],FEdy[k],i,j,FEdx[k],FEdy[k]);
			}
		}
		fprintf(poutfile,"\n");
	}

	//-------------------- FE structure calc -------------------------
	fprintf(poutfile,"#----- FE Structure-----\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"#FE%i\n",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(poutfile,"dFE%id%idX0full := subs(booltonum, evalb(",i,j);
			for (k = 0; k <= maxderiv; ++k) //ranges from 0 to maxderiv
			{
				fprintf(poutfile,"dFE%id%id%i0full <> 0",i,j,k);
				if (k != maxderiv)
				{
					fprintf(poutfile," or ");
				}
			}
			fprintf(poutfile,")):\n");

			fprintf(poutfile,"dFE%id%id0Xfull := subs(booltonum, evalb(",i,j);
			for (k = 1; k <= maxderiv; ++k)  //ranges from 1 to maxderiv
			{
				fprintf(poutfile,"dFE%id%id0%ifull <> 0",i,j,k);
				if (k != maxderiv)
				{
					fprintf(poutfile," or ");
				}
			}
			fprintf(poutfile,")):\n");
			fprintf(poutfile,"\n");
		}

		fprintf(poutfile,"dFE%ix := ",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(poutfile,"dFE%id%idX0full",i,j);
			if ((j+1) != FEn)
			{
				fprintf(poutfile," + ");
			}
		}
		fprintf(poutfile,":\n");
		fprintf(poutfile,"dFE%iy := ",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(poutfile,"dFE%id%id0Xfull",i,j);
			if ((j+1) != FEn)
			{
				fprintf(poutfile," + ");
			}
		}
		fprintf(poutfile,":\n");
		fprintf(poutfile,"\n");
	}

	fprintf(poutfile,"#Combine\n");
	fprintf(poutfile,"dFExfull := ");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"dFE%ix",i);
		if ((i+1) != FEn)
		{
			fprintf(poutfile," + ");
		}
	}
	fprintf(poutfile,";\n");
	fprintf(poutfile,"dFEyfull := ");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"dFE%iy",i);
		if ((i+1) != FEn)
		{
			fprintf(poutfile," + ");
		}
	}
	fprintf(poutfile,";\n");
	fprintf(poutfile,"\n");


	//-------------------- File Initialization --------------------
	fprintf(poutfile,"#:::::::::::::::::::: File Initialization ::::::::::::::::::::\n");
	fprintf(poutfile,"with(CodeGeneration);\n\n");

	fprintf(poutfile,"splitpath := FileTools:-SplitPath(currentdir()):\n\n");

	fprintf(poutfile,"splitFEpath := subsop(numelems(splitpath)-1=\"Funcs\",numelems(splitpath)=\"Defn_FEout.c\",splitpath):\n");
	fprintf(poutfile,"FEfileout := FileTools:-JoinPath(splitFEpath);\n");
	fprintf(poutfile,"if FileTools[Exists](FEfileout) then FileTools:-Remove(FEfileout) end if;\n\n");

	fprintf(poutfile,"splitJSpath := subsop(numelems(splitpath)-1=\"Funcs\",numelems(splitpath)=\"Defn_JSout.c\",splitpath):\n");
	fprintf(poutfile,"JSfileout := FileTools:-JoinPath(splitJSpath);\n");
	fprintf(poutfile,"if FileTools[Exists](JSfileout) then FileTools:-Remove(JSfileout) end if;\n\n");


	//-------------------- Export --------------------
	fprintf(poutfile,"#:::::::::::::::::::: Export ::::::::::::::::::::\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"#----- FE%i -----\n",i);
		fprintf(poutfile,"FE%iout := codegen[makeproc](convert(FE%ifull,float),FEvar):\n",i,i);
		fprintf(poutfile,"C(FE%iout,output=FEfileout,precision=double);\n\n",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(poutfile,"#FE%i du%i\n",i,j);
			for (k = 0; k < FEdlength; ++k)
			{
				fprintf(poutfile,"dFE%id%id%i%iout := codegen[makeproc](convert(dFE%id%id%i%ifull,float),FEvar):\n",i,j,FEdx[k],FEdy[k],i,j,FEdx[k],FEdy[k]);
				fprintf(poutfile,"C(dFE%id%id%i%iout,output=FEfileout,precision=double);\n",i,j,FEdx[k],FEdy[k]);
			}
			fprintf(poutfile,"\n");
		}
	}

	fprintf(poutfile,"#----- FE Structure -----\n");
	fprintf(poutfile,"dFExout := codegen[makeproc](dFExfull):\n");
	fprintf(poutfile,"C(dFExout,output=JSfileout);\n");
	fprintf(poutfile,"dFEyout := codegen[makeproc](dFEyfull):\n");
	fprintf(poutfile,"C(dFEyout,output=JSfileout);\n\n");

	for (i = 0; i < FEn; ++i)
	{
		fprintf(poutfile,"#FE%i\n",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(poutfile,"dFE%id%idX0out := codegen[makeproc](dFE%id%idX0full):\n",i,j,i,j);
			fprintf(poutfile,"C(dFE%id%idX0out,output=JSfileout);\n",i,j);
			fprintf(poutfile,"dFE%id%id0Xout := codegen[makeproc](dFE%id%id0Xfull):\n",i,j,i,j);
			fprintf(poutfile,"C(dFE%id%id0Xout,output=JSfileout);\n",i,j);
			fprintf(poutfile,"\n");
		}
	}

	//Cleanup
	fclose(poutfile);



	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- MapleREAD ICBC Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sprintf(qoutname,"./Gens/MapleREAD_ICBCexport.txt");
	qoutfile = fopen(qoutname,"w");
	if(qoutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file


	//-------------------- Export Variables -------------------------
	fprintf(qoutfile,"\n#:::::::::::::::::::: Export Variables ::::::::::::::::::::\n");
	fprintf(qoutfile,"BCsub1export := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(qoutfile,"diff(u%i(x,y),x) = u%id10, diff(u%i(x,y),y) = u%id01",j,j,j,j);
		if ((j+1) != FEn)
		{
			fprintf(qoutfile,", ");
		}
	}
	fprintf(qoutfile,"};\n");
	fprintf(qoutfile,"BCsub0export := {");
	for (j = 0; j < FEn; ++j)
	{
		fprintf(qoutfile,"u%i(x,y) = u%id00",j,j);
		if ((j+1) != FEn)
		{
			fprintf(qoutfile,", ");
		}
	}
	fprintf(qoutfile,"};\n\n");

	fprintf(qoutfile,"ICvar := [x::numeric,y::numeric,r_H::numeric,M::numeric,alpha::numeric,beta::numeric,delta::numeric]:\n\n");
	fprintf(qoutfile,"BCvar := [x::numeric,y::numeric,r_H::numeric,M::numeric,alpha::numeric,beta::numeric,");
	for (j = 0 ; j < FEn; ++j)
	{
		fprintf(qoutfile,"u%id00::numeric,u%id10::numeric,u%id01::numeric",j,j,j);
		if ((j+1) != FEn)
		{
			fprintf(qoutfile,",");
		}
	}
	fprintf(qoutfile,"]:\n\n");

	//-------------------- Export Equations --------------------
	fprintf(qoutfile,"#:::::::::::::::::::: Export Equations ::::::::::::::::::::\n");
	fprintf(qoutfile,"#----- IC -----\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(qoutfile,"IC%ifull := IC%i:\n",i,i);
	}
	fprintf(qoutfile,"\n");

	fprintf(qoutfile,"#----- BC -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(qoutfile,"#BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			fprintf(qoutfile,"BC%i%c%ifull := subs(BCsub1export, BCsub0export, BC%i%c%i):\n",i,BCdim[bc],BCloc[bc],i,BCdim[bc],BCloc[bc]);
		}
	}
	fprintf(qoutfile,"\n");

	//-------------------- BC deriv calc -------------------------
	fprintf(qoutfile,"#----- dBC -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		for (i = 0; i < FEn; ++i)
		{
			for (j = 0; j < FEn; ++j)
			{
				fprintf(qoutfile,"#BC%i%c%i du%i\n",i,BCdim[bc],BCloc[bc],j);
				for (k = 0; k < BCdlength; ++k)
				{
					fprintf(qoutfile,"dBC%i%c%id%id%i%ifull := simplify(diff(BC%i%c%ifull,u%id%i%i)):\n",i,BCdim[bc],BCloc[bc],j,BCdx[k],BCdy[k],i,BCdim[bc],BCloc[bc],j,BCdx[k],BCdy[k]);
				}
			}
			fprintf(qoutfile,"\n");
		}
	}

	//-------------------- BC structure calc -------------------------
	fprintf(qoutfile,"#----- BC Structure-----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(qoutfile,"#BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			for (j = 0; j < FEn; ++j)
			{
				fprintf(qoutfile,"dBC%i%c%id%idX0full := subs(booltonum, evalb(",i,BCdim[bc],BCloc[bc],j);
				for (k = 0; k <= 1; ++k) //ranges from 0 to 1
				{
					fprintf(qoutfile,"dBC%i%c%id%id%i0full <> 0",i,BCdim[bc],BCloc[bc],j,k);
					if ((k+1) != (1+1))
					{
						fprintf(qoutfile," or ");
					}
				}
				fprintf(qoutfile,")):\n");

				fprintf(qoutfile,"dBC%i%c%id%id0Xfull := subs(booltonum, evalb(",i,BCdim[bc],BCloc[bc],j);
				for (k = 1; k <= 1; ++k) //only a single k = 1
				{
					fprintf(qoutfile,"dBC%i%c%id%id0%ifull <> 0",i,BCdim[bc],BCloc[bc],j,k);
					if ((k+1) != 2)
					{
						fprintf(qoutfile," or ");
					}
				}
				fprintf(qoutfile,")):\n");
				fprintf(qoutfile,"\n");
			}

			fprintf(qoutfile,"dBC%i%c%ix := ",i,BCdim[bc],BCloc[bc]);
			for (j = 0; j < FEn; ++j)
			{
				fprintf(qoutfile,"dBC%i%c%id%idX0full",i,BCdim[bc],BCloc[bc],j);
				if ((j+1) != FEn)
				{
					fprintf(qoutfile," + ");
				}
			}
			fprintf(qoutfile,":\n");
			fprintf(qoutfile,"dBC%i%c%iy := ",i,BCdim[bc],BCloc[bc]);
			for (j = 0; j < FEn; ++j)
			{
				fprintf(qoutfile,"dBC%i%c%id%id0Xfull",i,BCdim[bc],BCloc[bc],j);
				if ((j+1) != FEn)
				{
					fprintf(qoutfile," + ");
				}
			}
			fprintf(qoutfile,":\n");
			fprintf(qoutfile,"\n");
		}

		fprintf(qoutfile,"dBC%c%ix := ",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			fprintf(qoutfile,"dBC%i%c%ix",i,BCdim[bc],BCloc[bc]);
			if ((i+1) != FEn)
			{
				fprintf(qoutfile," + ");
			}
		}
		fprintf(qoutfile,":\n");
		fprintf(qoutfile,"dBC%c%iy := ",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			fprintf(qoutfile,"dBC%i%c%iy",i,BCdim[bc],BCloc[bc]);
			if ((i+1) != FEn)
			{
				fprintf(qoutfile," + ");
			}
		}
		fprintf(qoutfile,":\n");
		fprintf(qoutfile,"\n");
	}

	fprintf(qoutfile,"#Combine\n");
	fprintf(qoutfile,"#BCX\n");
	fprintf(qoutfile,"dBCXxfull := dBCX0x + dBCX1x;\n");
	fprintf(qoutfile,"dBCXyfull := dBCX0y + dBCX1y;\n");
	fprintf(qoutfile,"\n");

	fprintf(qoutfile,"#BCY\n");
	fprintf(qoutfile,"dBCYxfull := dBCY1x + dBCY0x;\n");
	fprintf(qoutfile,"dBCYyfull := dBCY1y + dBCY0y;\n");
	fprintf(qoutfile,"\n");


	//-------------------- File Initialization --------------------
	fprintf(qoutfile,"#:::::::::::::::::::: File Initialization ::::::::::::::::::::\n");
	fprintf(qoutfile,"with(CodeGeneration);\n\n");

	fprintf(qoutfile,"splitpath := FileTools:-SplitPath(currentdir()):\n");
	fprintf(qoutfile,"splitICBCpath := subsop(numelems(splitpath)-1=\"Funcs\",numelems(splitpath)=\"Defn_ICBCout.c\",splitpath):\n");
	fprintf(qoutfile,"ICBCfileout := FileTools:-JoinPath(splitICBCpath);\n\n");
	fprintf(qoutfile,"if FileTools[Exists](ICBCfileout) then FileTools:-Remove(ICBCfileout) end if;\n\n");


	//-------------------- Export --------------------
	fprintf(qoutfile,"#:::::::::::::::::::: Export ::::::::::::::::::::\n");
	fprintf(qoutfile,"#----- IC -----\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(qoutfile,"IC%iout := codegen[makeproc](convert(IC%ifull,float),ICvar):\n",i,i);
		fprintf(qoutfile,"C(IC%iout,output=ICBCfileout,precision=double);\n",i);
	}
	fprintf(qoutfile,"\n");

	fprintf(qoutfile,"#----- BC -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(qoutfile,"#BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			fprintf(qoutfile,"BC%i%c%iout := codegen[makeproc](convert(BC%i%c%ifull,float),BCvar):\n",i,BCdim[bc],BCloc[bc],i,BCdim[bc],BCloc[bc]);
			fprintf(qoutfile,"C(BC%i%c%iout,output=ICBCfileout,precision=double);\n\n",i,BCdim[bc],BCloc[bc]);
			for (j = 0; j < FEn; ++j)
			{
				for (k = 0; k < BCdlength; ++k)
				{
					fprintf(qoutfile,"dBC%i%c%id%id%i%iout := codegen[makeproc](convert(dBC%i%c%id%id%i%ifull,float),BCvar):\n",i,BCdim[bc],BCloc[bc],j,BCdx[k],BCdy[k],i,BCdim[bc],BCloc[bc],j,BCdx[k],BCdy[k]);
					fprintf(qoutfile,"C(dBC%i%c%id%id%i%iout,output=ICBCfileout,precision=double);\n",i,BCdim[bc],BCloc[bc],j,BCdx[k],BCdy[k]);
				}
				fprintf(qoutfile,"\n");
			}
		}
	}

	fprintf(qoutfile,"#----- BC structure -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(qoutfile,"#BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			for (j = 0; j < FEn; ++j)
			{
				fprintf(qoutfile,"dBC%i%c%id%idX0out := codegen[makeproc](dBC%i%c%id%idX0full):\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j);
				fprintf(qoutfile,"C(dBC%i%c%id%idX0out,output=JSfileout);\n",i,BCdim[bc],BCloc[bc],j);

				fprintf(qoutfile,"dBC%i%c%id%id0Xout := codegen[makeproc](dBC%i%c%id%id0Xfull):\n",i,BCdim[bc],BCloc[bc],j,i,BCdim[bc],BCloc[bc],j);
				fprintf(qoutfile,"C(dBC%i%c%id%id0Xout,output=JSfileout);\n",i,BCdim[bc],BCloc[bc],j);
				fprintf(qoutfile,"\n");
			}
		}
	}

	fprintf(qoutfile,"#----- BC X -----\n");
	fprintf(qoutfile,"dBCXxout := codegen[makeproc](dBCXxfull):\n");
	fprintf(qoutfile,"C(dBCXxout,output=JSfileout);\n");
	fprintf(qoutfile,"dBCXyout := codegen[makeproc](dBCXyfull):\n");
	fprintf(qoutfile,"C(dBCXyout,output=JSfileout);\n\n");

	fprintf(qoutfile,"#----- BC Y -----\n");
	fprintf(qoutfile,"dBCYxout := codegen[makeproc](dBCYxfull):\n");
	fprintf(qoutfile,"C(dBCYxout,output=JSfileout);\n");
	fprintf(qoutfile,"dBCYyout := codegen[makeproc](dBCYyfull):\n");
	fprintf(qoutfile,"C(dBCYyout,output=JSfileout);\n\n");


	//Cleanup
	fclose(qoutfile);



	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- Jacobian Structure Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sprintf(JSoutname,"./Funcs/Decl_JSFuncs.c");
	JSoutfile = fopen(JSoutname,"w");
	if(JSoutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file

	fprintf(JSoutfile,"//-------------------- C Import from FE and ICBC --------------------\n\n");


	//----- Jacobian Structure Declaration -----
	fprintf(JSoutfile,"//----- Jacobian Structure Declaration -----\n");
	fprintf(JSoutfile,"int dFExout ();\n");
	fprintf(JSoutfile,"int dFEyout ();\n");
	fprintf(JSoutfile,"\n");
	fprintf(JSoutfile,"int dBCXxout ();\n");
	fprintf(JSoutfile,"int dBCXyout ();\n");
	fprintf(JSoutfile,"\n");
	fprintf(JSoutfile,"int dBCYxout ();\n");
	fprintf(JSoutfile,"int dBCYyout ();\n");
	fprintf(JSoutfile,"\n");

	for (i = 0; i < FEn; ++i)
	{
		fprintf(JSoutfile,"//FE%i\n",i);
		for (j = 0; j < FEn; ++j)
		{
			fprintf(JSoutfile,"int dFE%id%idX0out ();\n",i,j);
			fprintf(JSoutfile,"int dFE%id%id0Xout ();\n",i,j);
			fprintf(JSoutfile,"\n");
		}
	}

	for (bc = 0; bc < BClength; ++bc)
	{
		fprintf(JSoutfile,"//BC %c%i\n",BCdim[bc],BCloc[bc]);
		for (i = 0; i < FEn; ++i)
		{
			for (j = 0; j < FEn; ++j)
			{
				fprintf(JSoutfile,"int dBC%i%c%id%idX0out ();\n",i,BCdim[bc],BCloc[bc],j);
				fprintf(JSoutfile,"int dBC%i%c%id%id0Xout ();\n",i,BCdim[bc],BCloc[bc],j);
				fprintf(JSoutfile,"\n");
			}
		}
	}


	//Cleanup
	fclose(JSoutfile);



	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//-------------------- C Funcs Header Declaration Export --------------------
	//-------------------------------------------------------------------------------------------------------------------------------------------
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sprintf(coutname,"./Funcs/Decl_CFuncs.c");
	coutfile = fopen(coutname,"w");
	if(coutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file

	fprintf(coutfile,"//-------------------- C Import from FE and ICBC --------------------\n\n");


	//----- NPcalc Declaration -----
	fprintf(coutfile,"//----- Calculates Newton polynomial coefficients -----\n");
	fprintf(coutfile,"int BVP_NPcalc(const int p, const int r, const int rn, const int I, int s[], double x[], const int dim, const int IJ, double y[], const int pn,");
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(coutfile," double Q%i_a[],",k);
	}
	for (k = 0; k <= maxderiv; ++k)
	{
		fprintf(coutfile," double *d%iydk,",k);
	}
	fprintf(coutfile," const double xk);\n");
	fprintf(coutfile,"\n");


	//----- FE Declaration -----
	fprintf(coutfile,"//----- FE Declaration -----\n");

	//Calculate declaration variable
	strcpy(vardecl,"(double xk, double yk, double r_H, double M, double alpha, double beta,");
	for (j = 0 ; j < FEn; ++j)
	{
		for (k = 0; k < FEdlength; ++k)
		{
			sprintf(vardecltemp," double u%id%i%ik",j,FEdx[k],FEdy[k]);
			if ((k+1) != FEdlength)
			{
				strcat(vardecltemp,",");
			}
			strcat(vardecl, vardecltemp);
		}
		if ((j+1) != FEn)
		{
			strcat(vardecl,",");
		}
	}
	strcat(vardecl,")");

	for (i = 0; i < FEn; ++i)
	{
		fprintf(coutfile,"//FE%i\n",i);
		fprintf(coutfile,"double FE%iout %s;\n",i,vardecl);

		for (j = 0; j < FEn; ++j)
		{
			for (k = 0; k < FEdlength; ++k)
			{
				fprintf(coutfile,"double dFE%id%id%i%iout %s;\n",i,j,FEdx[k],FEdy[k],vardecl);
			}
			fprintf(coutfile,"\n");
		}
	}
	fprintf(coutfile,"\n");


	//----- IC Declaration -----
	fprintf(coutfile,"//----- IC Declaration -----\n");
	for (i = 0; i < FEn; ++i)
	{
		fprintf(coutfile,"double IC%iout (double xk, double yk, double r_H, double M, double alpha, double beta, double delta);\n",i);
	}
	fprintf(coutfile,"\n\n");


	//----- BC Declaration -----
	fprintf(coutfile,"//----- BC Declaration -----\n");
	for (bc = 0; bc < BClength; ++bc)
	{
		for (i = 0; i < FEn; ++i)
		{
			fprintf(coutfile,"//BC%i%c%i\n",i,BCdim[bc],BCloc[bc]);
			fprintf(coutfile,"double BC%i%c%iout (double xk, double yk, double r_H, double M, double alpha, double beta,",i,BCdim[bc],BCloc[bc]);
			for (j = 0 ; j < FEn; ++j)
			{
				fprintf(coutfile," double u%id00k, double u%id10k, double u%id01k",j,j,j);
				if ((j+1) != FEn)
				{
					fprintf(coutfile,",");
				}
			}
			fprintf(coutfile,");\n\n");
			for (j = 0; j < FEn; ++j)
			{
				for (k = 0; k < BCdlength; ++k)
				{
					fprintf(coutfile,"double dBC%i%c%id%id%i%iout (double xk, double yk, double r_H, double M, double alpha, double beta,",i,BCdim[bc],BCloc[bc],j,BCdx[k],BCdy[k]);
					for (l = 0 ; l < FEn; ++l)
					{
						fprintf(coutfile," double u%id00k, double u%id10k, double u%id01k",l,l,l);
						if ((l+1) != FEn)
						{
							fprintf(coutfile,",");
						}
					}
					fprintf(coutfile,");\n");
				}
				fprintf(coutfile,"\n");
			}
		}
	}

	//Cleanup
	fclose(coutfile);



	free(FEdx);
	free(FEdy);

	return 0;
}
