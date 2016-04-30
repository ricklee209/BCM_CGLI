



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 

extern int NBC;

extern int NBC_plus;

void BCM_Reading_IBM
(
// =================================================== //
int myid,

int *NNBC,

int mp_switch,

double (*weight) = new double[NBC*8+1],

int (*GCindex) = new int[NBC*4+1],

int (*IPsur) = new int[NBC*4+1],

double (*Nor_D) = new double[NBC+1],

double (*Nvec) = new double[NBC*3+1]
// =================================================== //
)

{

#include "BCM.h"
#include "prm.h"

	
	int iNBC, Ntemp;

	char str[1024];

	

	FILE *fptr;
	char file_name[100];

	if(mp_switch == 1) sprintf(file_name,"IBM_plus""%0.5d"".dat",myid);  
	if(mp_switch == -1) sprintf(file_name,"IBM_minus""%0.5d"".dat",myid); 

	fptr = fopen(file_name,"r"); 

	for (iNBC = 1; iNBC <= *NNBC; iNBC++) {

		
		Ntemp = (iNBC-1)*4;

		fscanf(fptr,"%d\t%d\t%d\t%d\n",&GCindex[Ntemp+1],&GCindex[Ntemp+2],&GCindex[Ntemp+3],&GCindex[Ntemp+4]);

		Ntemp = (iNBC-1)*4;

		fscanf(fptr,"%d\t%d\t%d\t%d\n",&IPsur[Ntemp+1],&IPsur[Ntemp+2],&IPsur[Ntemp+3],&IPsur[Ntemp+4]);

		Ntemp = (iNBC-1)*3;

		fscanf(fptr,"%lf\t%lf\t%lf\n",&Nvec[Ntemp+1],&Nvec[Ntemp+2],&Nvec[Ntemp+3]);

		Ntemp = (iNBC-1)*8;

		fscanf(fptr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			   &Nor_D[iNBC],
			   &weight[Ntemp+1],&weight[Ntemp+2],&weight[Ntemp+3],&weight[Ntemp+4],
			   &weight[Ntemp+5],&weight[Ntemp+6],&weight[Ntemp+7],&weight[Ntemp+8]);
		
		//if (myid == 0 && mp_switch==-1) printf("%f\t%f\t%f\t%f\n",Nor_D[iNBC],weight[Ntemp+2],weight[Ntemp+3],weight[Ntemp+4]);

	}

	fclose(fptr);	

}