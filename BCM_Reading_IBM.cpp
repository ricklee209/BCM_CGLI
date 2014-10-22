



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

double (*N_dis) = new double[NBC*4+1],

int (*GCindex) = new int[NBC*4+1],

int (*IPsur) = new int[NBC*4+1]

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

		fscanf(fptr,"%d\t%d\t%d\t%d\t",&IPsur[Ntemp+1],&IPsur[Ntemp+2],&IPsur[Ntemp+3],&IPsur[Ntemp+4]);

		Ntemp = (iNBC-1)*8;

		fscanf(fptr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			   &weight[Ntemp+1],&weight[Ntemp+2],&weight[Ntemp+3],&weight[Ntemp+4],
			   &weight[Ntemp+5],&weight[Ntemp+6],&weight[Ntemp+7],&weight[Ntemp+8]);

		Ntemp = (iNBC-1)*4;

		fscanf(fptr,"%lf\t%lf\t%lf\t%lf\n",&N_dis[Ntemp+1],&N_dis[Ntemp+2],&N_dis[Ntemp+3],&N_dis[Ntemp+4]);

		
		//if (myid == 0 && mp_switch==-1) printf("%f\t%f\t%f\t%f\n",weight[Ntemp+1],weight[Ntemp+2],weight[Ntemp+3],weight[Ntemp+4]);

	}

	fclose(fptr);	

}