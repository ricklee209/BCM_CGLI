



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

double (*Dweight) = new double[NBC*8+1],

double (*Nweight) = new double[NBC*8+1],

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
			   &Dweight[Ntemp+1],&Dweight[Ntemp+2],&Dweight[Ntemp+3],&Dweight[Ntemp+4],
			   &Dweight[Ntemp+5],&Dweight[Ntemp+6],&Dweight[Ntemp+7],&Dweight[Ntemp+8]);

		fscanf(fptr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			   &Nweight[Ntemp+1],&Nweight[Ntemp+2],&Nweight[Ntemp+3],&Nweight[Ntemp+4],
			   &Nweight[Ntemp+5],&Nweight[Ntemp+6],&Nweight[Ntemp+7],&Nweight[Ntemp+8]);
		
		//if (myid == 0 && mp_switch==1) printf("%f\t%f\t%f\t%f\n",Nweight[Ntemp+1],Nweight[Ntemp+2],Nweight[Ntemp+3],Nweight[Ntemp+4]);

	}

	fclose(fptr);	

}