



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


#include "Resolution.h"

extern int Ncube;    
extern int N_wallcube;    

void BCM_Nusselt_Sphere 
(
// ================================================================================ //
int myid,
int ncube,

double Th,

int (*csl) = new int[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],


double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// ================================================================================ //
)

{

	
#include "BCM.h"
#include "prm.h"
	
MPI_Comm comm;
comm=MPI_COMM_WORLD;



double Ra = 1e9*Pr_L, Raidus = 0.2;


char file_name[100];

int icount, g_icount, iCp;

double rho, U, V, W, VV, P, T, temp;
double dis,xc,yc,zc;
double Cpre, Cpre_mean, g_Cpre_mean;
double Nu, Nu_mean, g_Nu_mean;


double (*Cp_suf)[181] = new double[3][181];
double (*gCp)[181] = new double[3][181];


#pragma omp parallel
for (iCp = 1; iCp <= 180; iCp++) {  

	Cp_suf[0][iCp] = 0.0;
	Cp_suf[1][iCp] = 0.0;
	Cp_suf[3][iCp] = 0.0;
	
	gCp[0][iCp] = 0.0;
	gCp[1][iCp] = 0.0;
	gCp[3][iCp] = 0.0;

	}   

	
	Cpre_mean = 0.0;
	icount = 0;


	for (icube = 1; icube < ncube; icube++) {  

		for (i = n_buffer; i < nxx; i++) {
			for (j = n_buffer; j < nyy; j++) {
				for (k = n_buffer; k < nzz; k++) {  

					if(FWS[icube][i][j][k] == IGHOST && csl[icube] == 0) {

						xc = Xcnt[icube][i];
						yc = Ycnt[icube][j];
						zc = Zcnt[icube][k];

						dis = sqrt(xc*xc+yc*yc+zc*zc);

						if ( (dis - Raidus) > 0.0 ) {

							rho = U1_[icube][i][j][k][0];
							U = U1_[icube][i][j][k][1]/rho;
							V = U1_[icube][i][j][k][2]/rho;
							W = U1_[icube][i][j][k][3]/rho;
							P = ( U1_[icube][i][j][k][4]-0.5*rho*(U*U+V*V+W*W) )*(K-1);
							T = P/R/rho;

							Cpre = (P - P0) / (0.5*rho0*U0*U0);

							Nu = 2*Raidus/(Th-T0)*(T-Th)/(dis-Raidus)/sqrt(sqrt(Ra));  // ---- Nu/Ra^(1/4)


							temp = acos(-xc/dis)*180.0/3.1415926;


							for (iCp = 1; iCp <= 180; iCp++) {  

								if (temp >= iCp*1.0 && temp <= iCp+1.0 ) {

									Cp_suf[2][iCp] = Cp_suf[2][iCp] + 1.0;
									Cp_suf[1][iCp] = Cp_suf[1][iCp] + Nu;
									Cp_suf[0][iCp] = Cp_suf[0][iCp] + Cpre;

								}
							}

						}    // ---- if (xc*xc+yc*yc+zc*zc > 1.0) ---- //

					}    // ---- if(FWS[icube][i][j][k] == IGHOST) ---- //


				}
			}
		}

	}    // ---- for (icube = 1; icube < ncube; icube++) ---- //



	
	MPI_Reduce ((void*)Cp_suf,(void*)gCp, 543, MPI_DOUBLE, MPI_SUM, 0, comm );

	
	if(myid == 0) {

		FILE *fptr;

		fptr = fopen("Cp.dat","w"); 

		for (iCp = 1; iCp <= 179; iCp++) { 

			gCp[0][iCp] = gCp[0][iCp]/(gCp[1][iCp]+0.0000001);
			gCp[1][iCp] = gCp[1][iCp]/(gCp[2][iCp]+0.0000001);
			
			if( fabs(gCp[0][iCp])>0.0000001 ) fprintf(fptr,"%f\t%f\t%f\n",iCp*1.0-0.5,gCp[0][iCp],gCp[1][iCp]);
			
		}

		fclose(fptr);

	}

	
	delete []Cp_suf;
	delete []gCp;




}