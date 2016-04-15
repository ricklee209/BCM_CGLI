



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

void BCM_Mean_pressure_coefficient_Sphere 
(
// ================================================================================ //
int myid,
int ncube,

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



double Raidus = 0.005;


char file_name[100];

int icount, g_icount;

double rho, U, V, W, VV, P, T, temp;
double dis,xc,yc,zc;
double Cpre, Cpre_mean, g_Cpre_mean;





FILE *fptr;

sprintf(file_name,"Cp""%0.5d"".dat",myid);    
fptr = fopen(file_name,"w"); 



	icount = 0;
	for (icube = 1; icube < ncube; icube++) {  

		for (i = n_buffer; i < nxx; i++) {
			for (j = n_buffer; j < nyy; j++) {
				for (k = n_buffer; k < nzz; k++) {  

					xc = Xcnt[icube][i];
					yc = Ycnt[icube][j];
					zc = Zcnt[icube][k];

					dis = sqrt(xc*xc+yc*yc+zc*zc);


					if(FWS[icube][i][j][k] == IGHOST && (dis - Raidus) > 0.0 ) icount = icount + 1;


				}
			}
		}

	}    // ---- for (icube = 1; icube < ncube; icube++) ---- //

	

	Cpre_mean = 0.0;

//#pragma omp parallel for private(i,j,k,xc,yc,zc,dis,rho,U,V,W,P,T,Nu)

	for (icube = 1; icube < ncube; icube++) {  

		for (i = n_buffer; i < nxx; i++) {
			for (j = n_buffer; j < nyy; j++) {
				for (k = n_buffer; k < nzz; k++) {  

					if(FWS[icube][i][j][k] == IGHOST && csl[icube] == 0) {

						xc = Xcnt[icube][i];
						yc = Ycnt[icube][j];
						zc = Zcnt[icube][k];

						dis = sqrt(xc*xc+yc*yc+zc*zc);

						if ( (dis - Raidus) > 0.0 && yc > 0 ) {

							rho = U1_[icube][i][j][k][0];
							U = U1_[icube][i][j][k][1]/rho;
							V = U1_[icube][i][j][k][2]/rho;
							W = U1_[icube][i][j][k][3]/rho;
							P = ( U1_[icube][i][j][k][4]-0.5*rho*(U*U+V*V+W*W) )*(K-1);
							T = P/R/rho;

							Cpre = (P - P0) / (0.5*rho0*U0*U0);

							Cpre_mean = Cpre_mean + Cpre;

							temp = acos(-xc/dis)*180.0/3.1415926;

							//fprintf(fptr,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",xc,yc,zc,temp,Nu,T-Th,rho);

							fprintf(fptr,"%f\t%f\n",temp,Cpre);

						}    // ---- if (xc*xc+yc*yc+zc*zc > 1.0) ---- //

					}    // ---- if(FWS[icube][i][j][k] == IGHOST) ---- //


				}
			}
		}

	}    // ---- for (icube = 1; icube < ncube; icube++) ---- //

	
	fclose(fptr);



}