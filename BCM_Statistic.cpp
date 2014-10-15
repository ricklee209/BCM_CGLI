




#include <mpi.h>
#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;      

void BCM_Statistic
(
// =================================================== //
int myid,
int ncube,

int step,
int start_step,
int statistic_step,
int dp_step,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

double (*U1)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size], 
double (*U2)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size], 
double (*U3)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size], 
double (*U4)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size], 
double (*U5)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],

double (*Pall)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*VVall)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]

// =================================================== //
)

{

	
#include "BCM.h"
#include "prm.h"
#include "MPI_prm.h"

	double (*Pm)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*VVm)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*Pp)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];


	double rho,U,V,W,VV,P,C,T,H;

	int Average_step = statistic_step;

	#pragma omp parallel for private(i,j,k,rho,U,V,W,VV,P)	
	for (icube = 1; icube < ncube; icube++) { 
		for (i = n_buffer-1; i <= nxx; i++) {
			for (j = n_buffer-1; j <= nyy; j++) {
				for (k = n_buffer-1; k <= nzz; k++) {  

					/* flux parameter */
					rho = U1[icube][i][j][k];
					U = U2[icube][i][j][k]/rho;
					V = U3[icube][i][j][k]/rho;
					W = U4[icube][i][j][k]/rho;     
					VV = U*V+V*V+W*W;
					P = (U5[icube][i][j][k]-0.5*rho*VV)*(K-1)-101300;
					
					Pall[icube][i][j][k] = Pall[icube][i][j][k]+P;

					VVall[icube][i][j][k] = VVall[icube][i][j][k]+VV;

				}
			}
		}
	}


	if ( (step%Average_step) == 0) {

#pragma omp parallel for private(i,j,k)	
		for (icube = 1; icube < ncube; icube++) { 
			for (i = n_buffer-1; i <= nxx; i++) {
				for (j = n_buffer-1; j <= nyy; j++) {
					for (k = n_buffer-1; k <= nzz; k++) {  

						Pm[icube][i][j][k] = Pall[icube][i][j][k]/Average_step;

						VVm[icube][i][j][k] = VVall[icube][i][j][k]/Average_step;

					}
				}
			}
		}



#pragma omp parallel for private(i,j,k)	
		for (icube = 1; icube < ncube; icube++) { 
			for (i = n_buffer-1; i <= nxx; i++) {
				for (j = n_buffer-1; j <= nyy; j++) {
					for (k = n_buffer-1; k <= nzz; k++) {  

						Pall[icube][i][j][k] = 0;

						VVall[icube][i][j][k] = 0;

					}
				}
			}
		}


	}    // ---- if ( (step%statistic_step) == 0) ---- //




	
// =========================== Output averaged result =========================== //

	if ( step > (Average_step+start_step+2) && step%dp_step == 0) {

		double (*dPout)[X_size][Y_size][Z_size] = new double[MPI_Ncube][X_size][Y_size][Z_size];
		double (*Pmout)[X_size][Y_size][Z_size] = new double[MPI_Ncube][X_size][Y_size][Z_size];
		double (*VVmout)[X_size][Y_size][Z_size] = new double[MPI_Ncube][X_size][Y_size][Z_size];
	
#pragma omp parallel for private(i,j,k,rho,U,V,W,VV,P)	
		for (icube = 1; icube < ncube; icube++) { 
			for (i = n_buffer-1; i <= nxx; i++) {
				for (j = n_buffer-1; j <= nyy; j++) {
					for (k = n_buffer-1; k <= nzz; k++) {  

						rho = U1[icube][i][j][k];
						U = U2[icube][i][j][k]/rho;
						V = U3[icube][i][j][k]/rho;
						W = U4[icube][i][j][k]/rho;     
						VV = U*V+V*V+W*W;
						P = (U5[icube][i][j][k]-0.5*rho*VV)*(K-1)-101300;

						Pp[icube][i][j][k] = P-Pm[icube][i][j][k];

					}
				}
			}
		}

		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];
		MPI_Request requ1, reqps1;

		int nx_out = NcubeX+2;
		int ny_out = NcubeY+2;
		int nz_out = NcubeZ+2;

		int ncube_out = MPI_Ncube;

		int x_gcount[np], x_gdisp[np], q_gcount[np], q_gdisp[np];

		istart = 0;
		x_gdisp[0] = 0; 
		q_gdisp[0] = 0; 

		for (i = 0; i < np; i++) { 

		icount = 0;
		for (icube = 0; icube < MPI_Ncube; icube++) {

			if (rank_map[0][icube] == i) {

				icount = icount+1;
				istart = rank_map[1][icube];

			}
		}

		x_gcount[i] = icount*X_size;
		q_gcount[i] = icount*X_size*Y_size*Z_size;

		if (i < (np-1)) {
			
			x_gdisp[i+1] = x_gdisp[i]+istart*X_size;
			q_gdisp[i+1] = q_gdisp[i]+istart*X_size*Y_size*Z_size;
			}

		}
	
	
		idest = 0;

		icount = q_gcount[myid];

		MPI_Gatherv((void *)&Pp[1][0][0][0], icount, MPI_DOUBLE, (void *)&dPout[0][0][0][0], q_gcount, q_gdisp, MPI_DOUBLE, idest, comm);
		
		MPI_Gatherv((void *)&Pm[1][0][0][0], icount, MPI_DOUBLE, (void *)&Pmout[0][0][0][0], q_gcount, q_gdisp, MPI_DOUBLE, idest, comm);

		MPI_Gatherv((void *)&VVm[1][0][0][0], icount, MPI_DOUBLE, (void *)&VVmout[0][0][0][0], q_gcount, q_gdisp, MPI_DOUBLE, idest, comm);

		if (myid == idest) {

			char data[100];
			FILE *fptr;
			double temp = 1.0;
			int itemp = 1;

			sprintf(data,"qdP""%0.5d_%0.5d"".q",step, idest);
			fptr = fopen(data,"wb");

			fwrite(&ncube_out, sizeof(int), 1,fptr);

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				fwrite(&nx_out, sizeof(int), 1,fptr);
				fwrite(&ny_out, sizeof(int), 1,fptr);
				fwrite(&nz_out, sizeof(int), 1,fptr);
				fwrite(&itemp, sizeof(int), 1,fptr);

			}


			for (icube = 0; icube < MPI_Ncube; icube++)  {
				for (i = n_buffer-1; i <= nx+1; i++) {
					for (j = n_buffer-1; j <= ny+1; j++) { 
						for (k = n_buffer-1; k <= nz+1; k++) { 

							fwrite(&dPout[icube][k][j][i],sizeof(double),1,fptr);

						}
					}
				}

			}





			sprintf(data,"qPm""%0.5d_%0.5d"".q",step, idest);
			fptr = fopen(data,"wb");

			fwrite(&ncube_out, sizeof(int), 1,fptr);

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				fwrite(&nx_out, sizeof(int), 1,fptr);
				fwrite(&ny_out, sizeof(int), 1,fptr);
				fwrite(&nz_out, sizeof(int), 1,fptr);
				fwrite(&itemp, sizeof(int), 1,fptr);

			}

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				for (i = n_buffer-1; i <= nx+1; i++) {
					for (j = n_buffer-1; j <= ny+1; j++) { 
						for (k = n_buffer-1; k <= nz+1; k++) { 

							fwrite(&Pmout[icube][k][j][i],sizeof(double),1,fptr);

						}
					}
				}

			}



			sprintf(data,"qVVm""%0.5d_%0.5d"".q",step, idest);
			fptr = fopen(data,"wb");

			fwrite(&ncube_out, sizeof(int), 1,fptr);

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				fwrite(&nx_out, sizeof(int), 1,fptr);
				fwrite(&ny_out, sizeof(int), 1,fptr);
				fwrite(&nz_out, sizeof(int), 1,fptr);
				fwrite(&itemp, sizeof(int), 1,fptr);
				
			}

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				for (i = n_buffer-1; i <= nx+1; i++) {
					for (j = n_buffer-1; j <= ny+1; j++) { 
						for (k = n_buffer-1; k <= nz+1; k++) { 

							fwrite(&VVmout[icube][k][j][i],sizeof(double),1,fptr);

						}
					}
				}

			}

			
			fclose(fptr);

		}    // ---- if (myid == idest) ---- //

		delete [] dPout;
		delete [] Pmout;
		delete [] VVmout;

	}    // ---- if ( step%dp_step == 0) ---- //

// =========================== Output averaged result =========================== //
	
	delete [] Pm;
	delete [] VVm;
	delete [] Pp;
	

}
