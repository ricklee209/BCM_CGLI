



#include <mpi.h>
#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>

#include "Resolution.h"

extern int Ncube;    
extern int N_wallcube;    

void BCM_Initial_condition
(
// ================================================================================ //
int myid,
int ncube,
int n_wallcube,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
)
{
#include "BCM.h"
#include "MPI_prm.h"
#include "prm.h"

	double rho,U,V,W,VV,P,C,T,h,H;
	
	int iblank;

// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //	


	for (icube = 1; icube < ncube; icube++) { 
#pragma omp parallel for private(j,k,rho,U,V,W,VV,P,iblank)
		
		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {
				for (k = 0; k <= nzzz; k++) {  
				
						
					iblank = FWS[icube][i][j][k];

					if (iblank > ISOLID ) {

						rho = rho0;
						U = U0;
						V = V0;
						W = W0;
						VV = U*U+V*V+W*W;
						P = P0;

						U1[icube][i][j][k][0] = rho;
						U1[icube][i][j][k][1] = rho*U;
						U1[icube][i][j][k][2] = rho*V;
						U1[icube][i][j][k][3] = rho*W;
						U1[icube][i][j][k][4] = P/(K-1)+0.5*rho*VV;
						
					}
						
					else {
						
						U1[icube][i][j][k][0] = rho0;
						U1[icube][i][j][k][1] = 0;
						U1[icube][i][j][k][2] = 0;
						U1[icube][i][j][k][3] = 0;
						U1[icube][i][j][k][4] = P0/(K-1);

					}
						
				}
			}
		}
	}
	

// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //		
	/*
	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Status istat[8];
	MPI_Offset disp;
	
	int nx_out = NcubeX+2;
	int ny_out = NcubeY+2;
	int nz_out = NcubeZ+2;

	double (*Solution) = new double[(ncube-1)*nx_out*ny_out*nz_out*5+4*(ncube-1)];

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
		q_gcount[i] = icount*X_size*Y_size*Z_size*Ndim;

		if (i < (np-1)) {
			
			x_gdisp[i+1] = x_gdisp[i]+istart*nx_out*ny_out*nz_out*3;
			q_gdisp[i+1] = q_gdisp[i]+istart*nx_out*ny_out*nz_out*5+istart*4;

		}

	}
	

	
	char data[100];

	sprintf(data,"initial.q");
	
	MPI_File fh_initial;

	MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_initial ) ; 

	
	disp = (MPI_Ncube*3+1)*sizeof(int)+q_gdisp[myid]*sizeof(double);

	MPI_File_set_view(fh_initial, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	MPI_File_read(fh_initial, Solution, (ncube-1)*nx_out*ny_out*nz_out*5+(ncube-1)*4, MPI_DOUBLE, MPI_STATUS_IGNORE);

	
	icount = -1;
	for (icube = 1; icube < ncube; icube++) {  

		icount = icount + 1;
		Solution[icount] = 1.0;
		
		icount = icount + 1;
		Solution[icount] = 1.0;
		
		icount = icount + 1;
		Solution[icount] = 1.0;
		
		icount = icount + 1;
		Solution[icount] = 1.0;

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					U1[icube][k][j][i][0] = Solution[icount];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					U1[icube][k][j][i][1] = Solution[icount];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					U1[icube][k][j][i][2] = Solution[icount];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					U1[icube][k][j][i][3] = Solution[icount];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					U1[icube][k][j][i][4] = Solution[icount];

				}
			}
		}
	
	}
	
	MPI_File_close( &fh_initial );

	delete []Solution;
	*/
	
// --------------------------------------------------------------------------------------- //
// --------------------------------------------------------------------------------------- //		
	
}