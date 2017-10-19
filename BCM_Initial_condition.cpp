



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

	int switch_initial,

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
                    
                    
                    U1q[icube][i][j][k][0] = U1[icube][i][j][k][0];
					U1q[icube][i][j][k][1] = U1[icube][i][j][k][1];
					U1q[icube][i][j][k][2] = U1[icube][i][j][k][2];
					U1q[icube][i][j][k][3] = U1[icube][i][j][k][3];
					U1q[icube][i][j][k][4] = U1[icube][i][j][k][4];
                    
                    

				}
			}
		}
	}

	

	// --------------------------------------------------------------------------------------- //
	// --------------------------------------------------------------------------------------- //		

	if (switch_initial == 1) {

		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];
		MPI_Offset disp;

		MPI_Offset nx_out = NcubeX+2;
		MPI_Offset ny_out = NcubeY+2;
		MPI_Offset nz_out = NcubeZ+2;

		MPI_Offset int_size = sizeof(int);

		MPI_Offset float_size = sizeof(float);

		float (*Solution) = new float[(ncube-1)*nx_out*ny_out*nz_out*5+4*(ncube-1)];

		MPI_Offset x_gcount[np], x_gdisp[np], q_gcount[np], q_gdisp[np];

		istart = 0;
		x_gdisp[np-1] = 0; 
		q_gdisp[np-1] = 0; 

		for (i = 0; i < np-1; i++) { 

            for (icube = 0; icube < MPI_Ncube; icube++) {

                if (rank_map[0][icube] == i) {

                    icount = icube;
                    break;
                    
                }

            }

            x_gdisp[i] = (MPI_Offset)icount*(MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*3;
            q_gdisp[i] = (MPI_Offset)icount*(MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*5+(MPI_Offset)icount*4;

            
        }



		char data[100];

		sprintf(data,"initial.q");

		MPI_File fh_initial;

		MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_initial ) ; 




		disp =  ((MPI_Offset)MPI_Ncube*3)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+q_gdisp[myid]*(MPI_Offset)sizeof(float);



		//MPI_File_set_view(fh_initial, disp, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);

		//MPI_File_read(fh_initial, Solution, (ncube-1)*nx_out*ny_out*nz_out*5+(ncube-1)*4, MPI_FLOAT, MPI_STATUS_IGNORE);

		MPI_File_read_at_all(fh_initial, disp, Solution, (ncube-1)*nx_out*ny_out*nz_out*5+(ncube-1)*4, MPI_FLOAT, MPI_STATUS_IGNORE);


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

	// --------------------------------------------------------------------------------------- //
	// --------------------------------------------------------------------------------------- //



		
		sprintf(data,"Qinitial.q");

		MPI_File fh_previous;

		MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_previous ) ; 

		disp =  ((MPI_Offset)MPI_Ncube*3)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+q_gdisp[myid]*(MPI_Offset)sizeof(float);

		MPI_File_read_at_all(fh_previous, disp, Solution, (ncube-1)*nx_out*ny_out*nz_out*5+(ncube-1)*4, MPI_FLOAT, MPI_STATUS_IGNORE);


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
						U1q[icube][k][j][i][0] = Solution[icount];

					}
				}
			}

			for (i = n_buffer-1; i <= nx+1; i++) {
				for (j = n_buffer-1; j <= ny+1; j++) { 
					for (k = n_buffer-1; k <= nz+1; k++) { 

						icount = icount + 1;
						U1q[icube][k][j][i][1] = Solution[icount];

					}
				}
			}

			for (i = n_buffer-1; i <= nx+1; i++) {
				for (j = n_buffer-1; j <= ny+1; j++) { 
					for (k = n_buffer-1; k <= nz+1; k++) { 

						icount = icount + 1;
						U1q[icube][k][j][i][2] = Solution[icount];

					}
				}
			}

			for (i = n_buffer-1; i <= nx+1; i++) {
				for (j = n_buffer-1; j <= ny+1; j++) { 
					for (k = n_buffer-1; k <= nz+1; k++) { 

						icount = icount + 1;
						U1q[icube][k][j][i][3] = Solution[icount];

					}
				}
			}

			for (i = n_buffer-1; i <= nx+1; i++) {
				for (j = n_buffer-1; j <= ny+1; j++) { 
					for (k = n_buffer-1; k <= nz+1; k++) { 

						icount = icount + 1;
						U1q[icube][k][j][i][4] = Solution[icount];

					}
				}
			}

		}

		MPI_File_close( &fh_previous );


	
	// ----------------------------------------------------------------- //
	// ------------------------- Previous_step ------------------------- //







	

	// ------------------------- Previous_step ------------------------- //
	// ----------------------------------------------------------------- //

		
	
		
		delete []Solution;


	
	}    // ---- if (switch_initial == 1) ---- //

	
	for (icube = 1; icube < ncube; icube++) {  

#pragma omp parallel for private(j,k)

		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {
				for (k = 0; k <= nzzz; k++) {  

					U1_[icube][i][j][k][0] = U1[icube][i][j][k][0]; 
					U1_[icube][i][j][k][1] = U1[icube][i][j][k][1]; 
					U1_[icube][i][j][k][2] = U1[icube][i][j][k][2]; 
					U1_[icube][i][j][k][3] = U1[icube][i][j][k][3]; 
					U1_[icube][i][j][k][4] = U1[icube][i][j][k][4]; 

					
				}
			}
		}

	}



}