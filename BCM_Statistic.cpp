




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

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U2_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U3_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U4_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]
// =================================================== //
)

{


// ====== Multi-zone structured scalar function file format ====== //

	
#include "BCM.h"
#include "prm.h"
#include "MPI_prm.h"

	double (*U_mean)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*T_mean)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*U_int)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*T_int)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];


	double rho,U,V,W,VV,P,C,T,H;

	int Average_step = statistic_step;

	int num_variable_output = 4; // output four values //
    
    // U_mean, T_mean, U_intensity, T_intensity //

	int itemp = 0;
    
    itemp = step-start_step+1.0;    // ==== averaged steps ==== //

    
	#pragma omp parallel for private(i,j,k,rho,U,V,W,VV,P,T)	
	for (icube = 1; icube < ncube; icube++) { 
		for (i = n_buffer-1; i <= nxx; i++) {
			for (j = n_buffer-1; j <= nyy; j++) {
				for (k = n_buffer-1; k <= nzz; k++) {  

					/* flux parameter */
					rho = U1[icube][i][j][k][0];
					U = U1[icube][i][j][k][1]/rho;
					
					U1_sum[icube][i][j][k] = U1_sum[icube][i][j][k]+U;
                    U2_sum[icube][i][j][k] = U2_sum[icube][i][j][k]+U*U;
                    
                    

				}
			}
		}
	}

#pragma omp parallel for private(i,j,k)	
	for (icube = 1; icube < ncube; icube++) { 
		for (i = n_buffer-1; i <= nxx; i++) {
			for (j = n_buffer-1; j <= nyy; j++) {
				for (k = n_buffer-1; k <= nzz; k++) {  

					U_mean[icube][i][j][k] = U1_sum[icube][i][j][k]/itemp;
                    
					U_int[icube][i][j][k] = sqrt(fabs(U2_sum[icube][i][j][k]/itemp - U_mean[icube][i][j][k]*U_mean[icube][i][j][k]));

				}
			}
		}
	}
    
    
    #pragma omp parallel for private(i,j,k,rho,U,V,W,VV,P,T)	
	for (icube = 1; icube < ncube; icube++) { 
		for (i = n_buffer-1; i <= nxx; i++) {
			for (j = n_buffer-1; j <= nyy; j++) {
				for (k = n_buffer-1; k <= nzz; k++) {  

					/* flux parameter */
					rho = U1[icube][i][j][k][0];
					U = U1[icube][i][j][k][1]/rho;
					V = U1[icube][i][j][k][2]/rho;
					W = U1[icube][i][j][k][3]/rho;     
					VV = U*U+V*V+W*W;
					P = (U1[icube][i][j][k][4]-0.5*rho*VV)*(K-1);
                    T = P/rho/R;
					
					U3_sum[icube][i][j][k] = U3_sum[icube][i][j][k]+T;
                    U4_sum[icube][i][j][k] = U4_sum[icube][i][j][k]+T*T;

				}
			}
		}
	}

#pragma omp parallel for private(i,j,k)	
	for (icube = 1; icube < ncube; icube++) { 
		for (i = n_buffer-1; i <= nxx; i++) {
			for (j = n_buffer-1; j <= nyy; j++) {
				for (k = n_buffer-1; k <= nzz; k++) {  

					T_mean[icube][i][j][k] = U3_sum[icube][i][j][k]/itemp;
                    
					T_int[icube][i][j][k] = sqrt(fabs(U4_sum[icube][i][j][k]/itemp - T_mean[icube][i][j][k]*T_mean[icube][i][j][k]));

                    // H = (U4_sum[icube][i][j][k]/itemp - T_mean[icube][i][j][k]*T_mean[icube][i][j][k]);
                    
                    // if(T_mean[icube][i][j][k] < 298.0591)
                    // printf("%f\t%f\t%f\n",T_mean[icube][i][j][k],U3_sum[icube][i][j][k],itemp*1.0);
                    
				}
			}
		}
	}



	
// =========================== Output averaged result =========================== //

	if ( step%dp_step == 0 ) {

		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];
		MPI_Offset disp;

		int nx_out = NcubeX+2;
		int ny_out = NcubeY+2;
		int nz_out = NcubeZ+2;

		int ncube_out = MPI_Ncube;

		int q_gdisp[np];

		istart = 0;
		q_gdisp[np-1] = 0; 
	
		for (i = 0; i < np-1; i++) { 

			for (icube = 1; icube <= MPI_Ncube; icube++) {

				if (rank_map[0][icube] == i) {

					icount = icube;
					break;
                
				}

			}

			q_gdisp[i] = (MPI_Offset)icount*(MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*num_variable_output;

		}


		char data[100];
		
		float (*Solution) = new float[(ncube-1)*nx_out*ny_out*nz_out*num_variable_output];

		sprintf(data,"qBCM_function_""%0.5d"".q",step);

		MPI_File fh1_function;

		MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1_function ) ; 



		if (myid == 0) {

			FILE *fptr_solution_function;

			fptr_solution_function = fopen(data,"wb");

			fwrite(&ncube_out, sizeof(int), 1,fptr_solution_function);

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				fwrite(&nx_out, sizeof(int), 1,fptr_solution_function);
				fwrite(&ny_out, sizeof(int), 1,fptr_solution_function);
				fwrite(&nz_out, sizeof(int), 1,fptr_solution_function);
				fwrite(&num_variable_output, sizeof(int), 1,fptr_solution_function);

			}

			
			fclose(fptr_solution_function);

		}


		
		icount = -1;
		for (icube = 1; icube < ncube; icube++) {  

            for (k = n_buffer-1; k <= nz+1; k++) {
                for (j = n_buffer-1; j <= ny+1; j++) { 
                    for (i = n_buffer-1; i <= nx+1; i++) { 

						icount = icount + 1;
						Solution[icount] = U_mean[icube][i][j][k];

					}
				}
			}


            for (k = n_buffer-1; k <= nz+1; k++) {
                for (j = n_buffer-1; j <= ny+1; j++) { 
                    for (i = n_buffer-1; i <= nx+1; i++) { 

						icount = icount + 1;
						Solution[icount] = T_mean[icube][i][j][k];

					}
				}
			}
                
            for (k = n_buffer-1; k <= nz+1; k++) {
                for (j = n_buffer-1; j <= ny+1; j++) { 
                    for (i = n_buffer-1; i <= nx+1; i++) { 

						icount = icount + 1;
						Solution[icount] = U_int[icube][i][j][k];

					}
				}
			}


            for (k = n_buffer-1; k <= nz+1; k++) {
                for (j = n_buffer-1; j <= ny+1; j++) { 
                    for (i = n_buffer-1; i <= nx+1; i++) { 

						icount = icount + 1;
						Solution[icount] = T_int[icube][i][j][k];

					}
				}
			}

		}
		
		disp =  ((MPI_Offset)MPI_Ncube*4)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+q_gdisp[myid]*(MPI_Offset)sizeof(float);


		MPI_File_write_at_all(fh1_function, disp, Solution, (ncube-1)*nx_out*ny_out*nz_out*num_variable_output, MPI_FLOAT, MPI_STATUS_IGNORE);

		MPI_File_close( &fh1_function );


		delete [] Solution;

	}  // ---- if ( step%dp_step == 0 ) ---- //


// =========================== Output averaged result =========================== //
	
	delete [] U_mean;
	delete [] T_mean;
	delete [] U_int;
	delete [] T_int;
	

}
