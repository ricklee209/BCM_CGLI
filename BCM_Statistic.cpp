




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

double (*Pall)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*VVall)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]

// =================================================== //
)

{


// ====== Multi-zone structured scalar function file format ====== //

	
#include "BCM.h"
#include "prm.h"
#include "MPI_prm.h"

	double (*Pm)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*VVm)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];

	double (*Pp)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size];


	double rho,U,V,W,VV,P,C,T,H;

	int Average_step = statistic_step;

	int num_variable_output = 2; // output dp and Vmean two values //

	int itemp = 1;

	#pragma omp parallel for private(i,j,k,rho,U,V,W,VV,P)	
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
					
					Pall[icube][i][j][k] = Pall[icube][i][j][k]+P;

					VVall[icube][i][j][k] = VVall[icube][i][j][k]+sqrt(VV);

				}
			}
		}
	}


#pragma omp parallel for private(i,j,k)	
	for (icube = 1; icube < ncube; icube++) { 
		for (i = n_buffer-1; i <= nxx; i++) {
			for (j = n_buffer-1; j <= nyy; j++) {
				for (k = n_buffer-1; k <= nzz; k++) {  

					Pm[icube][i][j][k] = Pall[icube][i][j][k]/(step*1.0-start_step*1.0+1.0);

					VVm[icube][i][j][k] = VVall[icube][i][j][k]/(step*1.0-start_step*1.0+1.0);

				}
			}
		}
	}



	
// =========================== Output averaged result =========================== //

	if ( step%dp_step == 0 ) {


#pragma omp parallel for private(i,j,k,rho,U,V,W,VV,P)	
		for (icube = 1; icube < ncube; icube++) { 
			for (i = n_buffer-1; i <= nxx; i++) {
				for (j = n_buffer-1; j <= nyy; j++) {
					for (k = n_buffer-1; k <= nzz; k++) {  

						rho = U1[icube][i][j][k][0];
						U = U1[icube][i][j][k][1]/rho;
						V = U1[icube][i][j][k][2]/rho;
						W = U1[icube][i][j][k][3]/rho;     
						VV = U*U+V*V+W*W;
						P = (U1[icube][i][j][k][4]-0.5*rho*VV)*(K-1);

						Pp[icube][i][j][k] = P-Pm[icube][i][j][k];

					}
				}
			}
		}

		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];
		MPI_Offset disp;

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

				q_gdisp[i+1] = q_gdisp[i]+(MPI_Offset)icount*(MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*num_variable_output;


			}

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

			for (i = n_buffer-1; i <= nx+1; i++) {
				for (j = n_buffer-1; j <= ny+1; j++) { 
					for (k = n_buffer-1; k <= nz+1; k++) { 

						icount = icount + 1;
						Solution[icount] = Pp[icube][k][j][i];

					}
				}
			}


			for (i = n_buffer-1; i <= nx+1; i++) {
				for (j = n_buffer-1; j <= ny+1; j++) { 
					for (k = n_buffer-1; k <= nz+1; k++) { 

						icount = icount + 1;
						Solution[icount] = VVm[icube][k][j][i];

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
	
	delete [] Pm;
	delete [] VVm;
	delete [] Pp;
	

}
