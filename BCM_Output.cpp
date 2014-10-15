



#include <mpi.h>
#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;   
void BCM_Output
(
// ============================================================================ //
int myid,
int ncube,

int step,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*cube_size) = new double[Ncube],


double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size]

// ============================================================================ //
)

{
	#include "BCM.h"
	#include "MPI_prm.h"

	
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
		q_gcount[i] = icount*X_size*Y_size*Z_size*Ndim;

		if (i < (np-1)) {
			
			x_gdisp[i+1] = x_gdisp[i]+istart*nx_out*ny_out*nz_out*3;
			q_gdisp[i+1] = q_gdisp[i]+istart*nx_out*ny_out*nz_out*5+istart*4;

		}

		//printf("q_gdisp == %d\t%d\t%d\t%d\t\%d\n",myid,i,istart,istart,x_gdisp[i]);

	}
	
	
	idest = 0;
	
	icount = x_gcount[myid];
	
	char data[100];
	FILE *fptr;
	sprintf(data,"BCMdata""%0.5d"".x",1);
	MPI_File fh0;


	MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh0 ) ; 

	//if (myid == idest) {

		FILE *fptr_xyz;

		fptr_xyz = fopen(data,"wb");
		
		fwrite(&ncube_out, sizeof(int), 1,fptr_xyz);

		for (icube = 0; icube < MPI_Ncube; icube++)  {

			fwrite(&nx_out, sizeof(int), 1,fptr_xyz);
			fwrite(&ny_out, sizeof(int), 1,fptr_xyz);
			fwrite(&nz_out, sizeof(int), 1,fptr_xyz);
			
		}

		fclose(fptr_xyz);

	//}

	double (*Grid) = new double[(ncube-1)*nx_out*ny_out*nz_out*3];

	icount = -1;
	for (icube = 1; icube < ncube; icube++) {  

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount+1;

					Grid[icount] = Xcnt[icube][k];
					
				}
			}
		}


		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount+1;

					Grid[icount] = Ycnt[icube][j];

				}
			}
		}


		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount+1;

					Grid[icount] = Zcnt[icube][i];
					
				}
			}
		}

	}

	
	disp = (MPI_Ncube*3+1)*sizeof(int)+x_gdisp[myid]*sizeof(double);
	
	MPI_File_set_view(fh0, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	MPI_File_write(fh0, Grid, (ncube-1)*nx_out*ny_out*nz_out*3, MPI_DOUBLE, MPI_STATUS_IGNORE);



	delete []Grid;

	MPI_File_close( &fh0 );



	double (*Solution) = new double[(ncube-1)*nx_out*ny_out*nz_out*5+4*(ncube-1)];

	sprintf(data,"qBCMdata""%0.5d"".q",step);
	MPI_File fh1;

	MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1 ) ; 

	//if (myid == idest) {

		FILE *fptr_solution;

		fptr_solution = fopen(data,"wb");
		
		fwrite(&ncube_out, sizeof(int), 1,fptr_solution);

		for (icube = 0; icube < MPI_Ncube; icube++)  {

			fwrite(&nx_out, sizeof(int), 1,fptr_solution);
			fwrite(&ny_out, sizeof(int), 1,fptr_solution);
			fwrite(&nz_out, sizeof(int), 1,fptr_solution);

		}

		fclose(fptr_solution);

	//}



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
					Solution[icount] = U1[icube][k][j][i][0];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					Solution[icount] = U1[icube][k][j][i][1];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					Solution[icount] = U1[icube][k][j][i][2];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					Solution[icount] = U1[icube][k][j][i][3];

				}
			}
		}

		for (i = n_buffer-1; i <= nx+1; i++) {
			for (j = n_buffer-1; j <= ny+1; j++) { 
				for (k = n_buffer-1; k <= nz+1; k++) { 

					icount = icount + 1;
					Solution[icount] = U1[icube][k][j][i][4];

				}
			}
		}
	
	}
	
	disp = (MPI_Ncube*3+1)*sizeof(int)+q_gdisp[myid]*sizeof(double);

	MPI_File_set_view(fh1, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

	MPI_File_write(fh1, Solution, (ncube-1)*nx_out*ny_out*nz_out*5+(ncube-1)*4, MPI_DOUBLE, MPI_STATUS_IGNORE);





	delete []Solution;

	MPI_File_close( &fh1 );


}
