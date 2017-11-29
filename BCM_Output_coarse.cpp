



#include <mpi.h>
#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"
#include "Pre_selection.h"

extern int Ncube;   
void BCM_Output_coarse
(
// ============================================================================ //
int myid,
int ncube,

int step,

int switch_output,

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
    
    int CR_switch = 2;

	
	int nx_out = NcubeX/CR_switch+2;
	int ny_out = NcubeY/CR_switch+2;
	int nz_out = NcubeZ/CR_switch+2;



	int ncube_out = MPI_Ncube;

	
	int int_size = sizeof(int);
	
	int float_size = sizeof(float);

	int double_dize = sizeof(double);

	MPI_Offset  x_gcount[np], x_gdisp[np], q_gcount[np], q_gdisp[np];


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

        // printf("x_gdisp == %d\t%d\t%d\t\%d\n",myid,i,icube,x_gdisp[i]);
        
	}
    
	
	
	idest = 0;
	
	char data[100];

	if ( switch_output == 1 ) {

		sprintf(data,"BCMdata_coarse""%0.5d"".x",1);
		MPI_File fh0;


		MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh0 ) ; 

		if (myid == idest) {

			FILE *fptr_xyz;

			fptr_xyz = fopen(data,"wb");

			fwrite(&ncube_out, sizeof(int), 1,fptr_xyz);

			for (icube = 0; icube < MPI_Ncube; icube++)  {

				fwrite(&nx_out, sizeof(int), 1,fptr_xyz);
				fwrite(&ny_out, sizeof(int), 1,fptr_xyz);
				fwrite(&nz_out, sizeof(int), 1,fptr_xyz);

			}

			fclose(fptr_xyz);

		}

		float (*Grid) = new float[(ncube-1)*nx_out*ny_out*nz_out*3];

		icount = -1;
        
		for (icube = 1; icube < ncube; icube++) {  

            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 

                        icount = icount+1;

                        Grid[icount] = 0.5*( Xcnt[icube][2*i]+Xcnt[icube][2*i+1] );

                    }
                }
            }


            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 

                        icount = icount+1;

                        Grid[icount] = 0.5*( Ycnt[icube][2*j]+Ycnt[icube][2*j+1] );

                    }
                }
            }


            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 

                        icount = icount+1;

                        Grid[icount] = 0.5*( Zcnt[icube][2*k]+Zcnt[icube][2*k+1] );

                    }
                }
            }
            
		}


		disp = ((MPI_Offset)MPI_Ncube*3)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+x_gdisp[myid]*(MPI_Offset)sizeof(float);

		MPI_File_write_at_all(fh0, disp, Grid, (ncube-1)*nx_out*ny_out*nz_out*3, MPI_FLOAT, MPI_STATUS_IGNORE);

		switch_output = 0;

		delete []Grid;

		MPI_File_close( &fh0 );

	}





	float (*Solution) = new float[(ncube-1)*nx_out*ny_out*nz_out*5+4*(ncube-1)];

	sprintf(data,"qBCMdata_coarse""%0.5d"".q",step);
	MPI_File fh1;

	MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1 ) ; 

	if (myid == 0) {

		FILE *fptr_solution;

		fptr_solution = fopen(data,"wb");
		
		fwrite(&ncube_out, sizeof(int), 1,fptr_solution);

		for (icube = 0; icube < MPI_Ncube; icube++)  {

			fwrite(&nx_out, sizeof(int), 1,fptr_solution);
			fwrite(&ny_out, sizeof(int), 1,fptr_solution);
			fwrite(&nz_out, sizeof(int), 1,fptr_solution);
			
		}

		fclose(fptr_solution);

	}



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
    
        for (k = 0; k <= nz_out-1; k++) {
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                icount = icount + 1;
                
                Solution[icount] = 0.125*\
                    (U1[icube][i][j][k][0]+U1[icube][i+1][j][k][0]+\
                     U1[icube][i][j+1][k][0]+U1[icube][i][j][k+1][0]+\
                     U1[icube][i+1][j+1][k][0]+U1[icube][i+1][j][k+1][0]+\
                     U1[icube][i][j+1][k+1][0]+U1[icube][i+1][j+1][k+1][0]);

                }
            }
        }

        for (k = 0; k <= nz_out-1; k++) {
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                    icount = icount + 1;
                    
                Solution[icount] = 0.125*\
                    (U1[icube][i][j][k][1]+U1[icube][i+1][j][k][1]+\
                     U1[icube][i][j+1][k][1]+U1[icube][i][j][k+1][1]+\
                     U1[icube][i+1][j+1][k][1]+U1[icube][i+1][j][k+1][1]+\
                     U1[icube][i][j+1][k+1][1]+U1[icube][i+1][j+1][k+1][1]);

                }
            }
        }

        for (k = 0; k <= nz_out-1; k++) {
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                    icount = icount + 1;
                    
                Solution[icount] = 0.125*\
                    (U1[icube][i][j][k][2]+U1[icube][i+1][j][k][2]+\
                     U1[icube][i][j+1][k][2]+U1[icube][i][j][k+1][2]+\
                     U1[icube][i+1][j+1][k][2]+U1[icube][i+1][j][k+1][2]+\
                     U1[icube][i][j+1][k+1][2]+U1[icube][i+1][j+1][k+1][2]);


                }
            }
        }

        for (k = 0; k <= nz_out-1; k++) {
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                    icount = icount + 1;
                    
                Solution[icount] = 0.125*\
                    (U1[icube][i][j][k][3]+U1[icube][i+1][j][k][3]+\
                     U1[icube][i][j+1][k][3]+U1[icube][i][j][k+1][3]+\
                     U1[icube][i+1][j+1][k][3]+U1[icube][i+1][j][k+1][3]+\
                     U1[icube][i][j+1][k+1][3]+U1[icube][i+1][j+1][k+1][3]);


                }
            }
        }

        for (k = 0; k <= nz_out-1; k++) {
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                    icount = icount + 1;
                    
                Solution[icount] = 0.125*\
                    (U1[icube][i][j][k][4]+U1[icube][i+1][j][k][4]+\
                     U1[icube][i][j+1][k][4]+U1[icube][i][j][k+1][4]+\
                     U1[icube][i+1][j+1][k][4]+U1[icube][i+1][j][k+1][4]+\
                     U1[icube][i][j+1][k+1][4]+U1[icube][i+1][j+1][k+1][4]);


                }
            }
        } 
        
        
	
	}
	

	disp =  ((MPI_Offset)MPI_Ncube*3)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+q_gdisp[myid]*(MPI_Offset)sizeof(float);

	MPI_File_write_at_all(fh1, disp, Solution, (ncube-1)*nx_out*ny_out*nz_out*5+(ncube-1)*4, MPI_FLOAT, MPI_STATUS_IGNORE);

	MPI_File_close( &fh1 );


// ----------------------------------------------------------------- //
// ------------------------- Previous_step ------------------------- //



	

	delete []Solution;


}