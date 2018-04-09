#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;    
extern int N_wallcube;    

void BCM_Slice_output
(
// =================================================== //
int myid,
int ncube,

double Xp,

int step,

int slice_normal,

char slice_name[100],

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*cube_size) = new double[Ncube],

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size]
// ================================================== //
)

{
    #include "BCM.h"
	#include "MPI_prm.h"
    
    
    MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Status istat[8];
	MPI_Offset disp;
    
    int nx_out = NcubeX+1;
	int ny_out = NcubeY+1;
	int nz_out = NcubeZ+1;
    
    int int_size = sizeof(int);
	int float_size = sizeof(float);
	int double_dize = sizeof(double);
    
    int ii,jj,kk;
    
    double XX, YY, ZZ, SML;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double rho,U,V,W,VV,P,C,T,E;
    
    int io,jo,ko;
    double u1,u2,u3,u4,u5;
    double gu1,gu2,gu3,gu4,gu5;

    int i_count = 0;
    int g_count = 0;
    int max_g_count = 0;
    
    int ncube_local, ncube_out;
    
    int x_gcount[np];
    
    MPI_Offset  x_gdisp[np], g_gcount[np], q_gdisp[np];
    
    
    for (icube = 1; icube < ncube; icube++) {  

        dx = cube_size[icube];
    
        xmin = Xcube[icube];
        ymin = Ycube[icube];
        zmin = Zcube[icube];
        
        xmax = xmin+dx;
        ymax = ymin+dx;
        zmax = zmin+dx;
        
        if( (Xp < zmax) && (Xp >= zmin) ) {
            
            i_count = i_count + 1;
            
        }

    }
    
    ncube_local = i_count;
    
    MPI_Allreduce ((void*)&i_count,(void*)&g_count,1,MPI_INT,MPI_SUM,comm );
    
    ncube_out = g_count;
    
    idest = 0;
    
    itag=400;
    
    if(myid > 0) {
        
        MPI_Send(&i_count, 1, MPI_INT, idest, itag, comm);
        
    }
    else {
        
        for (i = 1; i <= np-1; i++) {
            
            MPI_Recv(&x_gcount[i], 1, MPI_INT, i, itag, comm, istat);
            
        }
        
        x_gcount[0] = ncube_local;
        
    }
    
    MPI_Bcast( x_gcount, np, MPI_INT, 0, comm );
    
    
    icount = 0;        
    x_gdisp[0] = 0; 
    q_gdisp[0] = 0;
	
	for (i = 1; i <= np-1; i++) { 
    
        icount = icount+x_gcount[i-1];
    
        x_gdisp[i] = (MPI_Offset)icount*(MPI_Offset)nx_out*(MPI_Offset)ny_out*2;
		q_gdisp[i] = (MPI_Offset)icount*(MPI_Offset)nx_out*(MPI_Offset)ny_out*(MPI_Offset)nz_out*5+(MPI_Offset)icount*4;

	}
    
    
    
    // --------------------------------------------------- //
    // -------------------- Grid file -------------------- //
    
	MPI_File fh0;
    
    MPI_File_open( MPI_COMM_WORLD, slice_name,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh0 ) ; 

    if (myid == idest) {

        FILE *fptr_xyz;

        fptr_xyz = fopen(slice_name,"wb");

        fwrite(&ncube_out, sizeof(int), 1,fptr_xyz);

        for (icube = 0; icube < ncube_out; icube++)  {

            fwrite(&nx_out, sizeof(int), 1,fptr_xyz);
            fwrite(&ny_out, sizeof(int), 1,fptr_xyz);

        }

        fclose(fptr_xyz);

	}
    
    
    float (*Grid) = new float[ncube_local*nx_out*ny_out*2];
    
    icount = -1;
    

    for (icube = 1; icube < ncube; icube++) {  
        
        dx = cube_size[icube];
    
        xmin = Xcube[icube];
        ymin = Ycube[icube];
        zmin = Zcube[icube];
        
        xmax = xmin+dx;
        ymax = ymin+dx;
        zmax = zmin+dx;
        
        if( (Xp < zmax) && (Xp >= zmin) ) {
        
        
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                    icount = icount+1;

                    Grid[icount] = 0.5*( Xcnt[icube][i+1]+Xcnt[icube][i+2] );

                }
            }
       
            for (j = 0; j <= ny_out-1; j++) { 
                for (i = 0; i <= nx_out-1; i++) { 

                    icount = icount+1;

                    Grid[icount] = 0.5*( Ycnt[icube][j+1]+Ycnt[icube][j+2] );

                }
            }

        }    // ==== if( (Xp <= zmax) && (Xp >= zmin) ) ==== //
        
    }
    
    
    // printf("%d\t%d\t%d\t%d\t%d\n",myid,ncube_local,ncube_out,icount,(icount+1)/17/17/17/3);
    
    
    disp = ((MPI_Offset)ncube_out*2)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+x_gdisp[myid]*(MPI_Offset)sizeof(float);

    MPI_File_write_at_all(fh0, disp, Grid, ncube_local*nx_out*ny_out*2, MPI_FLOAT, MPI_STATUS_IGNORE);

    delete []Grid;

    MPI_File_close( &fh0 );
    
    
    // -------------------- Grid file -------------------- //
    // --------------------------------------------------- //
    
    
    // ------------------------------------------------ //
    // -------------------- Q file -------------------- //
    
    char data[100];
    
    sprintf(data,"qslice""%0.5d"".q",step);
	MPI_File fh1;

	MPI_File_open( MPI_COMM_WORLD, data,  MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh1 ) ; 
    
    if (myid == idest) {

        FILE *fptr_xyz;

        fptr_xyz = fopen(data,"wb");

        fwrite(&ncube_out, sizeof(int), 1,fptr_xyz);

        for (icube = 0; icube < ncube_out; icube++)  {

            fwrite(&nx_out, sizeof(int), 1,fptr_xyz);
            fwrite(&ny_out, sizeof(int), 1,fptr_xyz);
            fwrite(&nz_out, sizeof(int), 1,fptr_xyz);

        }

        fclose(fptr_xyz);

	}
    
    
    float (*Solution) = new float[ncube_local*nx_out*ny_out*nz_out*5+4*(ncube-1)];
    
    
    icount = -1;
    

    for (icube = 1; icube < ncube; icube++) {  
        
        dx = cube_size[icube];
    
        xmin = Xcube[icube];
        ymin = Ycube[icube];
        zmin = Zcube[icube];
        
        xmax = xmin+dx;
        ymax = ymin+dx;
        zmax = zmin+dx;
        
        if( (Xp < zmax) && (Xp >= zmin) ) {
            
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
                    
                    ii = i+1;
                    jj = j+1;
                    kk = k+1;

                    icount = icount + 1;
                    
                    Solution[icount] = 0.125*\
                        (U1[icube][ii  ][jj  ][kk  ][0]+U1[icube][ii+1][jj  ][kk  ][0]+\
                         U1[icube][ii  ][jj+1][kk  ][0]+U1[icube][ii  ][jj  ][kk+1][0]+\
                         U1[icube][ii+1][jj+1][kk  ][0]+U1[icube][ii+1][jj  ][kk+1][0]+\
                         U1[icube][ii  ][jj+1][kk+1][0]+U1[icube][ii+1][jj+1][kk+1][0]);

                    }
                }
            }

            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 
                    
                    ii = i+1;
                    jj = j+1;
                    kk = k+1;

                    icount = icount + 1;
                    
                    Solution[icount] = 0.125*\
                        (U1[icube][ii  ][jj  ][kk  ][1]+U1[icube][ii+1][jj  ][kk  ][1]+\
                         U1[icube][ii  ][jj+1][kk  ][1]+U1[icube][ii  ][jj  ][kk+1][1]+\
                         U1[icube][ii+1][jj+1][kk  ][1]+U1[icube][ii+1][jj  ][kk+1][1]+\
                         U1[icube][ii  ][jj+1][kk+1][1]+U1[icube][ii+1][jj+1][kk+1][1]);

                    }
                }
            }

            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 
                    
                    ii = i+1;
                    jj = j+1;
                    kk = k+1;

                    icount = icount + 1;
                    
                    Solution[icount] = 0.125*\
                        (U1[icube][ii  ][jj  ][kk  ][2]+U1[icube][ii+1][jj  ][kk  ][2]+\
                         U1[icube][ii  ][jj+1][kk  ][2]+U1[icube][ii  ][jj  ][kk+1][2]+\
                         U1[icube][ii+1][jj+1][kk  ][2]+U1[icube][ii+1][jj  ][kk+1][2]+\
                         U1[icube][ii  ][jj+1][kk+1][2]+U1[icube][ii+1][jj+1][kk+1][2]);


                    }
                }
            }

            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 
                    
                    ii = i+1;
                    jj = j+1;
                    kk = k+1;

                    icount = icount + 1;
                    
                    Solution[icount] = 0.125*\
                        (U1[icube][ii  ][jj  ][kk  ][3]+U1[icube][ii+1][jj  ][kk  ][3]+\
                         U1[icube][ii  ][jj+1][kk  ][3]+U1[icube][ii  ][jj  ][kk+1][3]+\
                         U1[icube][ii+1][jj+1][kk  ][3]+U1[icube][ii+1][jj  ][kk+1][3]+\
                         U1[icube][ii  ][jj+1][kk+1][3]+U1[icube][ii+1][jj+1][kk+1][3]);


                    }
                }
            }

            for (k = 0; k <= nz_out-1; k++) {
                for (j = 0; j <= ny_out-1; j++) { 
                    for (i = 0; i <= nx_out-1; i++) { 
                    
                    ii = i+1;
                    jj = j+1;
                    kk = k+1;

                    icount = icount + 1;
                    
                    Solution[icount] = 0.125*\
                        (U1[icube][ii  ][jj  ][kk  ][4]+U1[icube][ii+1][jj  ][kk  ][4]+\
                         U1[icube][ii  ][jj+1][kk  ][4]+U1[icube][ii  ][jj  ][kk+1][4]+\
                         U1[icube][ii+1][jj+1][kk  ][4]+U1[icube][ii+1][jj  ][kk+1][4]+\
                         U1[icube][ii  ][jj+1][kk+1][4]+U1[icube][ii+1][jj+1][kk+1][4]);


                    }
                }
            } 
        
        }    // ==== if( (Xp <= zmax) && (Xp >= zmin) ) ==== //
        
    }
    
    
    disp = ((MPI_Offset)ncube_out*3)*(MPI_Offset)sizeof(int)+1*(MPI_Offset)sizeof(int)+q_gdisp[myid]*(MPI_Offset)sizeof(float);

    MPI_File_write_at_all(fh1, disp, Solution, ncube_local*nx_out*ny_out*nz_out*5+ncube_local*4, MPI_FLOAT, MPI_STATUS_IGNORE);

	MPI_File_close( &fh1 );
    
    delete []Solution;
    
    
    
    

    // -------------------- Q file -------------------- //
    // ------------------------------------------------ //
    
    
    
    
    
}