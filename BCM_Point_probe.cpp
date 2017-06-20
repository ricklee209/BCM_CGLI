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

void BCM_Point_probe
(
// =================================================== //
int myid,
int ncube,
int n_wallcube,

double Xp,
double Yp,
double Zp,

char probe_name[100],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*cube_size) = new double[Ncube],
int (*csl) = new int[Ncube],

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

int (*adj_number)[5][7] = new int[Ncube][5][7],

int (*wallcube) = new int[Ncube]
// ================================================== //
)

{

    #include "BCM.h"
    #include "prm.h"
    #include "MPI_prm.h"

    double XX, YY, ZZ, SML;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double rho,U,V,W,VV,P,C,T,E;
    
    int io,jo,ko;
    double u1,u2,u3,u4,u5;
    double gu1,gu2,gu3,gu4,gu5;

    int i_switch = 0;
    int g_switch = 0;
    
    io = -1;
    jo = -1;
    ko = -1;

#pragma omp parallel for private(i,j,k,dx,dy,dz,xmin,ymin,zmin,XX,YY,ZZ,io,jo,ko,\
rho,U,V,W,VV,E,P,T,u1,u2,u3,u4,u5)
    for (icube = 1; icube < ncube; icube++) {  

            dx = dy = dz = cube_size[icube]/NcubeX;
            
            xmin = Xcube[icube];
            ymin = Ycube[icube];
            zmin = Zcube[icube];

            for (i = n_buffer; i <= nx; i++) {
                for (j = n_buffer; j <= ny; j++) {
                    for (k = n_buffer; k <= nz; k++) {

                        XX = xmin + (i - n_buffer)*dx;
                        YY = ymin + (j - n_buffer)*dy;
                        ZZ = zmin + (k - n_buffer)*dz;
                        
                        if( Xp >= XX && Xp <= XX+dx && \
                            Yp >= YY && Yp <= YY+dy && \
                            Zp >= ZZ && Zp <= ZZ+dz ) {
                            
                            io = i;
                            jo = j;
                            ko = k;
                            
                            i_switch = 1;
                            
                            rho = U1[icube][i][j][k][0];
                            U = U1[icube][i][j][k][1]/rho;
                            V = U1[icube][i][j][k][2]/rho;
                            W = U1[icube][i][j][k][3]/rho;
                            VV = U*U+V*V+W*W;
                            E = U1[icube][i][j][k][4];
                            P = (E-0.5*rho*VV)*(K-1);
                            T = P/rho/R;
                            
                            u1 = P;
                            u2 = U;
                            u3 = V;
                            u4 = W;
                            u5 = T;
                            
                          }
                        
                        
                        
                }
             }
         }
         
    }

    
    MPI_Comm comm;
	comm=MPI_COMM_WORLD;

    MPI_Allreduce ((void*)&i_switch, (void*)&g_switch, 1, MPI_INT, MPI_SUM, comm );

    if( g_switch > 0 ) {

        MPI_Allreduce ((void*)&u1,(void*)&gu1,1,MPI_DOUBLE,MPI_SUM,comm );
        MPI_Allreduce ((void*)&u2,(void*)&gu2,1,MPI_DOUBLE,MPI_SUM,comm );
        MPI_Allreduce ((void*)&u3,(void*)&gu3,1,MPI_DOUBLE,MPI_SUM,comm );
        MPI_Allreduce ((void*)&u4,(void*)&gu4,1,MPI_DOUBLE,MPI_SUM,comm );
        MPI_Allreduce ((void*)&u5,(void*)&gu5,1,MPI_DOUBLE,MPI_SUM,comm );
        
        gu1 = gu1/g_switch;
        gu2 = gu2/g_switch;
        gu3 = gu3/g_switch;
        gu4 = gu4/g_switch;
        gu5 = gu5/g_switch;
        
        
        if( myid == 0) {
        
            FILE *fptr;
            fptr = fopen(probe_name,"a"); 
            fprintf(fptr,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",gu1,gu2,gu3,gu4,gu5);
            fclose(fptr);
        
        }
     
    }


}
 