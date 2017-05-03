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
extern int NXbc_l;
extern int NXbc_u;
extern int NYbc_l;
extern int NYbc_u;

void BCM_Abs_XY_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nXbc_l,
int nXbc_u,
int nYbc_l,
int nYbc_u,

double gXmax,
double gXmin,
double gYmax,
double gYmin,

double gdXmax,
double gdYmax,

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],
int (*Ybc_l) = new int[NYbc_l+1],
int (*Ybc_u) = new int[NYbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
)

{

#include "BCM.h"
#include "prm.h"
#include "Pre_selection.h"

    double rho,U,V,W,VV,P,C,T,h,H,E;


    double sr1,sr2,sr3,sr4,sr5;
    double sl1,sl2,sl3,sl4,sl5;
	
	double beta, C_plan, EE, E0;
	
	double xV_in_1, xV_out_1, xSigma_in, xSigma_out;
	double yV_in_1, yV_out_1, ySigma_in, ySigma_out;
	
	double xV_in_0, xV_out_0, xSigma_in_0, xSigma_out_0;
	double yV_in_0, yV_out_0, ySigma_in_0, ySigma_out_0;
    
    double n_abs_xm, n_abs_xp, n_abs_ym, n_abs_yp;
    double lenght_absXm, lenght_absYm;
    double lenght_absXp, lenght_absYp;
    double xmin, xmax, ymin, ymax;
    double abs1, abs2 , abs3, abs4;
    double XX, YY, SML;
    
    
    n_abs_xm = 0.8;
    n_abs_xp = 0.8;
    
    n_abs_ym = 0.8;
    n_abs_yp = 0.8;
    
    
    
    xV_in_1 = 0.0;
    xV_out_1 = 0.0;
    xSigma_in = 0.0;
    xSigma_out = 0.0;
    
    yV_in_1 = 0.0;
    yV_out_1 = 0.0;
    ySigma_in = 0.0;
    ySigma_out = 0.0;
    
    sr1 = sr2 = sr3 = sr4 = sr5 = 0.0;
    sl1 = sl2 = sl3 = sl4 = sl5 = 0.0;
     
	xV_in_0 = 1.15;
	xSigma_in_0 = 0.05;
	xV_out_0 = 1.25;
	xSigma_out_0 = 1.25;
	
	yV_in_0 = 1.25;
	ySigma_in_0 = 1.25;
	yV_out_0 = 1.25;
	ySigma_out_0 = 1.25;

    E0 = P0/(K-1.0)+0.5*rho0*U0*U0;
    

    lenght_absXm = gdXmax * (n_abs_xm*NcubeX-1.0);
    lenght_absYm = gdYmax * (n_abs_ym*NcubeY-1.0);
    
    lenght_absXp = gdXmax * (n_abs_xm*NcubeX-1.0);
    lenght_absYp = gdYmax * (n_abs_ym*NcubeY-1.0);
    
    SML = 1.0e-8;
    
    abs1 = gXmin + lenght_absXm + 0.5*gdXmax + SML;
    abs2 = gXmax - lenght_absXp - 0.5*gdXmax - SML;
    abs3 = gYmin + lenght_absYm + 0.5*gdYmax + SML;
    abs4 = gYmax - lenght_absYp - 0.5*gdYmax - SML;
    
	#pragma omp parallel for private(\
    dx,dy,dz,xmin,ymin,\
	i,j,k,\
	rho,U,V,W,VV,E,P,C,\
	beta,C_plan,\
	XX,YY,\
	xV_in_1,xSigma_in,\
	yV_in_1,ySigma_in,\
	xV_out_1,xSigma_out,\
	yV_out_1,ySigma_out,\
	sl1,sl2,sl3,sl4,sl5,\
	sr1,sr2,sr3,sr4,sr5\
	)

    for (icube = 1; icube < ncube; icube++) {  

		dx = dy = dz = cube_size[icube]/NcubeX;
        
        xmin = Xcube[icube]+0.5*dx;
        ymin = Ycube[icube]+0.5*dy;

		for (i = n_buffer; i <= nx; i++) {
			for (j = n_buffer; j <= ny; j++) {
				for (k = n_buffer; k <= nz; k++) {

                    rho = U1_[icube][i][j][k][0];
                    U = U1_[icube][i][j][k][1]/rho;
                    V = U1_[icube][i][j][k][2]/rho;
                    W = U1_[icube][i][j][k][3]/rho;
                    VV = U*U+V*V+W*W;
                    E = U1_[icube][i][j][k][4];
                    P = (E-0.5*rho*VV)*(K-1);
                    C = K*P/rho;

                    /* preconditioning */
                    beta = max(VV/C,e);

                    C_plan = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);
                    
                    XX = xmin + (i - n_buffer)*dx;
                    YY = ymin + (j - n_buffer)*dy;
                    
                    
                    
                    if( XX < abs1 ) {
                    
                        xV_in_1 = ((abs1-XX)/lenght_absXm)*((abs1-XX)/lenght_absXm)*((abs1-XX)/lenght_absXm);
                        xSigma_in = xV_in_1*xSigma_in_0/Char_D*20.0;
                        
                        xV_in_1 = xV_in_1*xV_in_0*C_plan;
                        

						//if( k == 10 && j == 10) {

							//printf("%d\t%f\t%f\t%f\n",i,(abs1-XX)/lenght_absXm, (abs1-XX)/dx,lenght_absXm/gdXmax);

						//}


						} 
					else if ( XX > abs2 ) {

                        xV_out_1 = ((XX-abs2)/lenght_absXp)*((XX-abs2)/lenght_absXp)*((XX-abs2)/lenght_absXp);
                        xSigma_out = xV_out_1*xSigma_out_0/Char_D*20.0;
                        
                        xV_out_1 = xV_out_1*xV_out_0*C_plan;


						//if( k == 10 && j == 10) {

							//printf("%d\t%f\t%f\t%f\n",i,(XX-abs2)/lenght_absXp, (XX-abs2)/dx,lenght_absXp/gdXmax);

						//}

                        
						} 
					else {
                    
                        xV_in_1 = 0.0;
                        xSigma_in = 0.0;
                        
                        xV_out_1 = 0.0;
                        xSigma_out = 0.0;
                        
                    }

					
                    C_plan = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);
                    
                    if( YY < abs3 ) {
                        
                        yV_in_1 = ((abs3-YY)/lenght_absYm)*((abs3-YY)/lenght_absYm)*((abs3-YY)/lenght_absYm);
                        ySigma_in = yV_in_1*ySigma_in_0/Char_D*20.0;
                        
                        yV_in_1 = yV_in_1*yV_in_0*C_plan;
                        
					    
						} 
					else if( YY > abs4  ) {
                    
                        yV_out_1 = ((YY-abs4)/lenght_absYp)*((YY-abs4)/lenght_absYp)*((YY-abs4)/lenght_absYp);
                        ySigma_out = yV_out_1*ySigma_out_0/Char_D*20.0;
                        
                        yV_out_1 = yV_out_1*yV_out_0*C_plan;
                        
                    
                    } else {
                        
                        yV_in_1 = 0.0;
                        ySigma_in = 0.0;
                    
                        yV_out_1 = 0.0;
                        ySigma_out = 0.0;
                        
                    }
                    
                    
                    sl1 = xV_in_1*( rho  - U1_[icube][i-1][j][k][0])/dx + xSigma_in*(rho   - rho0); 
                    sl2 = xV_in_1*( rho*U- U1_[icube][i-1][j][k][1])/dx + xSigma_in*(rho*U - rho0*U0);
                    sl3 = xV_in_1*( rho*V- U1_[icube][i-1][j][k][2])/dx + xSigma_in*(rho*V - rho0*V0);
                    sl4 = xV_in_1*( rho*W- U1_[icube][i-1][j][k][3])/dx + xSigma_in*(rho*W - rho0*W0);
                    sl5 = xV_in_1*( E    - U1_[icube][i-1][j][k][4])/dx + xSigma_in*(E     - E0);
                    
                    
                    sl1 = sl1 + yV_in_1*( rho  - U1_[icube][i][j+1][k][0])/dy + ySigma_in*(rho   - rho0); 
                    sl2 = sl2 + yV_in_1*( rho*U- U1_[icube][i][j+1][k][1])/dy + ySigma_in*(rho*U - rho0*U0);
                    sl3 = sl3 + yV_in_1*( rho*V- U1_[icube][i][j+1][k][2])/dy + ySigma_in*(rho*V - rho0*V0);
                    sl4 = sl4 + yV_in_1*( rho*W- U1_[icube][i][j+1][k][3])/dy + ySigma_in*(rho*W - rho0*W0);
                    sl5 = sl5 + yV_in_1*( E    - U1_[icube][i][j+1][k][4])/dy + ySigma_in*(E     - E0);
                    
                    
                    
                    
                    sr1 = xV_out_1*( rho  - U1_[icube][i-1][j][k][0])/dx + xSigma_out*(rho   - rho0); 
                    sr2 = xV_out_1*( rho*U- U1_[icube][i-1][j][k][1])/dx + xSigma_out*(rho*U - rho0*U0);
                    sr3 = xV_out_1*( rho*V- U1_[icube][i-1][j][k][2])/dx + xSigma_out*(rho*V - rho0*V0);
                    sr4 = xV_out_1*( rho*W- U1_[icube][i-1][j][k][3])/dx + xSigma_out*(rho*W - rho0*W0);
                    sr5 = xV_out_1*( E    - U1_[icube][i-1][j][k][4])/dx + xSigma_out*(E     - E0);
                    
                    
                    sr1 = sr1 + yV_out_1*( rho  - U1_[icube][i][j-1][k][0])/dy + ySigma_out*(rho   - rho0); 
                    sr2 = sr2 + yV_out_1*( rho*U- U1_[icube][i][j-1][k][1])/dy + ySigma_out*(rho*U - rho0*U0);
                    sr3 = sr3 + yV_out_1*( rho*V- U1_[icube][i][j-1][k][2])/dy + ySigma_out*(rho*V - rho0*V0);
                    sr4 = sr4 + yV_out_1*( rho*W- U1_[icube][i][j-1][k][3])/dy + ySigma_out*(rho*W - rho0*W0);
                    sr5 = sr5 + yV_out_1*( E    - U1_[icube][i][j-1][k][4])/dy + ySigma_out*(E     - E0);
                    
                                               
                    Fabs[icube][i][j][k][0] = -sr1-sl1;
                    Fabs[icube][i][j][k][1] = -sr2-sl2;
                    Fabs[icube][i][j][k][2] = -sr3-sl3;
                    Fabs[icube][i][j][k][3] = -sr4-sl4;
                    Fabs[icube][i][j][k][4] = -sr5-sl5;
                    

                }
            }
         }

      }
      
      
      
      #pragma omp parallel for private(iicube,j,k,rho,U,V,W,VV,P)
      
      for (icube = 1; icube <= nXbc_l; icube++) {  
		
            iicube = Xbc_l[icube];

            for (j = 2; j <= ny; j++) {
                for (k = 2; k <= nz; k++) {  

                    rho = U1_[iicube][2][j][k][0];
                    U = U1_[iicube][2][j][k][1]/rho;
                    V = U1_[iicube][2][j][k][2]/rho;
                    W = U1_[iicube][2][j][k][3]/rho;
                    VV = U*U+V*V+W*W;
                    //P = (U1_[iicube][2][j][k][4]-0.5*rho*VV)*(K-1);

                    rho = rho0;

                    U = Uin;
                    V = V0;
                    W = W0;
                    VV = U*U+V*V+W*W;
                    P = P0;


                    U1_[iicube][1][j][k][0] = rho;
                    U1_[iicube][1][j][k][1] = rho*U;
                    U1_[iicube][1][j][k][2] = rho*V;
                    U1_[iicube][1][j][k][3] = rho*W;
                    U1_[iicube][1][j][k][4] = P/(K-1)+0.5*rho*VV;

                    U1_[iicube][0][j][k][0] = rho;
                    U1_[iicube][0][j][k][1] = rho*U;
                    U1_[iicube][0][j][k][2] = rho*V;
                    U1_[iicube][0][j][k][3] = rho*W;
                    U1_[iicube][0][j][k][4] = P/(K-1)+0.5*rho*VV;

                }
            }
        }	
        #pragma omp barrier
        
        
        #pragma omp parallel for private(iicube,j,k)
        for (icube = 1; icube <= nXbc_u; icube++) {  

            iicube = Xbc_u[icube];
            
            for (j = 2; j <= ny; j++) {
                for (k = 2; k <= nz; k++) {  

                    U1_[iicube][nxx][j][k][0] = U1_[iicube][nx][j][k][0];
                    U1_[iicube][nxx][j][k][1] = U1_[iicube][nx][j][k][1];
                    U1_[iicube][nxx][j][k][2] = U1_[iicube][nx][j][k][2];
                    U1_[iicube][nxx][j][k][3] = U1_[iicube][nx][j][k][3];
                    U1_[iicube][nxx][j][k][4] = U1_[iicube][nx][j][k][4];

                    U1_[iicube][nxxx][j][k][0] = U1_[iicube][nx][j][k][0];
                    U1_[iicube][nxxx][j][k][1] = U1_[iicube][nx][j][k][1];
                    U1_[iicube][nxxx][j][k][2] = U1_[iicube][nx][j][k][2];
                    U1_[iicube][nxxx][j][k][3] = U1_[iicube][nx][j][k][3];
                    U1_[iicube][nxxx][j][k][4] = U1_[iicube][nx][j][k][4];

                }
            }

        }

        #pragma omp barrier        
          
        
        #pragma omp parallel for private(iicube,i,k)
        
        for (icube = 1; icube <= nYbc_l; icube++) {  

            iicube = Ybc_l[icube];

                for (i = 0; i <= nxxx; i++) {
                    for (k = 2; k <= nz; k++) {  

                        U1_[iicube][i][1][k][0] = U1_[iicube][i][2][k][0];
                        U1_[iicube][i][1][k][1] = U1_[iicube][i][2][k][1];
                        U1_[iicube][i][1][k][2] = U1_[iicube][i][2][k][2];
                        U1_[iicube][i][1][k][3] = U1_[iicube][i][2][k][3];
                        U1_[iicube][i][1][k][4] = U1_[iicube][i][2][k][4];

                        U1_[iicube][i][0][k][0] = U1_[iicube][i][2][k][0];
                        U1_[iicube][i][0][k][1] = U1_[iicube][i][2][k][1];
                        U1_[iicube][i][0][k][2] = U1_[iicube][i][2][k][2];
                        U1_[iicube][i][0][k][3] = U1_[iicube][i][2][k][3];
                        U1_[iicube][i][0][k][4] = U1_[iicube][i][2][k][4];

                    }
                }

        }		

        #pragma omp barrier        
        
    
    #pragma omp parallel for private(iicube,i,k)

        
	for (icube = 1; icube <= nYbc_u; icube++) {  

		iicube = Ybc_u[icube];

		for (i = 0; i <= nxxx; i++) {
			for (k = 2; k <= nz; k++) {  

				U1_[iicube][i][nyy][k][0] = U1_[iicube][i][ny][k][0];
				U1_[iicube][i][nyy][k][1] = U1_[iicube][i][ny][k][1];
				U1_[iicube][i][nyy][k][2] = U1_[iicube][i][ny][k][2];
				U1_[iicube][i][nyy][k][3] = U1_[iicube][i][ny][k][3];
				U1_[iicube][i][nyy][k][4] = U1_[iicube][i][ny][k][4];

				U1_[iicube][i][nyyy][k][0] = U1_[iicube][i][ny][k][0];
				U1_[iicube][i][nyyy][k][1] = U1_[iicube][i][ny][k][1];
				U1_[iicube][i][nyyy][k][2] = U1_[iicube][i][ny][k][2];
				U1_[iicube][i][nyyy][k][3] = U1_[iicube][i][ny][k][3];
				U1_[iicube][i][nyyy][k][4] = U1_[iicube][i][ny][k][4];

			}
		}
		
	}

    #pragma omp barrier       


}