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
extern int NZbc_l;
extern int NZbc_u;

void BCM_Abs_XYZ_boundary_condition
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
int nZbc_l,
int nZbc_u,

double gXmax,
double gXmin,
double gYmax,
double gYmin,
double gZmax,
double gZmin,


double gdXmax,
double gdYmax,
double gdZmax,

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],
int (*Ybc_l) = new int[NYbc_l+1],
int (*Ybc_u) = new int[NYbc_u+1],
int (*Zbc_l) = new int[NYbc_l+1],
int (*Zbc_u) = new int[NYbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*CFL_tau)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]
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
    double zV_in_1, zV_out_1, zSigma_in, zSigma_out;
	
	double xV_in_0, xV_out_0, xSigma_in_0, xSigma_out_0;
	double yV_in_0, yV_out_0, ySigma_in_0, ySigma_out_0;
    double zV_in_0, zV_out_0, zSigma_in_0, zSigma_out_0;
    
    double n_abs_xm, n_abs_xp, n_abs_ym, n_abs_yp, n_abs_zm, n_abs_zp;
    double lenght_absXm, lenght_absYm, lenght_absZm;
    double lenght_absXp, lenght_absYp, lenght_absZp;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double abs1, abs2 , abs3, abs4, abs5, abs6;
    double XX, YY, ZZ, SML;
    
    double U1,U2,U3,U4,U5;
    
    
    n_abs_xm = 0.4;
    n_abs_xp = 0.4;
    
    n_abs_ym = 0.4;
    n_abs_yp = 0.4;
    
    n_abs_zm = 0.4;
    n_abs_zp = 0.4;
    
    
    xV_in_1 = 0.0;
    xV_out_1 = 0.0;
    xSigma_in = 0.0;
    xSigma_out = 0.0;
    
    yV_in_1 = 0.0;
    yV_out_1 = 0.0;
    ySigma_in = 0.0;
    ySigma_out = 0.0;
    
    zV_in_1 = 0.0;
    zV_out_1 = 0.0;
    zSigma_in = 0.0;
    zSigma_out = 0.0;
    
    sr1 = sr2 = sr3 = sr4 = sr5 = 0.0;
    sl1 = sl2 = sl3 = sl4 = sl5 = 0.0;
     
	xV_in_0 = 1.15;
	xSigma_in_0 = 0.4;
	xV_out_0 = 1.25;
	xSigma_out_0 = 0.4;
	
	yV_in_0 = 0.0;
	ySigma_in_0 = 0.4;
	yV_out_0 = 0.0;
	ySigma_out_0 = 0.4;
    
    zV_in_0 = 0.0;
	zSigma_in_0 = 0.4;
	zV_out_0 = 0.0;
	zSigma_out_0 = 0.4;

    E0 = P0/(K-1.0)+0.5*rho0*U0*U0;
    

    lenght_absXm = gdXmax * (n_abs_xm*NcubeX-1.0);
    lenght_absYm = gdYmax * (n_abs_ym*NcubeY-1.0);
    lenght_absZm = gdZmax * (n_abs_zm*NcubeZ-1.0);
    
    lenght_absXp = gdXmax * (n_abs_xp*NcubeX-1.0);
    lenght_absYp = gdYmax * (n_abs_yp*NcubeY-1.0);
    lenght_absZp = gdZmax * (n_abs_zp*NcubeZ-1.0);
    
    SML = 1.0e-8;
    
    abs1 = gXmin + lenght_absXm + 0.5*gdXmax + SML;
    abs2 = gXmax - lenght_absXp - 0.5*gdXmax - SML;
    abs3 = gYmin + lenght_absYm + 0.5*gdYmax + SML;
    abs4 = gYmax - lenght_absYp - 0.5*gdYmax - SML;
    abs5 = gZmin + lenght_absZm + 0.5*gdZmax + SML;
    abs6 = gZmax - lenght_absZp - 0.5*gdZmax - SML;
    
    
    
    
    
    
    
    
    
      for (icube = 1; icube <= nXbc_l; icube++) {  
		
            iicube = Xbc_l[icube];
            
            dx = dy = dz = cube_size[iicube]/NcubeX;
            
            for (j = 2; j <= ny; j++) {
                for (k = 2; k <= nz; k++) {  
                
                    rho = U1_[iicube][2][j][k][0];
                    U = U1_[iicube][2][j][k][1]/rho;
                    V = U1_[iicube][2][j][k][2]/rho;
                    W = U1_[iicube][2][j][k][3]/rho;
                    VV = U*U+V*V+W*W;
                    P = (U1_[iicube][2][j][k][4]-0.5*rho*VV)*(K-1);

                    U = U0;
                    V = V0;
                    W = W0;
                    
                    rho = rho0;
                    VV = U*U+V*V+W*W;
                    
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

    
    
    
    
    
	#pragma omp parallel for private(\
    dx,dy,dz,xmin,ymin,zmin,\
	i,j,k,\
	rho,U,V,W,VV,E,P,C,\
	beta,C_plan,\
	XX,YY,ZZ,\
	xV_in_1,xSigma_in,\
	yV_in_1,ySigma_in,\
    zV_in_1,zSigma_in,\
	xV_out_1,xSigma_out,\
	yV_out_1,ySigma_out,\
    zV_out_1,zSigma_out,\
	sl1,sl2,sl3,sl4,sl5,\
	sr1,sr2,sr3,sr4,sr5,\
    U1,U2,U3,U4,U5\
	)

    for (icube = 1; icube < ncube; icube++) {  

		dx = dy = dz = cube_size[icube]/NcubeX;
        
        xmin = Xcube[icube]+0.5*dx;
        ymin = Ycube[icube]+0.5*dy;
        zmin = Zcube[icube]+0.5*dz;

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

                    XX = xmin + (i - n_buffer)*dx;
                    YY = ymin + (j - n_buffer)*dy;
                    ZZ = zmin + (k - n_buffer)*dz;
                    
                    C_plan = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);
                    
                    if( XX < abs1 ) {
                    
                        xV_in_1 = ((abs1-XX)/lenght_absXm)*((abs1-XX)/lenght_absXm)*((abs1-XX)/lenght_absXm);
                        xSigma_in = xV_in_1*xSigma_in_0/Char_D;
                        
                        xV_in_1 = xV_in_1*xV_in_0*C_plan;
                        
                        xV_out_1 = 0.0;
                        xSigma_out = 0.0;
                        
                        CFL_tau[icube][i][j][k] = -DTau_CFL;


						} 
					else if ( XX > abs2 ) {

                        xV_out_1 = ((XX-abs2)/lenght_absXp)*((XX-abs2)/lenght_absXp)*((XX-abs2)/lenght_absXp);
                        xSigma_out = xV_out_1*xSigma_out_0/Char_D;
                        
                        xV_out_1 = xV_out_1*xV_out_0*C_plan;

                        xV_in_1 = 0.0;
                        xSigma_in = 0.0;
                        
                        CFL_tau[icube][i][j][k] = -DTau_CFL;

                        
						} 
					else {
                    
                        xV_in_1 = 0.0;
                        xSigma_in = 0.0;
                        
                        xV_out_1 = 0.0;
                        xSigma_out = 0.0;
                        
                    }

					
                    // C_plan = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);
                    
                    if( YY < abs3 ) {
                        
                        yV_in_1 = ((abs3-YY)/lenght_absYm)*((abs3-YY)/lenght_absYm)*((abs3-YY)/lenght_absYm);
                        ySigma_in = yV_in_1*ySigma_in_0/Char_D;
                        
                        yV_in_1 = yV_in_1*yV_in_0*C_plan;
                        
                        yV_out_1 = 0.0;
                        ySigma_out = 0.0;
                        
                        CFL_tau[icube][i][j][k] = -DTau_CFL;
                        
					    
						} 
					else if( YY > abs4  ) {
                    
                        yV_out_1 = ((YY-abs4)/lenght_absYp)*((YY-abs4)/lenght_absYp)*((YY-abs4)/lenght_absYp);
                        ySigma_out = yV_out_1*ySigma_out_0/Char_D;
                        
                        yV_out_1 = yV_out_1*yV_out_0*C_plan;
                        
                        yV_in_1 = 0.0;
                        ySigma_in = 0.0;
                        
                        CFL_tau[icube][i][j][k] = -DTau_CFL;
                    
                    
                    } else {
                        
                        yV_in_1 = 0.0;
                        ySigma_in = 0.0;
                    
                        yV_out_1 = 0.0;
                        ySigma_out = 0.0;
                        
                    }
                    
                    
                    
                    // C_plan = 0.5*sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);
                    
                    if( ZZ < abs5 ) {
                        
                        zV_in_1 = ((abs5-ZZ)/lenght_absZm)*((abs5-ZZ)/lenght_absZm)*((abs5-ZZ)/lenght_absZm);
                        zSigma_in = zV_in_1*zSigma_in_0/Char_D;
                        
                        zV_in_1 = zV_in_1*zV_in_0*C_plan;
                        
                        zV_out_1 = 0.0;
                        zSigma_out = 0.0;
                        
                        CFL_tau[icube][i][j][k] = -DTau_CFL;
                        
                        // if Z is slip or no slip zV_in_1 et al. should be set to 0//
                        
					    
						} 
					else if( ZZ > abs6  ) {
                    
                        zV_out_1 = ((ZZ-abs6)/lenght_absZp)*((ZZ-abs6)/lenght_absZp)*((ZZ-abs6)/lenght_absZp);
                        zSigma_out = zV_out_1*zSigma_out_0/Char_D;
                        
                        zV_out_1 = zV_out_1*zV_out_0*C_plan;
                        
                        zV_in_1 = 0.0;
                        zSigma_in = 0.0;
                        
                        CFL_tau[icube][i][j][k] = -DTau_CFL;
                    
                    } else {
                        
                        zV_in_1 = 0.0;
                        zSigma_in = 0.0;
                    
                        zV_out_1 = 0.0;
                        zSigma_out = 0.0;
                        
                    }
                    
                    
                    sr1 = sr2 = sr3 = sr4 = sr5 = 0.0;
                    sl1 = sl2 = sl3 = sl4 = sl5 = 0.0;
                    
                    U1 = rho-rho0;
                    U2 = rho*U-rho0*U0;
                    U3 = rho*V-rho0*V0;
                    U4 = rho*W-rho0*W0;
                    U5 = E-E0;
                    
                    
                    sl1 = xV_in_1*( rho  - U1_[icube][i-1][j][k][0])/dx + xSigma_in*U1; 
                    sl2 = xV_in_1*( rho*U- U1_[icube][i-1][j][k][1])/dx + xSigma_in*U2;
                    sl3 = xV_in_1*( rho*V- U1_[icube][i-1][j][k][2])/dx + xSigma_in*U3;
                    sl4 = xV_in_1*( rho*W- U1_[icube][i-1][j][k][3])/dx + xSigma_in*U4;
                    sl5 = xV_in_1*( E    - U1_[icube][i-1][j][k][4])/dx + xSigma_in*U5;
                    
                    
                    sl1 = sl1 + yV_in_1*( rho  - U1_[icube][i][j+1][k][0])/dy + ySigma_in*U1; 
                    sl2 = sl2 + yV_in_1*( rho*U- U1_[icube][i][j+1][k][1])/dy + ySigma_in*U2;
                    sl3 = sl3 + yV_in_1*( rho*V- U1_[icube][i][j+1][k][2])/dy + ySigma_in*U3;
                    sl4 = sl4 + yV_in_1*( rho*W- U1_[icube][i][j+1][k][3])/dy + ySigma_in*U4;
                    sl5 = sl5 + yV_in_1*( E    - U1_[icube][i][j+1][k][4])/dy + ySigma_in*U5;
                    
                    sl1 = sl1 + zV_in_1*( rho  - U1_[icube][i][j][k+1][0])/dz + zSigma_in*U1; 
                    sl2 = sl2 + zV_in_1*( rho*U- U1_[icube][i][j][k+1][1])/dz + zSigma_in*U2;
                    sl3 = sl3 + zV_in_1*( rho*V- U1_[icube][i][j][k+1][2])/dz + zSigma_in*U3;
                    sl4 = sl4 + zV_in_1*( rho*W- U1_[icube][i][j][k+1][3])/dz + zSigma_in*U4;
                    sl5 = sl5 + zV_in_1*( E    - U1_[icube][i][j][k+1][4])/dz + zSigma_in*U5;
                    
                    
                    
                    sr1 = xV_out_1*( rho  - U1_[icube][i-1][j][k][0])/dx + xSigma_out*U1; 
                    sr2 = xV_out_1*( rho*U- U1_[icube][i-1][j][k][1])/dx + xSigma_out*U2;
                    sr3 = xV_out_1*( rho*V- U1_[icube][i-1][j][k][2])/dx + xSigma_out*U3;
                    sr4 = xV_out_1*( rho*W- U1_[icube][i-1][j][k][3])/dx + xSigma_out*U4;
                    sr5 = xV_out_1*( E    - U1_[icube][i-1][j][k][4])/dx + xSigma_out*U5;
                    
                    
                    sr1 = sr1 + yV_out_1*( rho  - U1_[icube][i][j-1][k][0])/dy + ySigma_out*U1; 
                    sr2 = sr2 + yV_out_1*( rho*U- U1_[icube][i][j-1][k][1])/dy + ySigma_out*U2;
                    sr3 = sr3 + yV_out_1*( rho*V- U1_[icube][i][j-1][k][2])/dy + ySigma_out*U3;
                    sr4 = sr4 + yV_out_1*( rho*W- U1_[icube][i][j-1][k][3])/dy + ySigma_out*U4;
                    sr5 = sr5 + yV_out_1*( E    - U1_[icube][i][j-1][k][4])/dy + ySigma_out*U5;
                    
                    sr1 = sr1 + zV_out_1*( rho  - U1_[icube][i][j][k-1][0])/dz + zSigma_out*U1; 
                    sr2 = sr2 + zV_out_1*( rho*U- U1_[icube][i][j][k-1][1])/dz + zSigma_out*U2;
                    sr3 = sr3 + zV_out_1*( rho*V- U1_[icube][i][j][k-1][2])/dz + zSigma_out*U3;
                    sr4 = sr4 + zV_out_1*( rho*W- U1_[icube][i][j][k-1][3])/dz + zSigma_out*U4;
                    sr5 = sr5 + zV_out_1*( E    - U1_[icube][i][j][k-1][4])/dz + zSigma_out*U5;
                    
                    
                    Fabs[icube][i][j][k][0] = -sr1-sl1;
                    Fabs[icube][i][j][k][1] = -sr2-sl2;
                    Fabs[icube][i][j][k][2] = -sr3-sl3;
                    Fabs[icube][i][j][k][3] = -sr4-sl4;
                    Fabs[icube][i][j][k][4] = -sr5-sl5;
                    

                }
            }
         }

      }
      
      
      

}