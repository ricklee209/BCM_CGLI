



#include <mpi.h>
#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


#include "Resolution.h"

extern int Ncube;  
extern int NXbc_l;
extern int NXbc_u;


void BCM_Abs_X_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nXbc_l,
int nXbc_u,

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
)

{

#include "BCM.h"
#include "prm.h"

	double rho,U,V,W,VV,P,C,T,h,H,E;


	double beta, C_plan;

	double V_in_0 = 1.15;
	
	double Sigma_in_0 = 0.1;

	double V_out_0 = 1.25;
	
	double Sigma_out_0 = 1.25;
	
	double V_in_1, V_out_1;
	
	double Sigma_in, Sigma_out;

	double e0;



#pragma omp parallel for private(iicube,i,j,k,dx,dy,dz,rho,U,V,W,VV,P,C,E,e0,beta,C_plan,V_in_1,Sigma_in) schedule(dynamic)

		for (icube = 1; icube <= nXbc_l; icube++) {  

			iicube = Xbc_l[icube];

			dx = dy = dz = cube_size[iicube]/(NcubeX*1.0);

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  


						rho = U1_[iicube][i][j][k][0];
						U = U1_[iicube][i][j][k][1]/rho;
						V = U1_[iicube][i][j][k][2]/rho;
						W = U1_[iicube][i][j][k][3]/rho;
						VV = U*U+V*V+W*W;
						E = U1_[iicube][i][j][k][4];
						P = (E-0.5*rho*VV)*(K-1);
						C = K*P/rho;

						e0 = P0/(K-1)+0.5*rho0*(U0*U0+V0*V0+W0*W0);

						/* preconditioning */
						beta = max(VV/C,e);

						C_plan = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);

						V_in_1 = pow( (nx-i+0.5)/(NcubeX*1.0), 3.0 )* V_in_0*C_plan;

						Sigma_in = Sigma_in_0*V_in_1/0.000045;

						Fabs[iicube][i][j][k][0] = -( V_in_1*( rho  -U1_[iicube][i-1][j][k][0] )/dx+Sigma_in*(rho  -rho0)    );
						Fabs[iicube][i][j][k][1] = -( V_in_1*( rho*U-U1_[iicube][i-1][j][k][1] )/dx+Sigma_in*(rho*U-rho0*U0) );
						Fabs[iicube][i][j][k][2] = -( V_in_1*( rho*V-U1_[iicube][i-1][j][k][2] )/dx+Sigma_in*(rho*V-rho0*V0) );
						Fabs[iicube][i][j][k][3] = -( V_in_1*( rho*W-U1_[iicube][i-1][j][k][3] )/dx+Sigma_in*(rho*W-rho0*W0) );
						Fabs[iicube][i][j][k][4] = -( V_in_1*( E    -U1_[iicube][i-1][j][k][4] )/dx+Sigma_in*(E    -e0)      );
						 

					}
				}
			}

		}

		
	//#pragma omp barrier

	

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



	
	
#pragma omp parallel for private(iicube,i,j,k,dx,dy,dz,rho,U,V,W,VV,P,C,E,e0,beta,C_plan,V_out_1,Sigma_out) schedule(dynamic)

		for (icube = 1; icube <= nXbc_u; icube++) {  

			iicube = Xbc_u[icube];

			dx = dy = dz = cube_size[iicube]/(NcubeX*1.0);

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j <= nyyy; j++) {
					for (k = 0; k <= nzzz; k++) {  


						rho = U1_[iicube][i][j][k][0];
						U = U1_[iicube][i][j][k][1]/rho;
						V = U1_[iicube][i][j][k][2]/rho;
						W = U1_[iicube][i][j][k][3]/rho;
						VV = U*U+V*V+W*W;
						E = U1_[iicube][i][j][k][4];
						P = (E-0.5*rho*VV)*(K-1);
						C = K*P/rho;

						e0 = P0/(K-1)+0.5*rho0*(U0*U0+V0*V0+W0*W0);

						/* preconditioning */
						beta = max(VV/C,e);

						C_plan = 0.5*sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);

						V_out_1 = pow( (i-n_buffer+0.5)/(NcubeX*1.0), 3.0 )* V_out_0*C_plan;

						Sigma_out = Sigma_out_0*V_out_1/0.000045;

						
						Fabs[iicube][i][j][k][0] = -( V_out_1*( rho  -U1_[iicube][i-1][j][k][0] )/dx+Sigma_out*(rho  -rho0)    );
						Fabs[iicube][i][j][k][1] = -( V_out_1*( rho*U-U1_[iicube][i-1][j][k][1] )/dx+Sigma_out*(rho*U-rho0*U0) );
						Fabs[iicube][i][j][k][2] = -( V_out_1*( rho*V-U1_[iicube][i-1][j][k][2] )/dx+Sigma_out*(rho*V-rho0*V0) );
						Fabs[iicube][i][j][k][3] = -( V_out_1*( rho*W-U1_[iicube][i-1][j][k][3] )/dx+Sigma_out*(rho*W-rho0*W0) );
						Fabs[iicube][i][j][k][4] = -( V_out_1*( E    -U1_[iicube][i-1][j][k][4] )/dx+Sigma_out*(E    -e0)      );
					}
				}
			}

	}

		 

	
	for (icube = 1; icube <= nXbc_u; icube++) {  

		iicube = Xbc_u[icube];
		
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {  

				rho = U1_[iicube][nx][j][k][0];
				U = U1_[iicube][nx][j][k][1]/rho;
				V = U1_[iicube][nx][j][k][2]/rho;
				W = U1_[iicube][nx][j][k][3]/rho;
				VV = U*U+V*V+W*W;
				P = 101300;

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







}
