



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;  
extern int NZbc_l;
extern int NZbc_u;


void BCM_Abs_Z_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nZbc_l,
int nZbc_u,

int (*Zbc_l) = new int[NZbc_l+1],
int (*Zbc_u) = new int[NZbc_u+1],


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


	double beta, C_plan;

	double V_in_0 = 0;
	
	double Sigma_in_0 = 1.25;

	double V_out_0 = 0;
	
	double Sigma_out_0 = -1.25;
	
	double V_in_1, V_out_1;
	
	double Sigma_in, Sigma_out;

	double e0;




#pragma omp parallel for private(iicube,i,j,k,dx,dy,dz,rho,U,V,W,VV,P,C,E,e0,beta,C_plan,V_in_1,Sigma_in) schedule(dynamic)

		for (icube = 1; icube <= nZbc_l; icube++) {  

			iicube = Zbc_l[icube];

			dx = dy = dz = cube_size[iicube]/(NcubeZ*1.0);

			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
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

						C_plan = 0.5*sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);

						V_in_1 = pow( (nz-k+0.5)/(NcubeZ*1.0), 3.0 )* V_in_0*C_plan;

						Sigma_in = Sigma_in_0*V_in_1/0.000045;

						Fabs[iicube][i][j][k][0] = -( V_in_1*( rho  -U1_[iicube][i][j][k-1][0] )/dy+Sigma_in*(rho  -rho0   ) );
						Fabs[iicube][i][j][k][1] = -( V_in_1*( rho*U-U1_[iicube][i][j][k-1][1] )/dy+Sigma_in*(rho*U-rho0*U0) );
						Fabs[iicube][i][j][k][2] = -( V_in_1*( rho*V-U1_[iicube][i][j][k-1][2] )/dy+Sigma_in*(rho*V-rho0*V0) );
						Fabs[iicube][i][j][k][3] = -( V_in_1*( rho*W-U1_[iicube][i][j][k-1][3] )/dy+Sigma_in*(rho*W-rho0*W0) );
						Fabs[iicube][i][j][k][4] = -( V_in_1*( E    -U1_[iicube][i][j][k-1][4] )/dy+Sigma_in*(E    -e0     ) );


					}
				}
			}

		}



	for (icube = 1; icube <= nZbc_l; icube++) {  

		iicube = Zbc_l[icube];

#pragma omp parallel for private(j)
		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {  

				U1_[iicube][i][j][1][0] = U1_[iicube][i][j][2][0];
				U1_[iicube][i][j][1][1] = U1_[iicube][i][j][2][1];
				U1_[iicube][i][j][1][2] = U1_[iicube][i][j][2][2];
				U1_[iicube][i][j][1][3] = U1_[iicube][i][j][2][3];
				U1_[iicube][i][j][1][4] = U1_[iicube][i][j][2][4];

				U1_[iicube][i][j][0][0] = U1_[iicube][i][j][2][0];
				U1_[iicube][i][j][0][1] = U1_[iicube][i][j][2][1];
				U1_[iicube][i][j][0][2] = U1_[iicube][i][j][2][2];
				U1_[iicube][i][j][0][3] = U1_[iicube][i][j][2][3];
				U1_[iicube][i][j][0][4] = U1_[iicube][i][j][2][4];

			}
		}

	}			




	
	
#pragma omp parallel for private(iicube,i,j,k,dx,dy,dz,rho,U,V,W,VV,P,C,E,e0,beta,C_plan,V_out_1,Sigma_out) schedule(dynamic)

		for (icube = 1; icube <= nZbc_u; icube++) {  

			iicube = Zbc_u[icube];

			dx = dy = dz = cube_size[iicube]/(NcubeZ*1.0);

			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
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

						C_plan = 0.5*sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);

						V_out_1 = pow( (k-n_buffer+0.5)/(NcubeZ*1.0), 3.0 )* V_out_0*C_plan;

						Sigma_out = Sigma_out_0*V_out_1/Char_D;

						
						Fabs[iicube][i][j][k][0] = -( V_out_1*( rho  -U1_[iicube][i][j][k-1][0] )/dy+Sigma_out*(rho  -rho0)    );
						Fabs[iicube][i][j][k][1] = -( V_out_1*( rho*U-U1_[iicube][i][j][k-1][1] )/dy+Sigma_out*(rho*U-rho0*U0) );
						Fabs[iicube][i][j][k][2] = -( V_out_1*( rho*V-U1_[iicube][i][j][k-1][2] )/dy+Sigma_out*(rho*V-rho0*V0) );
						Fabs[iicube][i][j][k][3] = -( V_out_1*( rho*W-U1_[iicube][i][j][k-1][3] )/dy+Sigma_out*(rho*W-rho0*W0) );
						Fabs[iicube][i][j][k][4] = -( V_out_1*( E    -U1_[iicube][i][j][k-1][4] )/dy+Sigma_out*(E    -e0)      );

					}
				}
			}

	}



		

	for (icube = 1; icube <= nZbc_u; icube++) {  

		iicube = Zbc_u[icube];

#pragma omp parallel for private(j)
		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {  

				U1_[iicube][i][j][nzz][0] = U1_[iicube][i][j][nz][0];
				U1_[iicube][i][j][nzz][1] = U1_[iicube][i][j][nz][1];
				U1_[iicube][i][j][nzz][2] = U1_[iicube][i][j][nz][2];
				U1_[iicube][i][j][nzz][3] = U1_[iicube][i][j][nz][3];
				U1_[iicube][i][j][nzz][4] = U1_[iicube][i][j][nz][4];

				U1_[iicube][i][j][nzzz][0] = U1_[iicube][i][j][nz][0];
				U1_[iicube][i][j][nzzz][1] = U1_[iicube][i][j][nz][1];
				U1_[iicube][i][j][nzzz][2] = U1_[iicube][i][j][nz][2];
				U1_[iicube][i][j][nzzz][3] = U1_[iicube][i][j][nz][3];
				U1_[iicube][i][j][nzzz][4] = U1_[iicube][i][j][nz][4];

			}
		}
	}





}