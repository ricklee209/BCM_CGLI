



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;  
extern int NYbc_l;
extern int NYbc_u;


void BCM_Abs_Y_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nYbc_l,
int nYbc_u,

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


	double beta, C_plan;

	double V_in_0 = 0;
	
	double Sigma_in_0 = 1.25;

	double V_out_0 = 0;
	
	double Sigma_out_0 = -1.25;
	
	double V_in_1, V_out_1;
	
	double Sigma_in, Sigma_out;

	double e0;




#pragma omp parallel for private(iicube,i,j,k,dx,dy,dz,rho,U,V,W,VV,P,C,E,e0,beta,C_plan,V_in_1,Sigma_in) schedule(dynamic)

		for (icube = 1; icube <= nYbc_l; icube++) {  

			iicube = Ybc_l[icube];

			dx = dy = dz = cube_size[iicube]/(NcubeY*1.0);

			for (i = 0; i <= nxxx; i++) {
				for (j = n_buffer; j <= ny; j++) {
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

						C_plan = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);

						V_in_1 = pow( (ny-j+0.5)/(NcubeY*1.0), 3.0 )* V_in_0*C_plan;

						Sigma_in = Sigma_in_0*V_in_1/Char_D;

						Fabs[iicube][i][j][k][0] = -( V_in_1*( rho  -U1_[iicube][i][j-1][k][0] )/dy+Sigma_in*(rho  -rho0   ) );
						Fabs[iicube][i][j][k][1] = -( V_in_1*( rho*U-U1_[iicube][i][j-1][k][1] )/dy+Sigma_in*(rho*U-rho0*U0) );
						Fabs[iicube][i][j][k][2] = -( V_in_1*( rho*V-U1_[iicube][i][j-1][k][2] )/dy+Sigma_in*(rho*V-rho0*V0) );
						Fabs[iicube][i][j][k][3] = -( V_in_1*( rho*W-U1_[iicube][i][j-1][k][3] )/dy+Sigma_in*(rho*W-rho0*W0) );
						Fabs[iicube][i][j][k][4] = -( V_in_1*( E    -U1_[iicube][i][j-1][k][4] )/dy+Sigma_in*(E    -e0     ) );


					}
				}
			}

		}




	for (icube = 1; icube <= nYbc_l; icube++) {  

		iicube = Ybc_l[icube];

#pragma omp parallel for private(k)
			for (i = 0; i <= nxxx; i++) {
				for (k = n_buffer; k <= nz; k++) {  

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






	
	
#pragma omp parallel for private(iicube,i,j,k,dx,dy,dz,rho,U,V,W,VV,P,C,E,e0,beta,C_plan,V_out_1,Sigma_out) schedule(dynamic)

		for (icube = 1; icube <= nYbc_u; icube++) {  

			iicube = Ybc_u[icube];

			dx = dy = dz = cube_size[iicube]/(NcubeY*1.0);

			for (i = 0; i <= nxxx; i++) {
				for (j = n_buffer; j <= ny; j++) {
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

						C_plan = 0.5*sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);

						V_out_1 = pow( (j-n_buffer+0.5)/(NcubeY*1.0), 3.0 )* V_out_0*C_plan;

						Sigma_out = Sigma_out_0*V_out_1/Char_D;

						
						Fabs[iicube][i][j][k][0] = -( V_out_1*( rho  -U1_[iicube][i][j-1][k][0] )/dy+Sigma_out*(rho  -rho0)    );
						Fabs[iicube][i][j][k][1] = -( V_out_1*( rho*U-U1_[iicube][i][j-1][k][1] )/dy+Sigma_out*(rho*U-rho0*U0) );
						Fabs[iicube][i][j][k][2] = -( V_out_1*( rho*V-U1_[iicube][i][j-1][k][2] )/dy+Sigma_out*(rho*V-rho0*V0) );
						Fabs[iicube][i][j][k][3] = -( V_out_1*( rho*W-U1_[iicube][i][j-1][k][3] )/dy+Sigma_out*(rho*W-rho0*W0) );
						Fabs[iicube][i][j][k][4] = -( V_out_1*( E    -U1_[iicube][i][j-1][k][4] )/dy+Sigma_out*(E    -e0)      );

					}
				}
			}

	}



	for (icube = 1; icube <= nYbc_u; icube++) {  

		iicube = Ybc_u[icube];

#pragma omp parallel for private(k)
		for (i = 0; i <= nxxx; i++) {
			for (k = n_buffer; k <= nz; k++) {  

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
	


}