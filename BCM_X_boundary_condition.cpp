



#include <mpi.h>
#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>


#include "Resolution.h"

extern int Ncube;  
extern int NXbc_l;
extern int NXbc_u;


void BCM_X_boundary_condition
(
// ================================================================================ //
int myid,

int nXbc_l,
int nXbc_u,

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
)

{

#include "BCM.h"
#include "prm.h"

	double rho,U,V,W,VV,P,C,T,h,H;
	


	for (icube = 1; icube <= nXbc_l; icube++) {  
		
		iicube = Xbc_l[icube];

#pragma omp parallel for private(k,rho,U,V,W,VV,P)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {  

				//rho = U1_[iicube][2][j][k][0];
				//U = U1_[iicube][2][j][k][1]/rho;
				//V = U1_[iicube][2][j][k][2]/rho;
				//W = U1_[iicube][2][j][k][3]/rho;
				//VV = U*U+V*V+W*W;
				//P = (U1_[iicube][2][j][k][4]-0.5*rho*VV)*(K-1);

				////rho = rho0;
				//rho = rho;
				//U = Uin;
				//V = V0;
				//W = W0;
				//VV = U*U+V*V+W*W;


				U1_[iicube][1][j][k][0] = U1_[iicube][2][j][k][0];
				U1_[iicube][1][j][k][1] = U1_[iicube][2][j][k][1];
				U1_[iicube][1][j][k][2] = U1_[iicube][2][j][k][2];
				U1_[iicube][1][j][k][3] = U1_[iicube][2][j][k][3];
				U1_[iicube][1][j][k][4] = U1_[iicube][2][j][k][4];

				U1_[iicube][0][j][k][0] = U1_[iicube][2][j][k][0];
				U1_[iicube][0][j][k][1] = U1_[iicube][2][j][k][1];
				U1_[iicube][0][j][k][2] = U1_[iicube][2][j][k][2];
				U1_[iicube][0][j][k][3] = U1_[iicube][2][j][k][3];
				U1_[iicube][0][j][k][4] = U1_[iicube][2][j][k][4];

			}
		}
	}			




	for (icube = 1; icube <= nXbc_u; icube++) {  

		iicube = Xbc_u[icube];

#pragma omp parallel for private(k,rho,U,V,W,VV)

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
				//U1_[iicube][nxx][j][k][4] = P/(K-1)+0.5*rho*VV;
				U1_[iicube][nxx][j][k][4] = U1_[iicube][nx][j][k][4];

				U1_[iicube][nxxx][j][k][0] = U1_[iicube][nx][j][k][0];
				U1_[iicube][nxxx][j][k][1] = U1_[iicube][nx][j][k][1];
				U1_[iicube][nxxx][j][k][2] = U1_[iicube][nx][j][k][2];
				U1_[iicube][nxxx][j][k][3] = U1_[iicube][nx][j][k][3];
				//U1_[iicube][nxxx][j][k][4] = P/(K-1)+0.5*rho*VV;
				U1_[iicube][nxxx][j][k][4] = U1_[iicube][nx][j][k][4];

			}
		}

	}



}