



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


void BCM_Z_boundary_condition
(
// ================================================================================ //
int nZbc_l,
int nZbc_u,

int (*Zbc_l) = new int[NZbc_l+1],
int (*Zbc_u) = new int[NZbc_u+1],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
)

{

#include "BCM.h"
#include "prm.h"

	double rho,U,V,W,VV,P,C,T,h,H;
	


#pragma omp parallel for private(iicube,i,j)
	for (icube = 1; icube <= nZbc_l; icube++) {  

		iicube = Zbc_l[icube];

			for (i = 2; i <= nx; i++) {
				for (j = 2; j <= ny; j++) {  

					U1_[iicube][i][j][1][0] = U1_[iicube][i][j][2][0];
					U1_[iicube][i][j][1][1] = U1_[iicube][i][j][2][1];
					U1_[iicube][i][j][1][2] = U1_[iicube][i][j][2][2];
					U1_[iicube][i][j][1][3] = U1_[iicube][i][j][2][3];
					U1_[iicube][i][j][1][4] = U1_[iicube][i][j][2][4];

					U1_[iicube][i][j][0][0] = U1_[iicube][i][j][3][0];
					U1_[iicube][i][j][0][1] = U1_[iicube][i][j][3][1];
					U1_[iicube][i][j][0][2] = U1_[iicube][i][j][3][2];
					U1_[iicube][i][j][0][3] = U1_[iicube][i][j][3][3];
					U1_[iicube][i][j][0][4] = U1_[iicube][i][j][3][4];

				}
			}

	}			



#pragma omp parallel for private(iicube,i,j)
	for (icube = 1; icube <= nZbc_u; icube++) {  

		iicube = Zbc_u[icube];

			for (i = 2; i <= nx; i++) {
				for (j = 2; j <= ny; j++) {  

					U1_[iicube][i][j][nzz][0] = U1_[iicube][i][j][nz][0];
					U1_[iicube][i][j][nzz][1] = U1_[iicube][i][j][nz][1];
					U1_[iicube][i][j][nzz][2] = U1_[iicube][i][j][nz][2];
					U1_[iicube][i][j][nzz][3] = U1_[iicube][i][j][nz][3];
					U1_[iicube][i][j][nzz][4] = U1_[iicube][i][j][nz][4];

					U1_[iicube][i][j][nzzz][0] = U1_[iicube][i][j][nz-1][0];
					U1_[iicube][i][j][nzzz][1] = U1_[iicube][i][j][nz-1][1];
					U1_[iicube][i][j][nzzz][2] = U1_[iicube][i][j][nz-1][2];
					U1_[iicube][i][j][nzzz][3] = U1_[iicube][i][j][nz-1][3];
					U1_[iicube][i][j][nzzz][4] = U1_[iicube][i][j][nz-1][4];

				}
			}
	}



}