



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#include "Resolution.h"

extern int Ncube;  
extern int NYbc_l;
extern int NYbc_u;


void BCM_Y_boundary_condition
(
// ================================================================================ //
int nYbc_l,
int nYbc_u,

int (*Ybc_l) = new int[NYbc_l+1],
int (*Ybc_u) = new int[NYbc_u+1],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
)

{

#include "BCM.h"
#include "prm.h"


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