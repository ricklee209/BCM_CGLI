#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"
#include "Pre_selection.h"

extern int Ncube;   

void BCM_ADM_filter
(
// ================================================================================ //

	int myid,

	int ncube,

	double (*cube_size) = new double[Ncube],

	double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*Roe_dis)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],

	double (*filter)[11][X_size][Y_size][Z_size][Ndim] = new double[Ncube][11][X_size][Y_size][Z_size][Ndim]


// ================================================================================ //
)

{


#include "BCM.h"
#include "prm.h"


	int ic;

	double tmp;


	double u1,u2,u3,u4,u5,

		u1i1,u2i1,u3i1,u4i1,
		u1i2,u2i2,u3i2,u4i2,
		u1i3,u2i3,u3i3,u4i3,

		u1_i1,u2_i1,u3_i1,u4_i1,
		u1_i2,u2_i2,u3_i2,u4_i2,
		u1_i3,u2_i3,u3_i3,u4_i3,

		u1j1,u2j1,u3j1,u4j1,
		u1j2,u2j2,u3j2,u4j2,
		u1j3,u2j3,u3j3,u4j3,

		u1_j1,u2_j1,u3_j1,u4_j1,
		u1_j2,u2_j2,u3_j2,u4_j2,
		u1_j3,u2_j3,u3_j3,u4_j3,

		u1k1,u2k1,u3k1,u4k1,
		u1k2,u2k2,u3k2,u4k2,
		u1k3,u2k3,u3k3,u4k3,

		u1_k1,u2_k1,u3_k1,u4_k1,
		u1_k2,u2_k2,u3_k2,u4_k2,
		u1_k3,u2_k3,u3_k3,u4_k3;

	double delta_I = 0.01;
  double ration_dissp = 0.7;
  double tolerance = 0.03; // 3 % 


	for(ic = 0; ic <= 10; ic++) {

		#pragma omp parallel for private(i,j,k)

		for (icube = 1; icube < ncube; icube++) {  

			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
					for (k = 0; k <= nzzz; k++) {


						filter[icube][ic][i][j][k][0] = U1_[icube][i][j][k][0];
						filter[icube][ic][i][j][k][1] = U1_[icube][i][j][k][1];
						filter[icube][ic][i][j][k][2] = U1_[icube][i][j][k][2];
						filter[icube][ic][i][j][k][3] = U1_[icube][i][j][k][3];
						filter[icube][ic][i][j][k][4] = U1_[icube][i][j][k][4];


					}
				}
			}

		}

	}




	for(ic = 4; ic <= 9; ic++) {

#pragma omp parallel for private(i,j,k,\
		u1,u2,u3,u4,\
		u1i1,u2i1,u3i1,u4i1,\
		u1_i1,u2_i1,u3_i1,u4_i1,\
		u1j1,u2j1,u3j1,u4j1,\
		u1_j1,u2_j1,u3_j1,u4_j1,\
		u1k1,u2k1,u3k1,u4k1,\
		u1_k1,u2_k1,u3_k1,u4_k1\
		)

		for (icube = 1; icube < ncube; icube++) {  

			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {

						u1 = filter[icube][ic][i][j][k][0];
						u2 = filter[icube][ic][i][j][k][1]/u1;
						u3 = filter[icube][ic][i][j][k][2]/u1;
						u4 = filter[icube][ic][i][j][k][3]/u1;


						u1_k1 = filter[icube][ic][i][j][k-1][0];
						u2_k1 = filter[icube][ic][i][j][k-1][1]/u1_k1;
						u3_k1 = filter[icube][ic][i][j][k-1][2]/u1_k1;
						u4_k1 = filter[icube][ic][i][j][k-1][3]/u1_k1;


						u1k1 = filter[icube][ic][i][j][k+1][0];
						u2k1 = filter[icube][ic][i][j][k+1][1]/u1k1;
						u3k1 = filter[icube][ic][i][j][k+1][2]/u1k1;
						u4k1 = filter[icube][ic][i][j][k+1][3]/u1k1;


						filter[icube][3][i][j][k][1] = 0.25*(u2_k1+2.0*u2+u2k1);
						filter[icube][3][i][j][k][2] = 0.25*(u3_k1+2.0*u3+u3k1);
						filter[icube][3][i][j][k][3] = 0.25*(u4_k1+2.0*u4+u4k1);
            
            
					}
				}
			}

      
			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {



						u1 = filter[icube][3][i][j][k][0];
						u2 = filter[icube][3][i][j][k][1]/u1;
						u3 = filter[icube][3][i][j][k][2]/u1;
						u4 = filter[icube][3][i][j][k][3]/u1;


						u1_j1 = filter[icube][3][i][j-1][k][0];
						u2_j1 = filter[icube][3][i][j-1][k][1]/u1_j1;
						u3_j1 = filter[icube][3][i][j-1][k][2]/u1_j1;
						u4_j1 = filter[icube][3][i][j-1][k][3]/u1_j1;

            
						u1j1 = filter[icube][3][i][j+1][k][0];
						u2j1 = filter[icube][3][i][j+1][k][1]/u1j1;
						u3j1 = filter[icube][3][i][j+1][k][2]/u1j1;
						u4j1 = filter[icube][3][i][j+1][k][3]/u1j1;


						filter[icube][2][i][j][k][1] = 0.25*(u2_j1+2.0*u2+u2j1);
						filter[icube][2][i][j][k][2] = 0.25*(u3_j1+2.0*u3+u3j1);
						filter[icube][2][i][j][k][3] = 0.25*(u4_j1+2.0*u4+u4j1);
            
            
					}
				}
			}

      
			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {



						u1 = filter[icube][2][i][j][k][0];
						u2 = filter[icube][2][i][j][k][1]/u1;
						u3 = filter[icube][2][i][j][k][2]/u1;
						u4 = filter[icube][2][i][j][k][3]/u1;


						u1_i1 = filter[icube][2][i-1][j][k][0];
						u2_i1 = filter[icube][2][i-1][j][k][1]/u1_i1;
						u3_i1 = filter[icube][2][i-1][j][k][2]/u1_i1;
						u4_i1 = filter[icube][2][i-1][j][k][3]/u1_i1;


						u1i1 = filter[icube][2][i+1][j][k][0];
						u2i1 = filter[icube][2][i+1][j][k][1]/u1i1;
						u3i1 = filter[icube][2][i+1][j][k][2]/u1i1;
						u4i1 = filter[icube][2][i+1][j][k][3]/u1i1;


						filter[icube][ic+1][i][j][k][1] = 0.25*(u2_i1+2.0*u2+u2i1);
						filter[icube][ic+1][i][j][k][2] = 0.25*(u3_i1+2.0*u3+u3i1);
						filter[icube][ic+1][i][j][k][3] = 0.25*(u4_i1+2.0*u4+u4i1);


					}
				}
			}

		}

	}



	#pragma omp parallel for private(i,j,k)

		for (icube = 1; icube < ncube; icube++) {  

			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {

						filter[icube][0][i][j][k][1] =   6.0 *filter[icube][4][i][j][k][1] \
														-15.0*filter[icube][5][i][j][k][1] \
														+20.0*filter[icube][6][i][j][k][1] \
														-15.0*filter[icube][7][i][j][k][1] \
														+6.0 *filter[icube][8][i][j][k][1] \
														-1.0 *filter[icube][9][i][j][k][1];

						filter[icube][0][i][j][k][2] =   6.0 *filter[icube][4][i][j][k][2] \
														-15.0*filter[icube][5][i][j][k][2] \
														+20.0*filter[icube][6][i][j][k][2] \
														-15.0*filter[icube][7][i][j][k][2] \
														+6.0 *filter[icube][8][i][j][k][2] \
														-1.0 *filter[icube][9][i][j][k][2];

						filter[icube][0][i][j][k][3] =   6.0 *filter[icube][4][i][j][k][3] \
														-15.0*filter[icube][5][i][j][k][3] \
														+20.0*filter[icube][6][i][j][k][3] \
														-15.0*filter[icube][7][i][j][k][3] \
														+6.0 *filter[icube][8][i][j][k][3] \
														-1.0 *filter[icube][9][i][j][k][3];



					}
				}
			}
		
		}





		

	for(ic = 4; ic <= 7; ic++) {

#pragma omp parallel for private(i,j,k,\
		u1,u2,u3,u4,\
		u1i1,u2i1,u3i1,u4i1,\
		u1_i1,u2_i1,u3_i1,u4_i1,\
		u1j1,u2j1,u3j1,u4j1,\
		u1_j1,u2_j1,u3_j1,u4_j1,\
		u1k1,u2k1,u3k1,u4k1,\
		u1_k1,u2_k1,u3_k1,u4_k1\
		)

		for (icube = 1; icube < ncube; icube++) {  

			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {
						
						u1 = filter[icube][ic][i][j][k][0];
						u2 = filter[icube][ic][i][j][k][1]/u1;
						u3 = filter[icube][ic][i][j][k][2]/u1;
						u4 = filter[icube][ic][i][j][k][3]/u1;


						u1_k1 = filter[icube][ic][i][j][k-1][0];
						u2_k1 = filter[icube][ic][i][j][k-1][1]/u1_k1;
						u3_k1 = filter[icube][ic][i][j][k-1][2]/u1_k1;
						u4_k1 = filter[icube][ic][i][j][k-1][3]/u1_k1;


						u1k1 = filter[icube][ic][i][j][k+1][0];
						u2k1 = filter[icube][ic][i][j][k+1][1]/u1k1;
						u3k1 = filter[icube][ic][i][j][k+1][2]/u1k1;
						u4k1 = filter[icube][ic][i][j][k+1][3]/u1k1;

						filter[icube][3][i][j][k][0] = 0.25*(u1_k1+2.0*u1+u1k1);
						filter[icube][3][i][j][k][1] = 0.25*(u2_k1+2.0*u2+u2k1);
						filter[icube][3][i][j][k][2] = 0.25*(u3_k1+2.0*u3+u3k1);
						filter[icube][3][i][j][k][3] = 0.25*(u4_k1+2.0*u4+u4k1);
						

						u1 = filter[icube][3][i][j][k][0];
						u2 = filter[icube][3][i][j][k][1]/u1;
						u3 = filter[icube][3][i][j][k][2]/u1;
						u4 = filter[icube][3][i][j][k][3]/u1;


						u1_j1 = filter[icube][3][i][j-1][k][0];
						u2_j1 = filter[icube][3][i][j-1][k][1]/u1_j1;
						u3_j1 = filter[icube][3][i][j-1][k][2]/u1_j1;
						u4_j1 = filter[icube][3][i][j-1][k][3]/u1_j1;


						u1j1 = filter[icube][3][i][j+1][k][0];
						u2j1 = filter[icube][3][i][j+1][k][1]/u1j1;
						u3j1 = filter[icube][3][i][j+1][k][2]/u1j1;
						u4j1 = filter[icube][3][i][j+1][k][3]/u1j1;



						filter[icube][2][i][j][k][0] = 0.25*(u1_j1+2.0*u1+u1j1);
						filter[icube][2][i][j][k][1] = 0.25*(u2_j1+2.0*u2+u2j1);
						filter[icube][2][i][j][k][2] = 0.25*(u3_j1+2.0*u3+u3j1);
						filter[icube][2][i][j][k][3] = 0.25*(u4_j1+2.0*u4+u4j1);

						
						u1 = filter[icube][2][i][j][k][0];
						u2 = filter[icube][2][i][j][k][1]/u1;
						u3 = filter[icube][2][i][j][k][2]/u1;
						u4 = filter[icube][2][i][j][k][3]/u1;


						u1_i1 = filter[icube][2][i-1][j][k][0];
						u2_i1 = filter[icube][2][i-1][j][k][1]/u1_i1;
						u3_i1 = filter[icube][2][i-1][j][k][2]/u1_i1;
						u4_i1 = filter[icube][2][i-1][j][k][3]/u1_i1;


						u1i1 = filter[icube][2][i+1][j][k][0];
						u2i1 = filter[icube][2][i+1][j][k][1]/u1i1;
						u3i1 = filter[icube][2][i+1][j][k][2]/u1i1;
						u4i1 = filter[icube][2][i+1][j][k][3]/u1i1;


						filter[icube][ic+1][i][j][k][0] = 0.25*(u1_i1+2.0*u1+u1i1);
						filter[icube][ic+1][i][j][k][1] = 0.25*(u2_i1+2.0*u2+u2i1);
						filter[icube][ic+1][i][j][k][2] = 0.25*(u3_i1+2.0*u3+u3i1);
						filter[icube][ic+1][i][j][k][3] = 0.25*(u4_i1+2.0*u4+u4i1);


					}
				}
			}

		}

	}




	
	#pragma omp parallel for private(i,j,k)

		for (icube = 1; icube < ncube; icube++) {  

			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {

						filter[icube][1][i][j][k][1] = 4.0*filter[icube][4][i][j][k][1] \
                                          -6.0*filter[icube][5][i][j][k][1] \
                                          +4.0*filter[icube][6][i][j][k][1] \
                                          -1.0*filter[icube][7][i][j][k][1];

						filter[icube][1][i][j][k][2] = 4.0*filter[icube][4][i][j][k][2] \
                                          -6.0*filter[icube][5][i][j][k][2] \
                                          +4.0*filter[icube][6][i][j][k][2] \
                                          -1.0*filter[icube][7][i][j][k][2];

						filter[icube][1][i][j][k][3] = 4.0*filter[icube][4][i][j][k][3] \
                                          -6.0*filter[icube][5][i][j][k][3] \
                                          +4.0*filter[icube][6][i][j][k][3] \
                                          -1.0*filter[icube][7][i][j][k][3];

					}
				}
			}
		
		}



		#pragma omp parallel for private(i,j,k,\
		u1,u2,u3,u4,\
		u1i1,u2i1,u3i1,\
		u1i2,u2i2,u3i2,\
		tmp)

		for (icube = 1; icube < ncube; icube++) {  

			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {
												
						u1 = U1_[icube][i][j][k][0];
						u2 = U1_[icube][i][j][k][1]/u1;
						u3 = U1_[icube][i][j][k][2]/u1;
						u4 = U1_[icube][i][j][k][3]/u1;

						u1i1 = u2-filter[icube][0][i][j][k][1];
						u2i1 = u3-filter[icube][0][i][j][k][2];
						u3i1 = u4-filter[icube][0][i][j][k][3];

						u1i2 = u2-filter[icube][1][i][j][k][1];
						u2i2 = u3-filter[icube][1][i][j][k][2];
						u3i2 = u4-filter[icube][1][i][j][k][3];


						tmp = ( u1i1*u1i1 + u2i1*u2i1 + u3i1*u3i1 )/ \
							  ( u1i2*u1i2 + u2i2*u2i2 + u3i2*u3i2 );


						if( tmp > (1.0+tolerance)*ration_dissp ) {
							
							Roe_dis[icube][i][j][k] = min(Roe_dis[icube][i][j][k]+delta_I, 1.0); 

						}

						if(tmp < (1.0-tolerance)*ration_dissp) {

							Roe_dis[icube][i][j][k] = max(Roe_dis[icube][i][j][k]-delta_I, 0.0); 

						}


					}
				}
			}
		
		}



}