




#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "Resolution.h"

extern int Ncube;    
extern int NBC_minus;

void BCM_Ghostcell_minus
(
// =============================================================================== //
int myid,

int *NNBC,

double Th,

double (*weight) = new double[NBC_minus*8+1],
int (*GCindex) = new int[NBC_minus*4+1],
int (*IPsur) = new int[NBC_minus*4+1],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// =============================================================================== //
)

{
	
#include "BCM.h"
#include "prm.h"
	
int gicube,gi,gj,gk;
int iNBC, Ntemp, ii, jj, kk, iw, i_switch;

double rho, U, V, W, VV, P, T, temp;

double u1[8], u2[8], u3[8], u4[8], u5[8];

double w[8];

double tmp1,tmp2,tmp3,tmp4;

double er_p = 0.9999;


// ---- GCindex => IPsur => weight ---- //
//
//#pragma omp parallel for private(Ntemp,icube,i,j,k,\
//	ii, jj, kk, iw, i_switch, \
//	gicube,gi,gj,gk,\
//	u1,u2,u3,u4,u5,\
//	rho,P,U,V,W,VV,temp\
//	)
//	

//#pragma omp parallel 


	for (iNBC = 1; iNBC <= *NNBC; iNBC++) {

		Ntemp = (iNBC-1)*4;

		icube = IPsur[Ntemp+1];
		i = IPsur[Ntemp+2];
		j = IPsur[Ntemp+3];
		k = IPsur[Ntemp+4];


		Ntemp = (iNBC-1)*8;


		for( iw = 0; iw < 8; iw ++ ) w[iw] = weight[Ntemp+1+iw];

		
		tmp1 = w[1];
		tmp2 = w[3];
		tmp3 = w[6];
		tmp4 = w[4];

		w[1] = tmp4;
		w[3] = tmp1;
		w[6] = tmp2;
		w[4] = tmp3;


		for (ii = i; ii <= i+1; ii++) {
			for (jj = j; jj <= j+1; jj++) {
				for (kk = k; kk <= k+1; kk++) {  

					rho = U1_[icube][ii][jj][kk][0];
					U = U1_[icube][ii][jj][kk][1]/rho;
					V = U1_[icube][ii][jj][kk][2]/rho;
					W = U1_[icube][ii][jj][kk][3]/rho;
					P = ( U1_[icube][ii][jj][kk][4]-0.5*rho*(U*U+V*V+W*W) )*(K-1);
					
					u1[(ii-i)*4 + 2*(jj-j) + (kk-k)] = rho;
					u2[(ii-i)*4 + 2*(jj-j) + (kk-k)] = U;
					u3[(ii-i)*4 + 2*(jj-j) + (kk-k)] = V;
					u4[(ii-i)*4 + 2*(jj-j) + (kk-k)] = W;
					u5[(ii-i)*4 + 2*(jj-j) + (kk-k)] = P;
					
				}
			}
		}

		
		
		Ntemp = (iNBC-1)*4;

		gicube = GCindex[Ntemp+1];
		gi = GCindex[Ntemp+2];
		gj = GCindex[Ntemp+3];
		gk = GCindex[Ntemp+4];


		rho = 0.0;
		U = 0.0;
		V = 0.0;
		W = 0.0;
		P = 0.0;

		i_switch = 0;


		
		for (ii = i; ii <= i+1; ii++) {
			for (jj = j; jj <= j+1; jj++) {
				for (kk = k; kk <= k+1; kk++) {  

					if( gi == ii && gj == jj && gk == kk ) {

						i_switch = 1;
						goto ourfor;
					
					}

				}
			}
		}
		
		ourfor:

		if( i_switch == 1 ) {

			iw = (gi-i)*4 + 2*(gj-j) + (gk-k);


			if (w[iw] > er_p) {

				rho = u1[iw] - rho0;
				P = u5[iw] - P0; 

			}
			else {

				for( ii = 0; ii < 8; ii ++ ) {

					if ( ii == iw ) continue;

					rho = rho + ( u1[ii] - rho0 ) * w[ii];
					P = P + ( u5[ii] - P0 ) * w[ii];

				}

				rho = rho/( 1.0-w[iw] );
				P = P/( 1.0-w[iw] );

			}



			for( ii = 0; ii < 8; ii ++ ) {

				if ( ii == iw ) continue;

				U = U + u2[ii]*w[ii];
				V = V + u3[ii]*w[ii];
				W = W + u4[ii]*w[ii];

			}

			rho = rho + rho0;
			P = P + P0;
			U = U / ( 2.0-w[iw] );
			V = V / ( 2.0-w[iw] );
			W = W / ( 2.0-w[iw] );

		}

		else {

			for( ii = 0; ii < 8; ii ++ ) {

				rho = rho + ( u1[ii] - rho0 ) * w[ii];
				U = U + u2[ii]*w[ii];
				V = V + u3[ii]*w[ii];
				W = W + u4[ii]*w[ii];
				P = P + ( u5[ii] - P0 ) * w[ii];

			}

			rho = rho + rho0;
			P = P + P0;
			U = 0.5*U;
			V = 0.5*V;
			W = 0.5*W;

		}    // ---- if( i_switch == 1 ) ---- //


		VV = U*U+V*V+W*W;

		U1_[gicube][gi][gj][gk][0] = rho;
		U1_[gicube][gi][gj][gk][1] = rho*U;
		U1_[gicube][gi][gj][gk][2] = rho*V;
		U1_[gicube][gi][gj][gk][3] = rho*W;
		U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

	}

	
}