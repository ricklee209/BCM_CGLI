




#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <mpi.h>

#include "Resolution.h"

extern int Ncube;    
extern int NBC_plus;

void BCM_Ghostcell_plus_Tem
(
// =============================================================================== //
int myid,

int *NNBC,

double Th,

double (*Dweight) = new double[NBC_plus*8+1],
double (*Nweight) = new double[NBC_plus*8+1],
int (*GCindex) = new int[NBC_plus*4+1],
int (*IPsur) = new int[NBC_plus*4+1],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// =============================================================================== //
)

{
	
#include "BCM.h"
#include "prm.h"
	
int gicube,gi,gj,gk;
int iNBC, Ntemp, ii, jj, kk, iw,  icase;

double rho, U, V, W, VV, P, T, temp;

double u1[8], u2[8], u3[8], u4[8], u5[8];

double Dw[8],Nw[8];



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


		for( iw = 0; iw < 8; iw ++ ) Dw[iw] = Dweight[Ntemp+1+iw];
		for( iw = 0; iw < 8; iw ++ ) Nw[iw] = Nweight[Ntemp+1+iw];

		Ntemp = (iNBC-1)*4;

		gicube = GCindex[Ntemp+1];
		gi = GCindex[Ntemp+2];
		gj = GCindex[Ntemp+3];
		gk = GCindex[Ntemp+4];



		icase = (gi-i)*4+(gj-j)*2+gk-k+1;

		
		for (ii = i; ii <= i+1; ii++) {
			for (jj = j; jj <= j+1; jj++) {
				for (kk = k; kk <= k+1; kk++) {  

					if( icase == (ii-i)*4 + 2*(jj-j) + (kk-k) + 1 ) {

						T = Th;
						U = 0.0;
						V = 0.0;
						W = 0.0;
						P = 0.0;

					}
					else {

						rho = U1_[icube][ii][jj][kk][0];
						U = U1_[icube][ii][jj][kk][1]/rho;
						V = U1_[icube][ii][jj][kk][2]/rho;
						W = U1_[icube][ii][jj][kk][3]/rho;
						P = ( U1_[icube][ii][jj][kk][4]-0.5*rho*(U*U+V*V+W*W) )*(K-1);
						T = P/rho/R;

					}

						u1[(ii-i)*4 + 2*(jj-j) + (kk-k)] = T;
						u2[(ii-i)*4 + 2*(jj-j) + (kk-k)] = U;
						u3[(ii-i)*4 + 2*(jj-j) + (kk-k)] = V;
						u4[(ii-i)*4 + 2*(jj-j) + (kk-k)] = W;
						u5[(ii-i)*4 + 2*(jj-j) + (kk-k)] = P;
					
					
				}
			}
		}


		T = 0.0;
		U = V = W = 0.0;
		P = 0.0;

		for( ii = 0; ii < 8; ii ++ ) {

			
			T = T + u1[ii]*Dw[ii];
			U = U + u2[ii]*Dw[ii];
			V = V + u3[ii]*Dw[ii];
			W = W + u4[ii]*Dw[ii];
			P = P + u5[ii]*Nw[ii];

		}

			


		VV = U*U+V*V+W*W;

		rho = P/T/R;

		U1_[gicube][gi][gj][gk][0] = rho;
		U1_[gicube][gi][gj][gk][1] = rho*U;
		U1_[gicube][gi][gj][gk][2] = rho*V;
		U1_[gicube][gi][gj][gk][3] = rho*W;
		U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

	}

	
}