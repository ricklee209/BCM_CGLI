




#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <mpi.h>

#include "Resolution.h"

extern int Nctbe;    
extern int NBC_plus;

void BCM_Ghostcell_plus
(
// =============================================================================== //
int myid,

int *NNBC,

double Th,

double (*weight) = new double[NBC_plus*8+1],
int (*GCindex) = new int[NBC_plus*4+1],
int (*IPsur) = new int[NBC_plus*4+1],

int (*FWS)[X_size][Y_size][Z_size] = new int[Nctbe][X_size][Y_size][Z_size],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Nctbe][X_size][Y_size][Z_size][Ndim]

// =============================================================================== //
)
	
{
	
#include "BCM.h"
#include "prm.h"

int gicube,gi,gj,gk;
int iNBC, Ntemp;

double rho, U, V, W, VV, P, T;

double r0,r1,r2,r3,r4,r5,r6,r7;
double u0,u1,u2,u3,u4,u5,u6,u7;
double v0,v1,v2,v3,v4,v5,v6,v7;
double w0,w1,w2,w3,w4,w5,w6,w7;
double vv0,vv1,vv2,vv3,vv4,vv5,vv6,vv7;
double p0,p1,p2,p3,p4,p5,p6,p7;
double t0,t1,t2,t3,t4,t5,t6,t7;

double wc1, wc2, wc3, wc4, wc5, wc6, wc7, wc8;

double er_p = 0.00001;

// ---- GCindex => IPsur => weight ---- //

#pragma omp parallel for private(Ntemp,icube,i,j,k,\
	r0, r1, r2, r3, r4, r5, r6, r7,\
	u0, u1, u2, u3, u4, u5, u6, u7,\
	v0, v1, v2, v3, v4, v5, v6, v7,\
	w0, w1, w2, w3, w4, w5, w6, w7,\
	p0, p1, p2, p3, p4, p5, p6, p7,\
	vv0, vv1, vv2, vv3, vv4, vv5, vv6, vv7,\
	t0, t1, t2, t3, t4, t5, t6, t7,\
	wc1,wc2,wc3,wc4,wc5,wc6,wc7,wc8,\
	gicube,gi,gj,gk,\
	rho,P,U,V,W,VV,T\
	)
	
	for (iNBC = 1; iNBC <= *NNBC; iNBC++) {

		Ntemp = (iNBC-1)*4;

		icube = IPsur[Ntemp+1];
		i = IPsur[Ntemp+2];
		j = IPsur[Ntemp+3];
		k = IPsur[Ntemp+4];
		
		r0 = U1_[icube][i][j][k][0];
		u0 = U1_[icube][i][j][k][1]/r0;
		v0 = U1_[icube][i][j][k][2]/r0;
		w0 = U1_[icube][i][j][k][3]/r0;
		vv0 = u0*u0+v0*v0+w0*w0;
		p0 = (U1_[icube][i][j][k][4]-0.5*r0*vv0)*(K-1);
		t0 = p0/r0/R;
	
		
		r1 = U1_[icube][i+1][j][k][0];
		u1 = U1_[icube][i+1][j][k][1]/r1;
		v1 = U1_[icube][i+1][j][k][2]/r1;
		w1 = U1_[icube][i+1][j][k][3]/r1;
		vv1 = u1*u1+v1*v1+w1*w1;
		p1 = (U1_[icube][i+1][j][k][4]-0.5*r1*vv1)*(K-1);
		t1 = p1/r1/R;

		r2 = U1_[icube][i][j+1][k][0];
		u2 = U1_[icube][i][j+1][k][1]/r2;
		v2 = U1_[icube][i][j+1][k][2]/r2;
		w2 = U1_[icube][i][j+1][k][3]/r2;
		vv2 = u2*u2+v2*v2+w2*w2;
		p2 = (U1_[icube][i][j+1][k][4]-0.5*r2*vv2)*(K-1);
		t2 = p2/r2/R;

		r3 = U1_[icube][i][j][k+1][0];
		u3 = U1_[icube][i][j][k+1][1]/r3;
		v3 = U1_[icube][i][j][k+1][2]/r3;
		w3 = U1_[icube][i][j][k+1][3]/r3;
		vv3 = u3*u3+v3*v3+w3*w3;
		p3 = (U1_[icube][i][j][k+1][4]-0.5*r3*vv3)*(K-1);
		t3 = p3/r3/R;

		r4 = U1_[icube][i+1][j+1][k][0];
		u4 = U1_[icube][i+1][j+1][k][1]/r4;
		v4 = U1_[icube][i+1][j+1][k][2]/r4;
		w4 = U1_[icube][i+1][j+1][k][3]/r4;
		vv4 = u4*u4+v4*v4+w4*w4;
		p4 = (U1_[icube][i+1][j+1][k][4]-0.5*r4*vv4)*(K-1);
		t4 = p4/r4/R;

		r5 = U1_[icube][i+1][j][k+1][0];
		u5 = U1_[icube][i+1][j][k+1][1]/r5;
		v5 = U1_[icube][i+1][j][k+1][2]/r5;
		w5 = U1_[icube][i+1][j][k+1][3]/r5;
		vv5 = u5*u5+v5*v5+w5*w5;
		p5 = (U1_[icube][i+1][j][k+1][4]-0.5*r5*vv5)*(K-1);
		t5 = p5/r5/R;

		r6 = U1_[icube][i][j+1][k+1][0];
		u6 = U1_[icube][i][j+1][k+1][1]/r6;
		v6 = U1_[icube][i][j+1][k+1][2]/r6;
		w6 = U1_[icube][i][j+1][k+1][3]/r6;
		vv6 = u6*u6+v6*v6+w6*w6;
		p6 = (U1_[icube][i][j+1][k+1][4]-0.5*r6*vv6)*(K-1);
		t6 = p6/r6/R;


		r7 = U1_[icube][i+1][j+1][k+1][0];
		u7 = U1_[icube][i+1][j+1][k+1][1]/r7;
		v7 = U1_[icube][i+1][j+1][k+1][2]/r7;
		w7 = U1_[icube][i+1][j+1][k+1][3]/r7;
		vv7 = u7*u7+v7*v7+w7*w7;
		p7 = (U1_[icube][i+1][j+1][k+1][4]-0.5*r7*vv7)*(K-1);
		t7 = p7/r7/R;


		Ntemp = (iNBC-1)*8;

		wc1 = weight[Ntemp+1];
		wc2 = weight[Ntemp+2];
		wc3 = weight[Ntemp+3];
		wc4 = weight[Ntemp+4];
		wc5 = weight[Ntemp+5];
		wc6 = weight[Ntemp+6];
		wc7 = weight[Ntemp+7];
		wc8 = weight[Ntemp+8];


		
		Ntemp = (iNBC-1)*4;

		gicube = GCindex[Ntemp+1];
		gi = GCindex[Ntemp+2];
		gj = GCindex[Ntemp+3];
		gk = GCindex[Ntemp+4];

		
		
		rho = wc1*(r0-rho0)+wc2*(r1-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0)+r0;

		P = wc1*(p0-p0)+wc2*(p1-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc7*(p6-p0)+wc8*(p7-p0)+p0;
		
		U = wc1*u0+wc2*u1+wc3*u2+wc4*u3+wc5*u4+wc6*u5+wc7*u6+wc8*u7;
		V = wc1*v0+wc2*v1+wc3*v2+wc4*v3+wc5*v4+wc6*v5+wc7*v6+wc8*v7;
		W = wc1*w0+wc2*w1+wc3*w2+wc4*w3+wc5*w4+wc6*w5+wc7*w6+wc8*w7;
		
		VV = U*U+V*V+W*W;

		U1_[gicube][gi][gj][gk][0] = rho;
		U1_[gicube][gi][gj][gk][1] = 0.5*rho*U;
		U1_[gicube][gi][gj][gk][2] = 0.5*rho*V;
		U1_[gicube][gi][gj][gk][3] = 0.5*rho*W;
		U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*0.25*VV;
		
		
		
		if ( gi==(i  ) & gj==(j  ) & gk==(k  ) ) { 
		
			rho = (wc2*(r1-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0))/(1-wc1+er_p)+r0;
			P = (wc2*(p1-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc7*(p6-p0)+wc8*(p7-p0))/(1-wc1+er_p)+p0;

			U = (wc2*u1+wc3*u2+wc4*u3+wc5*u4+wc6*u5+wc7*u6+wc8*u7)/(2-wc1);
			V = (wc2*v1+wc3*v2+wc4*v3+wc5*v4+wc6*v5+wc7*v6+wc8*v7)/(2-wc1);
			W = (wc2*w1+wc3*w2+wc4*w3+wc5*w4+wc6*w5+wc7*w6+wc8*w7)/(2-wc1);
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			}
			
		if ( gi==(i+1) & gj==(j  ) & gk==(k  ) ) { 
		
			rho = (wc1*(r0-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0))/(1-wc2+er_p)+r0;
			P = (wc1*(p0-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc7*(p6-p0)+wc8*(p7-p0))/(1-wc2+er_p)+p0;

			U = (wc1*u0+wc3*u2+wc4*u3+wc5*u4+wc6*u5+wc7*u6+wc8*u7)/(2-wc2);
			V = (wc1*v0+wc3*v2+wc4*v3+wc5*v4+wc6*v5+wc7*v6+wc8*v7)/(2-wc2);
			W = (wc1*w0+wc3*w2+wc4*w3+wc5*w4+wc6*w5+wc7*w6+wc8*w7)/(2-wc2);
			
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			}
			
		if ( gi==(i  ) & gj==(j+1) & gk==(k  ) ) { 
		
			rho = (wc1*(r0-rho0)+wc2*(r1-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0))/(1-wc3+er_p)+r0;
			P = (wc1*(p0-p0)+wc2*(p1-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc7*(p6-p0)+wc8*(p7-p0))/(1-wc3+er_p)+p0;

			U = (wc1*u0+wc2*u1+wc4*u3+wc5*u4+wc6*u5+wc7*u6+wc8*u7)/(2-wc3);
			V = (wc1*v0+wc2*v1+wc4*v3+wc5*v4+wc6*v5+wc7*v6+wc8*v7)/(2-wc3);
			W = (wc1*w0+wc2*w1+wc4*w3+wc5*w4+wc6*w5+wc7*w6+wc8*w7)/(2-wc3);
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			}
			
		if ( gi==(i  ) & gj==(j  ) & gk==(k+1) ) { 
		
			rho = (wc1*(r0-rho0)+wc2*(r1-rho0)+wc3*(r2-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0))/(1-wc4+er_p)+r0;
			P = (wc1*(p0-p0)+wc2*(p1-p0)+wc3*(p2-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc7*(p6-p0)+wc8*(p7-p0))/(1-wc4+er_p)+p0;

			U = (wc1*u0+wc2*u1+wc3*u2+wc5*u4+wc6*u5+wc7*u6+wc8*u7)/(2-wc4);
			V = (wc1*v0+wc2*v1+wc3*v2+wc5*v4+wc6*v5+wc7*v6+wc8*v7)/(2-wc4);
			W = (wc1*w0+wc2*w1+wc3*w2+wc5*w4+wc6*w5+wc7*w6+wc8*w7)/(2-wc4);
			
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			
			}
			
		if ( gi==(i+1) & gj==(j+1) & gk==(k  ) ) { 
			
			rho = (wc1*(r0-rho0)+wc2*(r1-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0))/(1-wc5+er_p)+r0;
			P = (wc1*(p0-p0)+wc2*(p1-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc6*(p5-p0)+wc7*(p6-p0)+wc8*(p7-p0))/(1-wc5+er_p)+p0;

			U = (wc1*u0+wc2*u1+wc3*u2+wc4*u3+wc6*u5+wc7*u6+wc8*u7)/(2-wc5);
			V = (wc1*v0+wc2*v1+wc3*v2+wc4*v3+wc6*v5+wc7*v6+wc8*v7)/(2-wc5);
			W = (wc1*w0+wc2*w1+wc3*w2+wc4*w3+wc6*w5+wc7*w6+wc8*w7)/(2-wc5);
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			}
			
		if ( gi==(i+1) & gj==(j  ) & gk==(k+1) ) { 
			
			rho = (wc1*(r0-rho0)+wc2*(r1-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc7*(r6-rho0)+wc8*(r7-rho0))/(1-wc6+er_p)+r0;
			P = (wc1*(p0-p0)+wc2*(p1-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc7*(p6-p0)+wc8*(p7-p0))/(1-wc6+er_p)+p0;

			U = (wc1*u0+wc2*u1+wc3*u2+wc4*u3+wc5*u4+wc7*u6+wc8*u7)/(2-wc6);
			V = (wc1*v0+wc2*v1+wc3*v2+wc4*v3+wc5*v4+wc7*v6+wc8*v7)/(2-wc6);
			W = (wc1*w0+wc2*w1+wc3*w2+wc4*w3+wc5*w4+wc7*w6+wc8*w7)/(2-wc6);
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			}
			
		if ( gi==(i  ) & gj==(j+1) & gk==(k+1) ) { 

			rho = (wc1*(r0-rho0)+wc2*(r1-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc8*(r7-rho0))/(1-wc7+er_p)+r0;
			P = (wc1*(p0-p0)+wc2*(p1-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc8*(p7-p0))/(1-wc7+er_p)+p0;

			U = (wc1*u0+wc2*u1+wc3*u2+wc4*u3+wc5*u4+wc6*u5+wc8*u7)/(2-wc7);
			V = (wc1*v0+wc2*v1+wc3*v2+wc4*v3+wc5*v4+wc6*v5+wc8*v7)/(2-wc7);
			W = (wc1*w0+wc2*w1+wc3*w2+wc4*w3+wc5*w4+wc6*w5+wc8*w7)/(2-wc7);
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i+1) & gj==(j+1) & gk==(k+1) ) {
			
			
			rho = (wc1*(r0-rho0)+wc2*(r1-rho0)+wc3*(r2-rho0)+wc4*(r3-rho0)+wc5*(r4-rho0)+wc6*(r5-rho0)+wc7*(r6-rho0))/(1-wc8+er_p)+r0;
			P = (wc1*(p0-p0)+wc2*(p1-p0)+wc3*(p2-p0)+wc4*(p3-p0)+wc5*(p4-p0)+wc6*(p5-p0)+wc7*(p6-p0))/(1-wc8+er_p)+p0;

			U = (wc1*u0+wc2*u1+wc3*u2+wc4*u3+wc5*u4+wc6*u5+wc7*u6)/(2-wc8);
			V = (wc1*v0+wc2*v1+wc3*v2+wc4*v3+wc5*v4+wc6*v5+wc7*v6)/(2-wc8);
			W = (wc1*w0+wc2*w1+wc3*w2+wc4*w3+wc5*w4+wc6*w5+wc7*w6)/(2-wc8);
			
			VV = U*U+V*V+W*W;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;


			}

		
	}




}