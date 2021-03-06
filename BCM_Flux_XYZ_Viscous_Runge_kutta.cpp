



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;   

void BCM_Flux_XYZ_Viscous_Runge_kutta 
	(
	// ================================================================================ //
	int myid,

	int ncube,

	int RK,

	double deltaT,

	double deltaTau,

	double e,


	int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

	int (*csl) = new int[Ncube],

	double (*cube_size) = new double[Ncube],

	double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*U1p1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*U1p2)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*Rku1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*Residual1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

	double (*er) = new double[8]

// ================================================================================ //
)

{



#include "BCM.h"
#include "prm.h"


	int IF,
		IFi1,IFi2,IFi3,IFi4,IFi5,
		IF_i1,IF_i2,IF_i3,IF_i4,IF_i5,

		IFj1,IFj2,IFj3,IFj4,IFj5,
		IF_j1,IF_j2,IF_j3,IF_j4,IF_j5,

		IFk1,IFk2,IFk3,IFk4,IFk5,
		IF_k1,IF_k2,IF_k3,IF_k4,IF_k5;

	int id1,id2,id3,id4,di5;


	double rho,U,V,W,VV,P,C,T,h,H;

	double 	irho,iU,iV,iW,iP,iT,
		jrho,jU,jV,jW,jP,jT,
		krho,kU,kV,kW,kP,kT,

		rhoi,Ui,Vi,Wi,Pi,Ti,
		rhoj,Uj,Vj,Wj,Pj,Tj,
		rhok,Uk,Vk,Wk,Pk,Tk,

		ijrho,ijU,ijV,ijW,ijP,ijT,
		jkrho,jkU,jkV,jkW,jkP,jkT,
		ikrho,ikU,ikV,ikW,ikP,ikT,

		irhoj,iUj,iVj,iWj,iPj,iTj,
		irhok,iUk,iVk,iWk,iPk,iTk,

		jrhoi,jUi,jVi,jWi,jPi,jTi,
		jrhok,jUk,jVk,jWk,jPk,jTk,

		krhoj,kUj,kVj,kWj,kPj,kTj,
		krhoi,kUi,kVi,kWi,kPi,kTi,

		rhoij,Uij,Vij,Wij,Pij,Tij,
		rhoik,Uik,Vik,Wik,Pik,Tik,
		rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,

		du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,
		dT_dx,dT_dy,dT_dz,

		mu_E,mu_T,Pr_E;

	double sigma_xx,sigma_xy,sigma_xz,
		sigma_yx,sigma_yy,sigma_yz,
		sigma_zx,sigma_zy,sigma_zz;


	double invXI,invET,invZT;

	double Tx,Ty,Tz,
		Ux,Uy,Uz,
		Vx,Vy,Vz,
		Wx,Wy,Wz;

	double LL1,LL2,LL3,LL4,LL5,
		LL1i,LL2i,LL3i,LL4i,LL5i,
		ML1,ML2,ML3,ML4,ML5,
		ML1j,ML2j,ML3j,ML4j,ML5j,
		NL1,NL2,NL3,NL4,NL5,
		NL1k,NL2k,NL3k,NL4k,NL5k,	
		MR1,MR2,MR3,MR4,MR5,

		ML1i,ML2i,ML3i,ML4i,ML5i,
		MR1i,MR2i,MR3i,MR4i,MR5i;





	double beta,S,deltaU,deltaP;

	double temp,temp1,temp2,temp3,temp4,temp5,temp6;


	double u1,u2,u3,u4,u5,

		u1i1,u2i1,u3i1,u4i1,u5i1,
		u1i2,u2i2,u3i2,u4i2,u5i2,
		u1i3,u2i3,u3i3,u4i3,u5i3,

		u1_i1,u2_i1,u3_i1,u4_i1,u5_i1,
		u1_i2,u2_i2,u3_i2,u4_i2,u5_i2,
		u1_i3,u2_i3,u3_i3,u4_i3,u5_i3,

		u1j1,u2j1,u3j1,u4j1,u5j1,
		u1j2,u2j2,u3j2,u4j2,u5j2,
		u1j3,u2j3,u3j3,u4j3,u5j3,

		u1_j1,u2_j1,u3_j1,u4_j1,u5_j1,
		u1_j2,u2_j2,u3_j2,u4_j2,u5_j2,
		u1_j3,u2_j3,u3_j3,u4_j3,u5_j3,

		u1k1,u2k1,u3k1,u4k1,u5k1,
		u1k2,u2k2,u3k2,u4k2,u5k2,
		u1k3,u2k3,u3k3,u4k3,u5k3,

		u1_k1,u2_k1,u3_k1,u4_k1,u5_k1,
		u1_k2,u2_k2,u3_k2,u4_k2,u5_k2,
		u1_k3,u2_k3,u3_k3,u4_k3,u5_k3;

	double uu1,uu2,uu3,uu4,uu5;

	double u1q,u2q,u3q,u4q,u5q;


	double _rho,_U,_V,_W,_VV,_P,_T,_C,_H;
	double rho_,U_,V_,W_,VV_,P_,T_,C_,H_;
	double dU1,dU2,dU3,dU4,dU5;
	double Fav1,Fav2,Fav3,Fav4,Fav5;

	double inFx1,inFx2,inFx3,inFx4,inFx5,
		inFx1i,inFx2i,inFx3i,inFx4i,inFx5i,
		inFy1,inFy2,inFy3,inFy4,inFy5,
		inFy1i,inFy2i,inFy3i,inFy4i,inFy5i,
		inFz1,inFz2,inFz3,inFz4,inFz5,
		inFz1i,inFz2i,inFz3i,inFz4i,inFz5i;

	double Rp1,Rp2,Rp3,Rp4,Rp5,
		Rf1,Rf2,Rf3,Rf4,Rf5,
		Rfx1,Rfx2,Rfx3,Rfx4,Rfx5,
		Rfy1,Rfy2,Rfy3,Rfy4,Rfy5,
		Rfz1,Rfz2,Rfz3,Rfz4,Rfz5,
		Rk1,Rk2,Rk3,Rk4,Rk5,
		vF1,vF2,vF3,vF4,vF5,
		Sx,Sy,Sz;

	double d11,d12,d13,d14,d15,
		d21,d22,d23,d24,d25,
		d31,d32,d33,d34,d35,
		d41,d42,d43,d44,d45,
		d51,d52,d53,d54,d55;

	double rhoold,Uold,Vold,Wold,VVold,Pold,Told;
	

	double dPx,dPy,dPz,dPxi,dPyi,dPzi;

	double e1 = 0;
	double e2 = 0;
	double e3 = 0;
	double e4 = 0;
	double e5 = 0;
	
	double e6 = 0;
	double e7 = 0;

	double theda_p, U_p, C_p;

	double A2 = -5.0/9.0;
	double A3 = -153.0/128.0;
	double B1 = 1.0/3.0;
	double B2 = 15.0/16.0;
	double B3 = 8.0/15.0;
	


#pragma omp parallel for private(\
	IF,\
	IFi1,IFi2,IFi3,IFi4,IFi5,\
	IF_i1,IF_i2,IF_i3,IF_i4,IF_i5,\
	IFj1,IFj2,IFj3,IFj4,IFj5,\
	IF_j1,IF_j2,IF_j3,IF_j4,IF_j5,\
	IFk1,IFk2,IFk3,IFk4,IFk5,\
	IF_k1,IF_k2,IF_k3,IF_k4,IF_k5,\
	id1,id2,id3,id4,di5,\
	rho,U,V,W,VV,P,C,T,h,H,\
	irho,iU,iV,iW,iP,iT,\
	jrho,jU,jV,jW,jP,jT,\
	krho,kU,kV,kW,kP,kT,\
	rhoi,Ui,Vi,Wi,Pi,Ti,\
	rhoj,Uj,Vj,Wj,Pj,Tj,\
	rhok,Uk,Vk,Wk,Pk,Tk,\
	ijrho,ijU,ijV,ijW,ijP,ijT,\
	jkrho,jkU,jkV,jkW,jkP,jkT,\
	ikrho,ikU,ikV,ikW,ikP,ikT,\
	irhoj,iUj,iVj,iWj,iPj,iTj,\
	irhok,iUk,iVk,iWk,iPk,iTk,\
	jrhoi,jUi,jVi,jWi,jPi,jTi,\
	jrhok,jUk,jVk,jWk,jPk,jTk,\
	krhoj,kUj,kVj,kWj,kPj,kTj,\
	krhoi,kUi,kVi,kWi,kPi,kTi,\
	rhoij,Uij,Vij,Wij,Pij,Tij,\
	rhoik,Uik,Vik,Wik,Pik,Tik,\
	rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,\
	du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,\
	dT_dx,dT_dy,dT_dz,\
	mu_E,mu_T,Pr_E,\
	sigma_xx,sigma_xy,sigma_xz,\
	sigma_yx,sigma_yy,sigma_yz,\
	sigma_zx,sigma_zy,sigma_zz,\
	invXI,invET,invZT,\
	Tx,Ty,Tz,\
	Ux,Uy,Uz,\
	Vx,Vy,Vz,\
	Wx,Wy,Wz,\
	LL1,LL2,LL3,LL4,LL5,\
	LL1i,LL2i,LL3i,LL4i,LL5i,\
	ML1,ML2,ML3,ML4,ML5,\
	ML1j,ML2j,ML3j,ML4j,ML5j,\
	NL1,NL2,NL3,NL4,NL5,\
	NL1k,NL2k,NL3k,NL4k,NL5k,\
	MR1,MR2,MR3,MR4,MR5,\
	ML1i,ML2i,ML3i,ML4i,ML5i,\
	MR1i,MR2i,MR3i,MR4i,MR5i,\
	beta,S,deltaU,deltaP,\
	temp,temp1,temp2,temp3,temp4,temp5,temp6,\
	u1,u2,u3,u4,u5,\
	u1i1,u2i1,u3i1,u4i1,u5i1,\
	u1i2,u2i2,u3i2,u4i2,u5i2,\
	u1i3,u2i3,u3i3,u4i3,u5i3,\
	u1_i1,u2_i1,u3_i1,u4_i1,u5_i1,\
	u1_i2,u2_i2,u3_i2,u4_i2,u5_i2,\
	u1_i3,u2_i3,u3_i3,u4_i3,u5_i3,\
	u1j1,u2j1,u3j1,u4j1,u5j1,\
	u1j2,u2j2,u3j2,u4j2,u5j2,\
	u1j3,u2j3,u3j3,u4j3,u5j3,\
	u1_j1,u2_j1,u3_j1,u4_j1,u5_j1,\
	u1_j2,u2_j2,u3_j2,u4_j2,u5_j2,\
	u1_j3,u2_j3,u3_j3,u4_j3,u5_j3,\
	u1k1,u2k1,u3k1,u4k1,u5k1,\
	u1k2,u2k2,u3k2,u4k2,u5k2,\
	u1k3,u2k3,u3k3,u4k3,u5k3,\
	u1_k1,u2_k1,u3_k1,u4_k1,u5_k1,\
	u1_k2,u2_k2,u3_k2,u4_k2,u5_k2,\
	u1_k3,u2_k3,u3_k3,u4_k3,u5_k3,\
	uu1,uu2,uu3,uu4,uu5,\
	u1q,u2q,u3q,u4q,u5q,\
	_rho,_U,_V,_W,_VV,_P,_T,_C,_H,\
	rho_,U_,V_,W_,VV_,P_,T_,C_,H_,\
	dU1,dU2,dU3,dU4,dU5,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	inFx1,inFx2,inFx3,inFx4,inFx5,\
	inFx1i,inFx2i,inFx3i,inFx4i,inFx5i,\
	inFy1,inFy2,inFy3,inFy4,inFy5,\
	inFy1i,inFy2i,inFy3i,inFy4i,inFy5i,\
	inFz1,inFz2,inFz3,inFz4,inFz5,\
	inFz1i,inFz2i,inFz3i,inFz4i,inFz5i,\
	Rp1,Rp2,Rp3,Rp4,Rp5,\
	Rf1,Rf2,Rf3,Rf4,Rf5,\
	Rfx1,Rfx2,Rfx3,Rfx4,Rfx5,\
	Rfy1,Rfy2,Rfy3,Rfy4,Rfy5,\
	Rfz1,Rfz2,Rfz3,Rfz4,Rfz5,\
	Rk1,Rk2,Rk3,Rk4,Rk5,\
	vF1,vF2,vF3,vF4,vF5,\
	Sx,Sy,Sz,\
	d11,d12,d13,d14,d15,\
	d21,d22,d23,d24,d25,\
	d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,\
	d51,d52,d53,d54,d55,\
	i,j,k,\
	dx,dy,dz,\
	dPx,dPy,dPz,dPxi,dPyi,dPzi\
	) 


	for (icube = 1; icube < ncube; icube++) {  

		dx = dy = dz = cube_size[icube]/NcubeX;
		invXI = invET = invZT = 1./dx;

		for (i = 2; i <= nx; i++) {
			for (j = 2; j <= ny; j++) {
				for (k = 2; k <= nz; k++) {


					uu1 = U1[icube][i][j][k][0];
					uu2 = U1[icube][i][j][k][1];
					uu3 = U1[icube][i][j][k][2];
					uu4 = U1[icube][i][j][k][3];
					uu5 = U1[icube][i][j][k][4];

					u1q = U1q[icube][i][j][k][0];
					u2q = U1q[icube][i][j][k][1];
					u3q = U1q[icube][i][j][k][2];
					u4q = U1q[icube][i][j][k][3];
					u5q = U1q[icube][i][j][k][4];

					IF = FWS[icube][i][j][k];
					u1 = U1_[icube][i][j][k][0];
					u2 = U1_[icube][i][j][k][1];
					u3 = U1_[icube][i][j][k][2];
					u4 = U1_[icube][i][j][k][3];
					u5 = U1_[icube][i][j][k][4];

					Rp1 = -(3*u1-4*uu1+u1q)/(2*deltaT);
					Rp2 = -(3*u2-4*uu2+u2q)/(2*deltaT);
					Rp3 = -(3*u3-4*uu3+u3q)/(2*deltaT);
					Rp4 = -(3*u4-4*uu4+u4q)/(2*deltaT);
					Rp5 = -(3*u5-4*uu5+u5q)/(2*deltaT);




					// -----------------------------------------------------------------//
					// -------------------------- Z-direction --------------------------//

					IF_k1 = FWS[icube][i][j][k-1];
					IF_k2 = FWS[icube][i][j][k-2];
					IFk1 = FWS[icube][i][j][k+1];
					IFk2 = FWS[icube][i][j][k+2];

					if(k > 2) {

						IF_k3 = FWS[icube][i][j][k-3];
						u1_k3 = U1_[icube][i][j][k-3][0];
						u2_k3 = U1_[icube][i][j][k-3][1];
						u3_k3 = U1_[icube][i][j][k-3][2];
						u4_k3 = U1_[icube][i][j][k-3][3];
						u5_k3 = U1_[icube][i][j][k-3][4];
					}

					u1_k2 = U1_[icube][i][j][k-2][0];
					u2_k2 = U1_[icube][i][j][k-2][1];
					u3_k2 = U1_[icube][i][j][k-2][2];
					u4_k2 = U1_[icube][i][j][k-2][3];
					u5_k2 = U1_[icube][i][j][k-2][4];

					
					u1_k1 = U1_[icube][i][j][k-1][0];
					u2_k1 = U1_[icube][i][j][k-1][1];
					u3_k1 = U1_[icube][i][j][k-1][2];
					u4_k1 = U1_[icube][i][j][k-1][3];
					u5_k1 = U1_[icube][i][j][k-1][4];

					u1k1 = U1_[icube][i][j][k+1][0];
					u2k1 = U1_[icube][i][j][k+1][1];
					u3k1 = U1_[icube][i][j][k+1][2];
					u4k1 = U1_[icube][i][j][k+1][3];
					u5k1 = U1_[icube][i][j][k+1][4];

					u1k2 = U1_[icube][i][j][k+2][0]; 
					u2k2 = U1_[icube][i][j][k+2][1]; 
					u3k2 = U1_[icube][i][j][k+2][2]; 
					u4k2 = U1_[icube][i][j][k+2][3]; 
					u5k2 = U1_[icube][i][j][k+2][4]; 


					if(k < nz) { 

						IFk3 = FWS[icube][i][j][k+3];

						u1k3 = U1_[icube][i][j][k+3][0];
						u2k3 = U1_[icube][i][j][k+3][1];
						u3k3 = U1_[icube][i][j][k+3][2];
						u4k3 = U1_[icube][i][j][k+3][3];
						u5k3 = U1_[icube][i][j][k+3][4];

						
					}

					// -------------------------- Z-direction --------------------------//
					// -----------------------------------------------------------------//




					// ------------------------------------------------------------------- //
					// ----------------------------- MUSCL-Z ----------------------------- //




































































					if( IFk1 == IFLUID  && IF_k1 == IFLUID) { 

						if ( k==2 | IF_k2 != IFLUID ) {

							ML1 = 1./6*(-1*u1_k2+5*u1_k1+2*u1);
							ML2 = 1./6*(-1*u2_k2+5*u2_k1+2*u2);
							ML3 = 1./6*(-1*u3_k2+5*u3_k1+2*u3);
							ML4 = 1./6*(-1*u4_k2+5*u4_k1+2*u4);
							ML5 = 1./6*(-1*u5_k2+5*u5_k1+2*u5);

							MR1 = 1./6*(2*u1_k1+5*u1-1*u1k1);
							MR2 = 1./6*(2*u2_k1+5*u2-1*u2k1);
							MR3 = 1./6*(2*u3_k1+5*u3-1*u3k1);
							MR4 = 1./6*(2*u4_k1+5*u4-1*u4k1);
							MR5 = 1./6*(2*u5_k1+5*u5-1*u5k1);

						}
						else {

							ML1 = 1./60*(2*u1_k3-13*u1_k2+47*u1_k1+27*u1-3*u1k1);
							ML2 = 1./60*(2*u2_k3-13*u2_k2+47*u2_k1+27*u2-3*u2k1);
							ML3 = 1./60*(2*u3_k3-13*u3_k2+47*u3_k1+27*u3-3*u3k1);
							ML4 = 1./60*(2*u4_k3-13*u4_k2+47*u4_k1+27*u4-3*u4k1);
							ML5 = 1./60*(2*u5_k3-13*u5_k2+47*u5_k1+27*u5-3*u5k1);

							MR1 = 1./60*(-3*u1_k2+27*u1_k1+47*u1-13*u1k1+2*u1k2);
							MR2 = 1./60*(-3*u2_k2+27*u2_k1+47*u2-13*u2k1+2*u2k2);
							MR3 = 1./60*(-3*u3_k2+27*u3_k1+47*u3-13*u3k1+2*u3k2);
							MR4 = 1./60*(-3*u4_k2+27*u4_k1+47*u4-13*u4k1+2*u4k2);
							MR5 = 1./60*(-3*u5_k2+27*u5_k1+47*u5-13*u5k1+2*u5k2);

						}
						
						
						
						
						if ( k==nx | IFk2 != IFLUID) {

							ML1i = 1./6*(-1*u1_k1+5*u1+2*u1k1);
							ML2i = 1./6*(-1*u2_k1+5*u2+2*u2k1);
							ML3i = 1./6*(-1*u3_k1+5*u3+2*u3k1);
							ML4i = 1./6*(-1*u4_k1+5*u4+2*u4k1);
							ML5i = 1./6*(-1*u5_k1+5*u5+2*u5k1);

							MR1i = 1./6*(2*u1+5*u1k1-1*u1k2);
							MR2i = 1./6*(2*u2+5*u2k1-1*u2k2);
							MR3i = 1./6*(2*u3+5*u3k1-1*u3k2);
							MR4i = 1./6*(2*u4+5*u4k1-1*u4k2);
							MR5i = 1./6*(2*u5+5*u5k1-1*u5k2);

						}
						else {

							ML1i = 1./60*(2*u1_k2-13*u1_k1+47*u1+27*u1k1-3*u1k2);
							ML2i = 1./60*(2*u2_k2-13*u2_k1+47*u2+27*u2k1-3*u2k2);
							ML3i = 1./60*(2*u3_k2-13*u3_k1+47*u3+27*u3k1-3*u3k2);
							ML4i = 1./60*(2*u4_k2-13*u4_k1+47*u4+27*u4k1-3*u4k2);
							ML5i = 1./60*(2*u5_k2-13*u5_k1+47*u5+27*u5k1-3*u5k2);

							MR1i = 1./60*(-3*u1_k1+27*u1+47*u1k1-13*u1k2+2*u1k3);
							MR2i = 1./60*(-3*u2_k1+27*u2+47*u2k1-13*u2k2+2*u2k3);
							MR3i = 1./60*(-3*u3_k1+27*u3+47*u3k1-13*u3k2+2*u3k3);
							MR4i = 1./60*(-3*u4_k1+27*u4+47*u4k1-13*u4k2+2*u4k3);
							MR5i = 1./60*(-3*u5_k1+27*u5+47*u5k1-13*u5k2+2*u5k3);

						}
						
						
						
						
						
					}
					else {

						ML1 = u1_k1;
						ML2 = u2_k1;
						ML3 = u3_k1;
						ML4 = u4_k1;
						ML5 = u5_k1;
						
						MR1 = u1;
						MR2 = u2;
						MR3 = u3;
						MR4 = u4;
						MR5 = u5;
						
						
						ML1i = u1;
						ML2i = u2;
						ML3i = u3;
						ML4i = u4;
						ML5i = u5;
						
						MR1i = u1k1;
						MR2i = u2k1;
						MR3i = u3k1;
						MR4i = u4k1;
						MR5i = u5k1;

					}

					// ----------------------------- MUSCL-Z ----------------------------- //
					// ------------------------------------------------------------------- //



					// ------------------------------------------------------------------ //
					// ----------------------------- Flux-Z ----------------------------- //

					/* -------------------------------- */
					/* ----------- Backward ----------- */


					/* jump dU */
					dU1 = MR1-ML1;
					dU2 = MR2-ML2;
					dU3 = MR3-ML3;
					dU4 = MR4-ML4;
					dU5 = MR5-ML5;

					/* lefr parameter */
					_rho = ML1;
					_U = ML2/_rho;
					_V = ML3/_rho;
					_W = ML4/_rho;
					_VV = _U*_U+_V*_V+_W*_W;
					_P = (ML5-0.5*_rho*_VV)*(K-1);
					_T = _P/_rho;
					_C = K*_P/_rho;
					_H = 0.5*_VV+_C/(K-1);

					/* right parameter */
					rho_ = MR1;
					U_ = MR2/rho_;
					V_ = MR3/rho_;
					W_ = MR4/rho_;
					VV_ = U_*U_+V_*V_+W_*W_;
					P_ = (MR5-0.5*rho_*VV_)*(K-1);
					T_ = P_/rho_;
					C_ = K*P_/rho_;
					H_ = 0.5*VV_+C_/(K-1);

					/* flux varibale */

					
					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					
					// temp5 = 1.0;
					// temp6 = 1.0;
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;
					VV = U*U+V*V+W*W;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
					P = rho*C/K;


					beta = max(VV/C,e);    // ---- theda ---- //
					
					temp = 0.5*(1+beta)*W;    // ---- U' ---- //

					S = 0.5*sqrt(4*beta*C+W*W*(1-beta)*(1-beta));   // ---- C' ---- //





					temp1 = (P_-_P)/rho/beta/C;

					temp2 = temp/S*(W_-_W);

					temp3 = 0.5*(1-beta)*W*temp/S;



					deltaU = (S-temp3-beta*fabs(W))*temp1+temp2;

					deltaP = temp/S*(P_-_P)+(S-fabs(W)+temp3)*rho*(W_-_W);

					/* artificial viscosity */
					Fav1 = fabs(W)*dU1+deltaU*rho;
					Fav2 = fabs(W)*dU2+deltaU*rho*U;
					Fav3 = fabs(W)*dU3+deltaU*rho*V;
					Fav4 = fabs(W)*dU4+deltaU*rho*W+deltaP;
					Fav5 = fabs(W)*dU5+deltaU*rho*H+deltaP*W;

					/* inviscid fluxes */

					inFz1 = 0.5*((_rho*_W+rho_*W_-Ep*Fav1));
					inFz2 = 0.5*((_rho*_U*_W+rho_*U_*W_)-Ep*Fav2);
					inFz3 = 0.5*((_rho*_V*_W+rho_*V_*W_)-Ep*Fav3);
					inFz4 = 0.5*((_rho*_W*_W+_P+rho_*W_*W_+P_)-Ep*Fav4);
					inFz5 = 0.5*((_W*(3.5*_P+0.5*_rho*_VV)+W_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

					dPz = 0.5*(_P+P_);
					/* ----------- Backward ----------- */
					/* -------------------------------- */


					/* --------------------------------- */
					/* ------------ Forward ------------ */


					/* lefr parameter */
					

					/* jump dU */
					dU1 = MR1i-ML1i;
					dU2 = MR2i-ML2i;
					dU3 = MR3i-ML3i;
					dU4 = MR4i-ML4i;
					dU5 = MR5i-ML5i;

					/* lefr parameter */
					_rho = ML1i;
					_U = ML2i/_rho;
					_V = ML3i/_rho;
					_W = ML4i/_rho;
					_VV = _U*_U+_V*_V+_W*_W;
					_P = (ML5i-0.5*_rho*_VV)*(K-1);
					_T = _P/_rho;
					_C = K*_P/_rho;
					_H = 0.5*_VV+_C/(K-1);

					/* right parameter */
					rho_ = MR1i;
					U_ = MR2i/rho_;
					V_ = MR3i/rho_;
					W_ = MR4i/rho_;
					VV_ = U_*U_+V_*V_+W_*W_;
					P_ = (MR5i-0.5*rho_*VV_)*(K-1);
					T_ = P_/rho_;
					C_ = K*P_/rho_;
					H_ = 0.5*VV_+C_/(K-1);

					/* flux varibale */

					
					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					
					// temp5 = 1.0;
					// temp6 = 1.0;
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;
					VV = U*U+V*V+W*W;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
					P = rho*C/K;


					beta = max(VV/C,e);    // ---- theda ---- //
					
					temp = 0.5*(1+beta)*W;    // ---- U' ---- //

					S = 0.5*sqrt(4*beta*C+W*W*(1-beta)*(1-beta));   // ---- C' ---- //





					temp1 = (P_-_P)/rho/beta/C;

					temp2 = temp/S*(W_-_W);

					temp3 = 0.5*(1-beta)*W*temp/S;



					deltaU = (S-temp3-beta*fabs(W))*temp1+temp2;

					deltaP = temp/S*(P_-_P)+(S-fabs(W)+temp3)*rho*(W_-_W);

					/* artificial viscosity */
					Fav1 = fabs(W)*dU1+deltaU*rho;
					Fav2 = fabs(W)*dU2+deltaU*rho*U;
					Fav3 = fabs(W)*dU3+deltaU*rho*V;
					Fav4 = fabs(W)*dU4+deltaU*rho*W+deltaP;
					Fav5 = fabs(W)*dU5+deltaU*rho*H+deltaP*W;

					/* inviscid fluxes */

					inFz1i = 0.5*((_rho*_W+rho_*W_-Ep*Fav1));
					inFz2i = 0.5*((_rho*_U*_W+rho_*U_*W_)-Ep*Fav2);
					inFz3i = 0.5*((_rho*_V*_W+rho_*V_*W_)-Ep*Fav3);
					inFz4i = 0.5*((_rho*_W*_W+_P+rho_*W_*W_+P_)-Ep*Fav4);
					inFz5i = 0.5*((_W*(3.5*_P+0.5*_rho*_VV)+W_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);


					dPzi = 0.5*(_P+P_);
					/* ------------ Forward ------------ */
					/* --------------------------------- */


					Rfz1 = inFz1i - inFz1;
					Rfz2 = inFz2i - inFz2;
					Rfz3 = inFz3i - inFz3;
					Rfz4 = inFz4i - inFz4;
					Rfz5 = inFz5i - inFz5;

					// ----------------------------- Flux-Z ----------------------------- //
					// ------------------------------------------------------------------ //






					// -----------------------------------------------------------------//
					// -------------------------- X-direction --------------------------//

					IF_i1 = FWS[icube][i-1][j][k];
					u1_i1 = U1_[icube][i-1][j][k][0];
					u2_i1 = U1_[icube][i-1][j][k][1];
					u3_i1 = U1_[icube][i-1][j][k][2];
					u4_i1 = U1_[icube][i-1][j][k][3];
					u5_i1 = U1_[icube][i-1][j][k][4];
					
					
					IF_i2 = FWS[icube][i-2][j][k];
					u1_i2 = U1_[icube][i-2][j][k][0];
					u2_i2 = U1_[icube][i-2][j][k][1];
					u3_i2 = U1_[icube][i-2][j][k][2];
					u4_i2 = U1_[icube][i-2][j][k][3];
					u5_i2 = U1_[icube][i-2][j][k][4];


					if(i > 2) {

						IF_i3 = FWS[icube][i-3][j][k];
						u1_i3 = U1_[icube][i-3][j][k][0];
						u2_i3 = U1_[icube][i-3][j][k][1];
						u3_i3 = U1_[icube][i-3][j][k][2];
						u4_i3 = U1_[icube][i-3][j][k][3];
						u5_i3 = U1_[icube][i-3][j][k][4];

					}

					IFi1 = FWS[icube][i+1][j][k];
					u1i1 = U1_[icube][i+1][j][k][0];
					u2i1 = U1_[icube][i+1][j][k][1];
					u3i1 = U1_[icube][i+1][j][k][2];
					u4i1 = U1_[icube][i+1][j][k][3];
					u5i1 = U1_[icube][i+1][j][k][4];

					IFi2 = FWS[icube][i+2][j][k];
					u1i2 = U1_[icube][i+2][j][k][0];
					u2i2 = U1_[icube][i+2][j][k][1];
					u3i2 = U1_[icube][i+2][j][k][2];
					u4i2 = U1_[icube][i+2][j][k][3];
					u5i2 = U1_[icube][i+2][j][k][4];


					if(i < nx)   {

						IFi3 = FWS[icube][i+3][j][k];
						u1i3 = U1_[icube][i+3][j][k][0];
						u2i3 = U1_[icube][i+3][j][k][1];
						u3i3 = U1_[icube][i+3][j][k][2];
						u4i3 = U1_[icube][i+3][j][k][3];
						u5i3 = U1_[icube][i+3][j][k][4];

					}


					// -------------------------- X-direction --------------------------//
					// -----------------------------------------------------------------//



					// ------------------------------------------------------------------- //
					// ----------------------------- MUSCL-X ----------------------------- //

































































					if( IFi1 == IFLUID  && IF_i1 == IFLUID) { 

						if ( i==2 | IF_i2 != IFLUID ) {

							ML1 = 1./6*(-1*u1_i2+5*u1_i1+2*u1);
							ML2 = 1./6*(-1*u2_i2+5*u2_i1+2*u2);
							ML3 = 1./6*(-1*u3_i2+5*u3_i1+2*u3);
							ML4 = 1./6*(-1*u4_i2+5*u4_i1+2*u4);
							ML5 = 1./6*(-1*u5_i2+5*u5_i1+2*u5);

							MR1 = 1./6*(2*u1_i1+5*u1-1*u1i1);
							MR2 = 1./6*(2*u2_i1+5*u2-1*u2i1);
							MR3 = 1./6*(2*u3_i1+5*u3-1*u3i1);
							MR4 = 1./6*(2*u4_i1+5*u4-1*u4i1);
							MR5 = 1./6*(2*u5_i1+5*u5-1*u5i1);

						}
						else {

							ML1 = 1./60*(2*u1_i3-13*u1_i2+47*u1_i1+27*u1-3*u1i1);
							ML2 = 1./60*(2*u2_i3-13*u2_i2+47*u2_i1+27*u2-3*u2i1);
							ML3 = 1./60*(2*u3_i3-13*u3_i2+47*u3_i1+27*u3-3*u3i1);
							ML4 = 1./60*(2*u4_i3-13*u4_i2+47*u4_i1+27*u4-3*u4i1);
							ML5 = 1./60*(2*u5_i3-13*u5_i2+47*u5_i1+27*u5-3*u5i1);

							MR1 = 1./60*(-3*u1_i2+27*u1_i1+47*u1-13*u1i1+2*u1i2);
							MR2 = 1./60*(-3*u2_i2+27*u2_i1+47*u2-13*u2i1+2*u2i2);
							MR3 = 1./60*(-3*u3_i2+27*u3_i1+47*u3-13*u3i1+2*u3i2);
							MR4 = 1./60*(-3*u4_i2+27*u4_i1+47*u4-13*u4i1+2*u4i2);
							MR5 = 1./60*(-3*u5_i2+27*u5_i1+47*u5-13*u5i1+2*u5i2);

						}
						
						
						
						
						if ( i==nx | IFi2 != IFLUID) {

							ML1i = 1./6*(-1*u1_i1+5*u1+2*u1i1);
							ML2i = 1./6*(-1*u2_i1+5*u2+2*u2i1);
							ML3i = 1./6*(-1*u3_i1+5*u3+2*u3i1);
							ML4i = 1./6*(-1*u4_i1+5*u4+2*u4i1);
							ML5i = 1./6*(-1*u5_i1+5*u5+2*u5i1);

							MR1i = 1./6*(2*u1+5*u1i1-1*u1i2);
							MR2i = 1./6*(2*u2+5*u2i1-1*u2i2);
							MR3i = 1./6*(2*u3+5*u3i1-1*u3i2);
							MR4i = 1./6*(2*u4+5*u4i1-1*u4i2);
							MR5i = 1./6*(2*u5+5*u5i1-1*u5i2);

						}
						else {

							ML1i = 1./60*(2*u1_i2-13*u1_i1+47*u1+27*u1i1-3*u1i2);
							ML2i = 1./60*(2*u2_i2-13*u2_i1+47*u2+27*u2i1-3*u2i2);
							ML3i = 1./60*(2*u3_i2-13*u3_i1+47*u3+27*u3i1-3*u3i2);
							ML4i = 1./60*(2*u4_i2-13*u4_i1+47*u4+27*u4i1-3*u4i2);
							ML5i = 1./60*(2*u5_i2-13*u5_i1+47*u5+27*u5i1-3*u5i2);

							MR1i = 1./60*(-3*u1_i1+27*u1+47*u1i1-13*u1i2+2*u1i3);
							MR2i = 1./60*(-3*u2_i1+27*u2+47*u2i1-13*u2i2+2*u2i3);
							MR3i = 1./60*(-3*u3_i1+27*u3+47*u3i1-13*u3i2+2*u3i3);
							MR4i = 1./60*(-3*u4_i1+27*u4+47*u4i1-13*u4i2+2*u4i3);
							MR5i = 1./60*(-3*u5_i1+27*u5+47*u5i1-13*u5i2+2*u5i3);

						}
						
						
						
						
						
					}
					else {

						ML1 = u1_i1;
						ML2 = u2_i1;
						ML3 = u3_i1;
						ML4 = u4_i1;
						ML5 = u5_i1;
						
						MR1 = u1;
						MR2 = u2;
						MR3 = u3;
						MR4 = u4;
						MR5 = u5;
						
						
						ML1i = u1;
						ML2i = u2;
						ML3i = u3;
						ML4i = u4;
						ML5i = u5;
						
						MR1i = u1i1;
						MR2i = u2i1;
						MR3i = u3i1;
						MR4i = u4i1;
						MR5i = u5i1;

					}

					// ----------------------------- MUSCL-X ----------------------------- //
					// ------------------------------------------------------------------- //




					// ------------------------------------------------------------------ //
					// ----------------------------- Flux-X ----------------------------- //



					/* -------------------------------- */
					/* ----------- Backward ----------- */

					/* jump dU */
					dU1 = MR1-ML1;
					dU2 = MR2-ML2;
					dU3 = MR3-ML3;
					dU4 = MR4-ML4;
					dU5 = MR5-ML5;

					/* lefr parameter */
					_rho = ML1;
					_U = ML2/_rho;
					_V = ML3/_rho;
					_W = ML4/_rho;
					_VV = _U*_U+_V*_V+_W*_W;
					_P = (ML5-0.5*_rho*_VV)*(K-1);
					_T = _P/_rho;
					_C = K*_P/_rho;
					_H = 0.5*_VV+_C/(K-1);

					/* right parameter */
					rho_ = MR1;
					U_ = MR2/rho_;
					V_ = MR3/rho_;
					W_ = MR4/rho_;
					VV_ = U_*U_+V_*V_+W_*W_;
					P_ = (MR5-0.5*rho_*VV_)*(K-1);
					T_ = P_/rho_;
					C_ = K*P_/rho_;
					H_ = 0.5*VV_+C_/(K-1);

					/* flux varibale */

					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					
					// temp5 = 1.0;
					// temp6 = 1.0;
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;
					VV = U*U+V*V+W*W;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
					P = rho*C/K;


					beta = max(VV/C,e);    // ---- theda ---- //
					
					temp = 0.5*(1+beta)*U;    // ---- U' ---- //

					S = 0.5*sqrt(4*beta*C+U*U*(1-beta)*(1-beta));   // ---- C' ---- //





					temp1 = (P_-_P)/rho/beta/C;

					temp2 = temp/S*(U_-_U);

					temp3 = 0.5*(1-beta)*U*temp/S;



					deltaU = (S-temp3-beta*fabs(U))*temp1+temp2;

					deltaP = temp/S*(P_-_P)+(S-fabs(U)+temp3)*rho*(U_-_U);

					/* artificial viscosity */
					Fav1 = fabs(U)*dU1+deltaU*rho;
					Fav2 = fabs(U)*dU2+deltaU*rho*U+deltaP;
					Fav3 = fabs(U)*dU3+deltaU*rho*V;
					Fav4 = fabs(U)*dU4+deltaU*rho*W;
					Fav5 = fabs(U)*dU5+deltaU*rho*H+deltaP*U;

					/* inviscid fluxes */
					inFx1 = 0.5*((_rho*_U+rho_*U_)-Ep*Fav1);
					inFx2 = 0.5*((_rho*_U*_U+_P+rho_*U_*U_+P_)-Ep*Fav2);
					inFx3 = 0.5*((_rho*_V*_U+rho_*V_*U_)-Ep*Fav3);
					inFx4 = 0.5*((_rho*_W*_U+rho_*W_*U_)-Ep*Fav4);
					inFx5 = 0.5*((_U*(3.5*_P+0.5*_rho*_VV)+U_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

					
					dPx = 0.5*(_P+P_);
					/* ----------- Backward ----------- */
					/* -------------------------------- */





					/* --------------------------------- */
					/* ------------ Forward ------------ */
					

					/* jump dU */
					dU1 = MR1i-ML1i;
					dU2 = MR2i-ML2i;
					dU3 = MR3i-ML3i;
					dU4 = MR4i-ML4i;
					dU5 = MR5i-ML5i;

					/* lefr parameter */
					_rho = ML1i;
					_U = ML2i/_rho;
					_V = ML3i/_rho;
					_W = ML4i/_rho;
					_VV = _U*_U+_V*_V+_W*_W;
					_P = (ML5i-0.5*_rho*_VV)*(K-1);
					_T = _P/_rho;
					_C = K*_P/_rho;
					_H = 0.5*_VV+_C/(K-1);

					/* right parameter */
					rho_ = MR1i;
					U_ = MR2i/rho_;
					V_ = MR3i/rho_;
					W_ = MR4i/rho_;
					VV_ = U_*U_+V_*V_+W_*W_;
					P_ = (MR5i-0.5*rho_*VV_)*(K-1);
					T_ = P_/rho_;
					C_ = K*P_/rho_;
					H_ = 0.5*VV_+C_/(K-1);

					/* flux varibale */

					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					
					// temp5 = 1.0;
					// temp6 = 1.0;
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;
					VV = U*U+V*V+W*W;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
					P = rho*C/K;


					beta = max(VV/C,e);    // ---- theda ---- //
					
					temp = 0.5*(1+beta)*U;    // ---- U' ---- //

					S = 0.5*sqrt(4*beta*C+U*U*(1-beta)*(1-beta));   // ---- C' ---- //






					temp1 = (P_-_P)/rho/beta/C;

					temp2 = temp/S*(U_-_U);

					temp3 = 0.5*(1-beta)*U*temp/S;



					deltaU = (S-temp3-beta*fabs(U))*temp1+temp2;



					deltaP = temp/S*(P_-_P)+(S-fabs(U)+temp3)*rho*(U_-_U);

					/* artificial viscosity */
					Fav1 = fabs(U)*dU1+deltaU*rho;
					Fav2 = fabs(U)*dU2+deltaU*rho*U+deltaP;
					Fav3 = fabs(U)*dU3+deltaU*rho*V;
					Fav4 = fabs(U)*dU4+deltaU*rho*W;
					Fav5 = fabs(U)*dU5+deltaU*rho*H+deltaP*U;

					/* inviscid fluxes */
					inFx1i = 0.5*((_rho*_U+rho_*U_)-Ep*Fav1);
					inFx2i = 0.5*((_rho*_U*_U+_P+rho_*U_*U_+P_)-Ep*Fav2);
					inFx3i = 0.5*((_rho*_V*_U+rho_*V_*U_)-Ep*Fav3);
					inFx4i = 0.5*((_rho*_W*_U+rho_*W_*U_)-Ep*Fav4);
					inFx5i = 0.5*((_U*(3.5*_P+0.5*_rho*_VV)+U_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

					
					dPxi = 0.5*(_P+P_);
					/* ------------ Forward ------------ */
					/* --------------------------------- */


					Rfx1 = inFx1i - inFx1;
					Rfx2 = inFx2i - inFx2;
					Rfx3 = inFx3i - inFx3;
					Rfx4 = inFx4i - inFx4;
					Rfx5 = inFx5i - inFx5;

					// ----------------------------- Flux-X ----------------------------- //
					// ------------------------------------------------------------------ //







					// -----------------------------------------------------------------//
					// -------------------------- Y-direction --------------------------//

					IF_j1 = FWS[icube][i][j-1][k];
					u1_j1 = U1_[icube][i][j-1][k][0];
					u2_j1 = U1_[icube][i][j-1][k][1];
					u3_j1 = U1_[icube][i][j-1][k][2];
					u4_j1 = U1_[icube][i][j-1][k][3];
					u5_j1 = U1_[icube][i][j-1][k][4];

					IF_j2 = FWS[icube][i][j-2][k];
					u1_j2 = U1_[icube][i][j-2][k][0];
					u2_j2 = U1_[icube][i][j-2][k][1];
					u3_j2 = U1_[icube][i][j-2][k][2];
					u4_j2 = U1_[icube][i][j-2][k][3];
					u5_j2 = U1_[icube][i][j-2][k][4];


					if(j > 2) {

						IF_j3 = FWS[icube][i][j-3][k];
						u1_j3 = U1_[icube][i][j-3][k][0];
						u2_j3 = U1_[icube][i][j-3][k][1];
						u3_j3 = U1_[icube][i][j-3][k][2];
						u4_j3 = U1_[icube][i][j-3][k][3];
						u5_j3 = U1_[icube][i][j-3][k][4];

					}

					IFj1 = FWS[icube][i][j+1][k];
					u1j1 = U1_[icube][i][j+1][k][0];
					u2j1 = U1_[icube][i][j+1][k][1];
					u3j1 = U1_[icube][i][j+1][k][2];
					u4j1 = U1_[icube][i][j+1][k][3];
					u5j1 = U1_[icube][i][j+1][k][4];

					IFj2 = FWS[icube][i][j+2][k];
					u1j2 = U1_[icube][i][j+2][k][0];
					u2j2 = U1_[icube][i][j+2][k][1];
					u3j2 = U1_[icube][i][j+2][k][2];
					u4j2 = U1_[icube][i][j+2][k][3];
					u5j2 = U1_[icube][i][j+2][k][4];

					if(j < ny)   {

						IFj3 = FWS[icube][i][j+3][k];
						u1j3 = U1_[icube][i][j+3][k][0];
						u2j3 = U1_[icube][i][j+3][k][1];
						u3j3 = U1_[icube][i][j+3][k][2];
						u4j3 = U1_[icube][i][j+3][k][3];
						u5j3 = U1_[icube][i][j+3][k][4];

					}

					// -------------------------- Y-direction --------------------------//
					// -----------------------------------------------------------------//


					// ------------------------------------------------------------------- //
					// ----------------------------- MUSCL-Y ----------------------------- //


































































					if( IFj1 == IFLUID  && IF_j1 == IFLUID) { 

						if ( j==2 | IF_j2 != IFLUID ) {

							ML1 = 1./6*(-1*u1_j2+5*u1_j1+2*u1);
							ML2 = 1./6*(-1*u2_j2+5*u2_j1+2*u2);
							ML3 = 1./6*(-1*u3_j2+5*u3_j1+2*u3);
							ML4 = 1./6*(-1*u4_j2+5*u4_j1+2*u4);
							ML5 = 1./6*(-1*u5_j2+5*u5_j1+2*u5);

							MR1 = 1./6*(2*u1_j1+5*u1-1*u1j1);
							MR2 = 1./6*(2*u2_j1+5*u2-1*u2j1);
							MR3 = 1./6*(2*u3_j1+5*u3-1*u3j1);
							MR4 = 1./6*(2*u4_j1+5*u4-1*u4j1);
							MR5 = 1./6*(2*u5_j1+5*u5-1*u5j1);

						}
						else {

							ML1 = 1./60*(2*u1_j3-13*u1_j2+47*u1_j1+27*u1-3*u1j1);
							ML2 = 1./60*(2*u2_j3-13*u2_j2+47*u2_j1+27*u2-3*u2j1);
							ML3 = 1./60*(2*u3_j3-13*u3_j2+47*u3_j1+27*u3-3*u3j1);
							ML4 = 1./60*(2*u4_j3-13*u4_j2+47*u4_j1+27*u4-3*u4j1);
							ML5 = 1./60*(2*u5_j3-13*u5_j2+47*u5_j1+27*u5-3*u5j1);

							MR1 = 1./60*(-3*u1_j2+27*u1_j1+47*u1-13*u1j1+2*u1j2);
							MR2 = 1./60*(-3*u2_j2+27*u2_j1+47*u2-13*u2j1+2*u2j2);
							MR3 = 1./60*(-3*u3_j2+27*u3_j1+47*u3-13*u3j1+2*u3j2);
							MR4 = 1./60*(-3*u4_j2+27*u4_j1+47*u4-13*u4j1+2*u4j2);
							MR5 = 1./60*(-3*u5_j2+27*u5_j1+47*u5-13*u5j1+2*u5j2);

						}
						
						
						
						
						if ( j==nx | IFj2 != IFLUID) {

							ML1i = 1./6*(-1*u1_j1+5*u1+2*u1j1);
							ML2i = 1./6*(-1*u2_j1+5*u2+2*u2j1);
							ML3i = 1./6*(-1*u3_j1+5*u3+2*u3j1);
							ML4i = 1./6*(-1*u4_j1+5*u4+2*u4j1);
							ML5i = 1./6*(-1*u5_j1+5*u5+2*u5j1);

							MR1i = 1./6*(2*u1+5*u1j1-1*u1j2);
							MR2i = 1./6*(2*u2+5*u2j1-1*u2j2);
							MR3i = 1./6*(2*u3+5*u3j1-1*u3j2);
							MR4i = 1./6*(2*u4+5*u4j1-1*u4j2);
							MR5i = 1./6*(2*u5+5*u5j1-1*u5j2);

						}
						else {

							ML1i = 1./60*(2*u1_j2-13*u1_j1+47*u1+27*u1j1-3*u1j2);
							ML2i = 1./60*(2*u2_j2-13*u2_j1+47*u2+27*u2j1-3*u2j2);
							ML3i = 1./60*(2*u3_j2-13*u3_j1+47*u3+27*u3j1-3*u3j2);
							ML4i = 1./60*(2*u4_j2-13*u4_j1+47*u4+27*u4j1-3*u4j2);
							ML5i = 1./60*(2*u5_j2-13*u5_j1+47*u5+27*u5j1-3*u5j2);

							MR1i = 1./60*(-3*u1_j1+27*u1+47*u1j1-13*u1j2+2*u1j3);
							MR2i = 1./60*(-3*u2_j1+27*u2+47*u2j1-13*u2j2+2*u2j3);
							MR3i = 1./60*(-3*u3_j1+27*u3+47*u3j1-13*u3j2+2*u3j3);
							MR4i = 1./60*(-3*u4_j1+27*u4+47*u4j1-13*u4j2+2*u4j3);
							MR5i = 1./60*(-3*u5_j1+27*u5+47*u5j1-13*u5j2+2*u5j3);

						}
						
						
						
						
						
					}
					else {

						ML1 = u1_j1;
						ML2 = u2_j1;
						ML3 = u3_j1;
						ML4 = u4_j1;
						ML5 = u5_j1;
						
						MR1 = u1;
						MR2 = u2;
						MR3 = u3;
						MR4 = u4;
						MR5 = u5;
						
						
						ML1i = u1;
						ML2i = u2;
						ML3i = u3;
						ML4i = u4;
						ML5i = u5;
						
						MR1i = u1j1;
						MR2i = u2j1;
						MR3i = u3j1;
						MR4i = u4j1;
						MR5i = u5j1;

					}


					// ----------------------------- MUSCL-Y ----------------------------- //
					// ------------------------------------------------------------------- //



					// ------------------------------------------------------------------ //
					// ----------------------------- Flux-Y ----------------------------- //

					/* -------------------------------- */
					/* ----------- Backward ----------- */

					/* jump dU */
					dU1 = MR1-ML1;
					dU2 = MR2-ML2;
					dU3 = MR3-ML3;
					dU4 = MR4-ML4;
					dU5 = MR5-ML5;

					/* lefr parameter */
					_rho = ML1;
					_U = ML2/_rho;
					_V = ML3/_rho;
					_W = ML4/_rho;
					_VV = _U*_U+_V*_V+_W*_W;
					_P = (ML5-0.5*_rho*_VV)*(K-1);
					_T = _P/_rho;
					_C = K*_P/_rho;
					_H = 0.5*_VV+_C/(K-1);

					/* right parameter */
					rho_ = MR1;
					U_ = MR2/rho_;
					V_ = MR3/rho_;
					W_ = MR4/rho_;
					VV_ = U_*U_+V_*V_+W_*W_;
					P_ = (MR5-0.5*rho_*VV_)*(K-1);
					T_ = P_/rho_;
					C_ = K*P_/rho_;
					H_ = 0.5*VV_+C_/(K-1);

					/* flux varibale */

					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					
					// temp5 = 1.0;
					// temp6 = 1.0;
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;
					VV = U*U+V*V+W*W;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
					P = rho*C/K;


					beta = max(VV/C,e);    // ---- theda ---- //
					
					temp = 0.5*(1+beta)*V;    // ---- U' ---- //

					S = 0.5*sqrt(4*beta*C+V*V*(1-beta)*(1-beta));   // ---- C' ---- //






					temp1 = (P_-_P)/rho/beta/C;

					temp2 = temp/S*(V_-_V);

					temp3 = 0.5*(1-beta)*V*temp/S;



					deltaU = (S-temp3-beta*fabs(V))*temp1+temp2;



					deltaP = temp/S*(P_-_P)+(S-fabs(V)+temp3)*rho*(V_-_V);

					/* artificial viscosity */
					Fav1 = fabs(V)*dU1+deltaU*rho;
					Fav2 = fabs(V)*dU2+deltaU*rho*U;
					Fav3 = fabs(V)*dU3+deltaU*rho*V+deltaP;
					Fav4 = fabs(V)*dU4+deltaU*rho*W;
					Fav5 = fabs(V)*dU5+deltaU*rho*H+deltaP*V;

					/* inviscid fluxes */
					inFy1 = 0.5*((_rho*_V+rho_*V_-Ep*Fav1));
					inFy2 = 0.5*((_rho*_U*_V+rho_*U_*V_)-Ep*Fav2);
					inFy3 = 0.5*((_rho*_V*_V+_P+rho_*V_*V_+P_)-Ep*Fav3);
					inFy4 = 0.5*((_rho*_W*_V+rho_*W_*V_)-Ep*Fav4);
					inFy5 = 0.5*((_V*(3.5*_P+0.5*_rho*_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

					dPy = 0.5*(_P+P_);


					/* ----------- Backward ----------- */
					/* -------------------------------- */


					/* --------------------------------- */
					/* ------------ Forward ------------ */
					
					/* jump dU */
					dU1 = MR1i-ML1i;
					dU2 = MR2i-ML2i;
					dU3 = MR3i-ML3i;
					dU4 = MR4i-ML4i;
					dU5 = MR5i-ML5i;

					/* lefr parameter */
					_rho = ML1i;
					_U = ML2i/_rho;
					_V = ML3i/_rho;
					_W = ML4i/_rho;
					_VV = _U*_U+_V*_V+_W*_W;
					_P = (ML5i-0.5*_rho*_VV)*(K-1);
					_T = _P/_rho;
					_C = K*_P/_rho;
					_H = 0.5*_VV+_C/(K-1);

					/* right parameter */
					rho_ = MR1i;
					U_ = MR2i/rho_;
					V_ = MR3i/rho_;
					W_ = MR4i/rho_;
					VV_ = U_*U_+V_*V_+W_*W_;
					P_ = (MR5i-0.5*rho_*VV_)*(K-1);
					T_ = P_/rho_;
					C_ = K*P_/rho_;
					H_ = 0.5*VV_+C_/(K-1);

					/* flux varibale */

					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					
					// temp5 = 1.0;
					// temp6 = 1.0;
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;
					VV = U*U+V*V+W*W;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
					P = rho*C/K;



					beta = max(VV/C,e);    // ---- theda ---- //
					
					temp = 0.5*(1+beta)*V;    // ---- U' ---- //

					S = 0.5*sqrt(4*beta*C+V*V*(1-beta)*(1-beta));   // ---- C' ---- //







					temp1 = (P_-_P)/rho/beta/C;

					temp2 = temp/S*(V_-_V);

					temp3 = 0.5*(1-beta)*V*temp/S;



					deltaU = (S-temp3-beta*fabs(V))*temp1+temp2;

					deltaP = temp/S*(P_-_P)+(S-fabs(V)+temp3)*rho*(V_-_V);


					/* artificial viscosity */
					Fav1 = fabs(V)*dU1+deltaU*rho;
					Fav2 = fabs(V)*dU2+deltaU*rho*U;
					Fav3 = fabs(V)*dU3+deltaU*rho*V+deltaP;
					Fav4 = fabs(V)*dU4+deltaU*rho*W;
					Fav5 = fabs(V)*dU5+deltaU*rho*H+deltaP*V;

					/* inviscid fluxes */
					inFy1i = 0.5*((_rho*_V+rho_*V_-Ep*Fav1));
					inFy2i = 0.5*((_rho*_U*_V+rho_*U_*V_)-Ep*Fav2);
					inFy3i = 0.5*((_rho*_V*_V+_P+rho_*V_*V_+P_)-Ep*Fav3);
					inFy4i = 0.5*((_rho*_W*_V+rho_*W_*V_)-Ep*Fav4);
					inFy5i = 0.5*((_V*(3.5*_P+0.5*_rho*_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);


					dPyi = 0.5*(_P+P_);
					/* ------------ Forward ------------ */
					/* --------------------------------- */



					Rfy1 = inFy1i - inFy1;
					Rfy2 = inFy2i - inFy2;
					Rfy3 = inFy3i - inFy3;
					Rfy4 = inFy4i - inFy4;
					Rfy5 = inFy5i - inFy5;

					// ----------------------------- Flux-Y ----------------------------- //
					// ------------------------------------------------------------------ //


					Rf1 = -(Rfx1/dx+Rfy1/dy+Rfz1/dz);
					Rf2 = -(Rfx2/dx+Rfy2/dy+Rfz2/dz);
					Rf3 = -(Rfx3/dx+Rfy3/dy+Rfz3/dz);
					Rf4 = -(Rfx4/dx+Rfy4/dy+Rfz4/dz);
					Rf5 = -(Rfx5/dx+Rfy5/dy+Rfz5/dz);







					// -------------------------------------------------------------------------//
					// ----------------------------- Viscous_term ----------------------------- //




					// ==== -1 1 0 ==== //
					irhoj = U1_[icube][i-1][j+1][k][0];
					iUj = U1_[icube][i-1][j+1][k][1]/irhoj;
					iVj = U1_[icube][i-1][j+1][k][2]/irhoj;
					iWj = U1_[icube][i-1][j+1][k][3]/irhoj;
					iPj = (U1_[icube][i-1][j+1][k][4]-0.5*irhoj*(iUj*iUj+iVj*iVj+iWj*iWj))*(K-1);
					iTj = iPj/(irhoj*R);

					// ==== -1 0 1 ==== //
					irhok = U1_[icube][i-1][j][k+1][0];
					iUk = U1_[icube][i-1][j][k+1][1]/irhok;
					iVk = U1_[icube][i-1][j][k+1][2]/irhok;
					iWk = U1_[icube][i-1][j][k+1][3]/irhok;
					iPk = (U1_[icube][i-1][j][k+1][4]-0.5*irhok*(iUk*iUk+iVk*iVk+iWk*iWk))*(K-1);
					iTk = iPk/(irhok*R);

					// ==== 1 -1 0 ==== //
					jrhoi = U1_[icube][i+1][j-1][k][0];
					jUi = U1_[icube][i+1][j-1][k][1]/jrhoi;
					jVi = U1_[icube][i+1][j-1][k][2]/jrhoi;
					jWi = U1_[icube][i+1][j-1][k][3]/jrhoi;
					jPi = (U1_[icube][i+1][j-1][k][4]-0.5*jrhoi*(jUi*jUi+jVi*jVi+jWi*jWi))*(K-1);
					jTi = jPi/(jrhoi*R);

					// ==== 1 0 -1 ==== //
					krhoi = U1_[icube][i+1][j][k-1][0];
					kUi = U1_[icube][i+1][j][k-1][1]/krhoi;
					kVi = U1_[icube][i+1][j][k-1][2]/krhoi;
					kWi = U1_[icube][i+1][j][k-1][3]/krhoi;
					kPi = (U1_[icube][i+1][j][k-1][4]-0.5*krhoi*(kUi*kUi+kVi*kVi+kWi*kWi))*(K-1);
					kTi = kPi/(krhoi*R);


					// ==== 0 -1 1 ==== //
					jrhok = U1_[icube][i][j-1][k+1][0];
					jUk = U1_[icube][i][j-1][k+1][1]/jrhok;
					jVk = U1_[icube][i][j-1][k+1][2]/jrhok;
					jWk = U1_[icube][i][j-1][k+1][3]/jrhok;
					jPk = (U1_[icube][i][j-1][k+1][4]-0.5*jrhok*(jUk*jUk+jVk*jVk+jWk*jWk))*(K-1);
					jTk = jPk/(jrhok*R);


					// ==== 0 1 -1 ==== //
					krhoj = U1_[icube][i][j+1][k-1][0];
					kUj = U1_[icube][i][j+1][k-1][1]/krhoj;
					kVj = U1_[icube][i][j+1][k-1][2]/krhoj;
					kWj = U1_[icube][i][j+1][k-1][3]/krhoj;
					kPj = (U1_[icube][i][j+1][k-1][4]-0.5*krhoj*(kUj*kUj+kVj*kVj+kWj*kWj))*(K-1);
					kTj = kPj/(krhoj*R);


					// ==== -1 -1 0 ==== //
					ijrho = U1_[icube][i-1][j-1][k][0];
					ijU = U1_[icube][i-1][j-1][k][1]/ijrho;
					ijV = U1_[icube][i-1][j-1][k][2]/ijrho;
					ijW = U1_[icube][i-1][j-1][k][3]/ijrho;
					ijP = (U1_[icube][i-1][j-1][k][4]-0.5*ijrho*(ijU*ijU+ijV*ijV+ijW*ijW))*(K-1);
					ijT = ijP/(ijrho*R);


					// ==== 0 -1 -1 ==== //
					jkrho = U1_[icube][i][j-1][k-1][0];
					jkU = U1_[icube][i][j-1][k-1][1]/jkrho;
					jkV = U1_[icube][i][j-1][k-1][2]/jkrho;
					jkW = U1_[icube][i][j-1][k-1][3]/jkrho;
					jkP = (U1_[icube][i][j-1][k-1][4]-0.5*jkrho*(jkU*jkU+jkV*jkV+jkW*jkW))*(K-1);
					jkT = jkP/(jkrho*R);


					// ==== -1 0 -1 ==== //
					ikrho = U1_[icube][i-1][j][k-1][0];
					ikU = U1_[icube][i-1][j][k-1][1]/ikrho;
					ikV = U1_[icube][i-1][j][k-1][2]/ikrho;
					ikW = U1_[icube][i-1][j][k-1][3]/ikrho;
					ikP = (U1_[icube][i-1][j][k-1][4]-0.5*ikrho*(ikU*ikU+ikV*ikV+ikW*ikW))*(K-1);
					ikT = ikP/(ikrho*R);


					// ==== 1 1 0 ==== //
					rhoij = U1_[icube][i+1][j+1][k][0];
					Uij = U1_[icube][i+1][j+1][k][1]/rhoij;
					Vij = U1_[icube][i+1][j+1][k][2]/rhoij;
					Wij = U1_[icube][i+1][j+1][k][3]/rhoij;
					Pij = (U1_[icube][i+1][j+1][k][4]-0.5*rhoij*(Uij*Uij+Vij*Vij+Wij*Wij))*(K-1);
					Tij = Pij/(rhoij*R);


					// ==== 1 0 1 ==== //
					rhoik = U1_[icube][i+1][j][k+1][0];
					Uik = U1_[icube][i+1][j][k+1][1]/rhoik;
					Vik = U1_[icube][i+1][j][k+1][2]/rhoik;
					Wik = U1_[icube][i+1][j][k+1][3]/rhoik;
					Pik = (U1_[icube][i+1][j][k+1][4]-0.5*rhoik*(Uik*Uik+Vik*Vik+Wik*Wik))*(K-1);
					Tik = Pik/(rhoik*R);


					// ==== 0 1 1 ==== //
					rhojk = U1_[icube][i][j+1][k+1][0];
					Ujk = U1_[icube][i][j+1][k+1][1]/rhojk;
					Vjk = U1_[icube][i][j+1][k+1][2]/rhojk;
					Wjk = U1_[icube][i][j+1][k+1][3]/rhojk;
					Pjk = (U1_[icube][i][j+1][k+1][4]-0.5*rhojk*(Ujk*Ujk+Vjk*Vjk+Wjk*Wjk))*(K-1);
					Tjk = Pjk/(rhojk*R);



					rho = u1;
					U = u2/u1;
					V = u3/u1;
					W = u4/u1;
					VV = U*U+V*V+W*W;
					P = (u5-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho/R;
					H = 0.5*VV+C/(K-1);



					// ==== -1 0 0 ==== //
					irho = u1_i1;
					iU = u2_i1/irho;
					iV = u3_i1/irho;
					iW = u4_i1/irho;
					iP = (u5_i1-0.5*irho*(iU*iU+iV*iV+iW*iW))*(K-1);
					iT = iP/(irho*R);
					
					
					
					// ==== 0 -1 0 ==== //
					jrho = u1_j1;
					jU = u2_j1/jrho;
					jV = u3_j1/jrho;
					jW = u4_j1/jrho;
					jP = (u5_j1-0.5*jrho*(jU*jU+jV*jV+jW*jW))*(K-1);
					jT = jP/(jrho*R);

					// ==== 0 0 -1 ==== //
					krho = u1_k1;
					kU = u2_k1/krho;
					kV = u3_k1/krho;
					kW = u4_k1/krho;
					kP = (u5_k1-0.5*krho*(kU*kU+kV*kV+kW*kW))*(K-1);
					kT = kP/(krho*R);


					// ==== 1 0 0 ==== //
					rhoi = u1i1;
					Ui = u2i1/rhoi;
					Vi = u3i1/rhoi;
					Wi = u4i1/rhoi;
					Pi = (u5i1-0.5*rhoi*(Ui*Ui+Vi*Vi+Wi*Wi))*(K-1);
					Ti = Pi/(rhoi*R);

					// ==== 0 1 0 ==== //
					rhoj = u1j1;
					Uj = u2j1/rhoj;
					Vj = u3j1/rhoj;
					Wj = u4j1/rhoj;
					Pj = (u5j1-0.5*rhoj*(Uj*Uj+Vj*Vj+Wj*Wj))*(K-1);
					Tj = Pj/(rhoj*R);

					// ==== 0 0 1  ==== //
					rhok = u1k1;
					Uk = u2k1/rhok;
					Vk = u3k1/rhok;
					Wk = u4k1/rhok;
					Pk = (u5k1-0.5*rhok*(Uk*Uk+Vk*Vk+Wk*Wk))*(K-1);
					Tk = Pk/(rhok*R);








						
					
					//--------------------------//
					//---- X viscous fluxes ----//



					/* -------------------------------- */
					/* ----------- Backward ----------- */
					
					
					/* viscous*/
					Tx = 0.5*(T+iT);

					temp = Tx/T0;

					mu_E = mu_L*temp*sqrt(temp)*(T0+110.)/(Tx+110.);

					Pr_E = Pr_L;


					/* X-direction */
					du_dx = (U-iU)*invXI;
					dv_dx = (V-iV)*invXI;
					dw_dx = (W-iW)*invXI;

					dT_dx = (T-iT)*invXI;


					/* Y-direction */
					du_dy = (Uj+iUj-jU-ijU)/4*invET;
					dv_dy = (Vj+iVj-jV-ijV)/4*invET;
					dw_dy = (Wj+iWj-jW-ijW)/4*invET;

					dT_dy = (Tj+iTj-jT-ijT)/4*invET;


					/* Z-direction */
					du_dz = (Uk+iUk-kU-ikU)/4*invZT;
					dv_dz = (Vk+iVk-kV-ikV)/4*invZT;
					dw_dz = (Wk+iWk-kW-ikW)/4*invZT;

					dT_dz = (Tk+iTk-kT-ikT)/4*invZT;


					Ux = 0.5*(U+iU);
					Vx = 0.5*(V+iV);
					Wx = 0.5*(W+iW);

					sigma_xx = mu_E*(2*du_dx-2.0/3.0*(du_dx+dv_dy+dw_dz));
					sigma_xy = mu_E*(du_dy+dv_dx);
					sigma_xz = mu_E*(du_dz+dw_dx);

					LL2 = sigma_xx;

					LL3 = sigma_xy;

					LL4 = sigma_xz;

					LL5 = Ux*sigma_xx+Vx*sigma_xy+Wx*sigma_xz+mu_E*Cv*K*dT_dx/(Pr_E);

					


					/* ----------- Backward ----------- */
					/* -------------------------------- */


					/* --------------------------------- */
					/* ------------ Forward ------------ */
					
					
					/* viscous*/
					Tx = 0.5*(T+Ti);

					temp = Tx/T0;

					mu_E = mu_L*temp*sqrt(temp)*(T0+110.)/(Tx+110.);

					Pr_E = Pr_L;
					
					
					/* X-direction */
					du_dx = (Ui-U)*invXI;
					dv_dx = (Vi-V)*invXI;
					dw_dx = (Wi-W)*invXI;

					dT_dx = (Ti-T)*invXI;


					/* Y-direction */
					du_dy = (Uj+Uij-jU-jUi)/4*invET;
					dv_dy = (Vj+Vij-jV-jVi)/4*invET;
					dw_dy = (Wj+Wij-jW-jWi)/4*invET;

					dT_dy = (Tj+Tij-jT-jTi)/4*invET;


					/* Z-direction */
					du_dz = (Uk+Uik-kU-kUi)/4*invZT;
					dv_dz = (Vk+Vik-kV-kVi)/4*invZT;
					dw_dz = (Wk+Wik-kW-kWi)/4*invZT;

					dT_dz = (Tk+Tik-kT-kTi)/4*invZT;


					Ux = 0.5*(U+Ui);
					Vx = 0.5*(V+Vi);
					Wx = 0.5*(W+Wi);

					sigma_xx = mu_E*(2*du_dx-2.0/3.0*(du_dx+dv_dy+dw_dz));
					sigma_xy = mu_E*(du_dy+dv_dx);
					sigma_xz = mu_E*(du_dz+dw_dx);

					LL2i = sigma_xx;

					LL3i = sigma_xy;

					LL4i = sigma_xz;

					LL5i = Ux*sigma_xx+Vx*sigma_xy+Wx*sigma_xz+mu_E*Cv*K*dT_dx/(Pr_E);


					/* ------------ Forward ------------ */
					/* --------------------------------- */

					
					//---- X viscous fluxes ----//
					//--------------------------//



					//--------------------------//
					//---- Y viscous fluxes ----//

					/* viscous*/
					Ty = 0.5*(T+jT);

					temp = Ty/T0;

					mu_E = mu_L*temp*sqrt(temp)*(T0+110.)/(Ty+110.);

					Pr_E = Pr_L;
					

					/* -------------------------------- */
					/* ----------- Backward ----------- */

					/* X-direction */
					du_dx = (Ui+jUi-iU-ijU)/4*invXI;
					dv_dx = (Vi+jVi-iV-ijV)/4*invXI;
					dw_dx = (Wi+jWi-iW-ijW)/4*invXI;

					dT_dx = (Ti+jTi-iT-ijT)/4*invXI;


					/* Y-direction */
					du_dy = (U-jU)*invET;
					dv_dy = (V-jV)*invET;
					dw_dy = (W-jW)*invET;

					dT_dy = (T-jT)*invET;


					/* Z-direction */
					du_dz = (Uk+jUk-kU-jkU)/4*invZT;
					dv_dz = (Vk+jVk-kV-jkV)/4*invZT;
					dw_dz = (Wk+jWk-kW-jkW)/4*invZT;

					dT_dz = (Tk+jTk-kT-jkT)/4*invZT;

					Uy = 0.5*(U+jU);
					Vy = 0.5*(V+jV);
					Wy = 0.5*(W+jW);

					sigma_yx = mu_E*(du_dy + dv_dx);
					sigma_yy = mu_E*(2*dv_dy-2.0/3.0*(du_dx+dv_dy+dw_dz));
					sigma_yz = mu_E*(dv_dz+dw_dy);

					ML2 = sigma_yx;

					ML3 = sigma_yy;

					ML4 = sigma_yz;

					ML5 = Uy*sigma_yx+Vy*sigma_yy+Wy*sigma_yz+mu_E*Cv*K*dT_dy/(Pr_E);

					/* ----------- Backward ----------- */
					/* -------------------------------- */



					/* --------------------------------- */
					/* ------------ Forward ------------ */

					Ty = 0.5*(T+Tj);

					temp = Ty/T0;

					mu_E = mu_L*temp*sqrt(temp)*(T0+110.)/(Ty+110.);

					Pr_E = Pr_L;

					/* X-direction */
					du_dx = (Ui+Uij-iU-iUj)/4*invXI;
					dv_dx = (Vi+Vij-iV-iVj)/4*invXI;
					dw_dx = (Wi+Wij-iW-iWj)/4*invXI;

					dT_dx = (Ti+Tij-iT-iTj)/4*invXI;


					/* Y-direction */
					du_dy = (Uj-U)*invET;
					dv_dy = (Vj-V)*invET;
					dw_dy = (Wj-W)*invET;

					dT_dy = (Tj-T)*invET;


					/* Z-direction */
					du_dz = (Uk+Ujk-kU-kUj)/4*invZT;
					dv_dz = (Vk+Vjk-kV-kVj)/4*invZT;
					dw_dz = (Wk+Wjk-kW-kWj)/4*invZT;

					dT_dz = (Tk+Tjk-kT-kTj)/4*invZT;

					Uy = 0.5*(U+Uj);
					Vy = 0.5*(V+Vj);
					Wy = 0.5*(W+Wj);

					sigma_yx = mu_E*(du_dy + dv_dx);
					sigma_yy = mu_E*(2*dv_dy-2.0/3.0*(du_dx+dv_dy+dw_dz));
					sigma_yz = mu_E*(dv_dz+dw_dy);

					ML2j = sigma_yx;

					ML3j = sigma_yy;

					ML4j = sigma_yz;

					ML5j = Uy*sigma_yx+Vy*sigma_yy+Wy*sigma_yz+mu_E*Cv*K*dT_dy/(Pr_E);


					/* ------------ Forward ------------ */
					/* --------------------------------- */

					
					//---- Y viscous fluxes ----//
					//--------------------------//



					//--------------------------//
					//---- Z viscous fluxes ----//

					Tz = 0.5*(T+kT);

					temp = Tz/T0;

					mu_E = mu_L*temp*sqrt(temp)*(T0+110.)/(Tz+110.);

					Pr_E = Pr_L;

					/* -------------------------------- */
					/* ----------- Backward ----------- */

					/* X-direction */
					du_dx = (Ui+kUi-iU-ikU)/4*invXI;
					dv_dx = (Vi+kVi-iV-ikV)/4*invXI;
					dw_dx = (Wi+kWi-iW-ikW)/4*invXI;

					dT_dx = (Ti+kTi-iT-ikT)/4*invXI;


					/* Y-direction */
					du_dy = (Uj+kUj-jU-jkU)/4*invET;
					dv_dy = (Vj+kVj-jV-jkV)/4*invET;
					dw_dy = (Wj+kWj-jW-jkW)/4*invET;

					dT_dy = (Tj+kTj-jT-jkT)/4*invET;


					/* Z-direction */
					du_dz = (U-kU)*invZT;
					dv_dz = (V-kV)*invZT;
					dw_dz = (W-kW)*invZT;

					dT_dz = (T-kT)*invZT;

					Uz = 0.5*(U+kU);
					Vz = 0.5*(V+kV);
					Wz = 0.5*(W+kW);

					sigma_zx = mu_E*(du_dz + dw_dx);
					sigma_zy = mu_E*(dv_dz + dw_dy);
					sigma_zz = mu_E*(2*dw_dz-2.0/3.0*(du_dx+dv_dy+dw_dz));


					NL2 = sigma_zx;

					NL3 = sigma_zy;;

					NL4 = sigma_zz;

					NL5 = Uz*sigma_zx+Vz*sigma_zy+Wz*sigma_zz+mu_E*Cv*K*dT_dz/(Pr_E);

					/* ----------- Backward ----------- */
					/* -------------------------------- */




					/* --------------------------------- */
					/* ------------ Forward ------------ */

					Tz = 0.5*(T+Tk);

					temp = Tz/T0;

					mu_E = mu_L*temp*sqrt(temp)*(T0+110.)/(Tz+110.);

					Pr_E = Pr_L;

					/* X-direction */
					du_dx = (Ui+Uik-iU-iUk)/4*invXI;
					dv_dx = (Vi+Vik-iV-iVk)/4*invXI;
					dw_dx = (Wi+Wik-iW-iWk)/4*invXI;

					dT_dx = (Ti+Tik-iT-iTk)/4*invXI;


					/* Y-direction */
					du_dy = (Uj+Ujk-jU-jUk)/4*invET;
					dv_dy = (Vj+Vjk-jV-jVk)/4*invET;
					dw_dy = (Wj+Wjk-jW-jWk)/4*invET;

					dT_dy = (Tj+Tjk-jT-jTk)/4*invET;


					/* Z-direction */
					du_dz = (Uk-U)*invZT;
					dv_dz = (Vk-V)*invZT;
					dw_dz = (Wk-W)*invZT;

					dT_dz = (Tk-T)*invZT;

					Uz = 0.5*(U+Uk);
					Vz = 0.5*(V+Vk);
					Wz = 0.5*(W+Wk);

					sigma_zx = mu_E*(du_dz + dw_dx);
					sigma_zy = mu_E*(dv_dz + dw_dy);
					sigma_zz = mu_E*(2*dw_dz-2.0/3.0*(du_dx+dv_dy+dw_dz));


					NL2k = sigma_zx;

					NL3k = sigma_zy;;

					NL4k = sigma_zz;

					NL5k = Uz*sigma_zx+Vz*sigma_zy+Wz*sigma_zz+mu_E*Cv*K*dT_dz/(Pr_E);

					/* ------------ Forward ------------ */
					/* --------------------------------- */

					
					//---- Z viscous fluxes ----//
					//--------------------------//



					vF2 = (LL2i-LL2)*invXI+(ML2j-ML2)*invET+(NL2k-NL2)*invZT;
					vF3 = (LL3i-LL3)*invXI+(ML3j-ML3)*invET+(NL3k-NL3)*invZT;
					vF4 = (LL4i-LL4)*invXI+(ML4j-ML4)*invET+(NL4k-NL4)*invZT;
					vF5 = (LL5i-LL5)*invXI+(ML5j-ML5)*invET+(NL5k-NL5)*invZT;

					
					// ----------------------------- Viscous_term ----------------------------- //
					// -------------------------------------------------------------------------//


					Rk1 = Rp1+Rf1;
					Rk2 = Rp2+Rf2+vF2;
					Rk3 = Rp3+Rf3+vF3;
					Rk4 = Rp4+Rf4+vF4;
					Rk5 = Rp5+Rf5+vF5;

					// ------------------------------------------------------------------------------------------- //
					// --------------------------------------- Runge-Kutta --------------------------------------- //

					rho = u1;
					U = u2/u1;
					V = u3/u1;
					W = u4/u1;
					VV = U*U+V*V+W*W;
					P = (u5-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho;
					H = 0.5*VV+C/(K-1);

					/* preconditioning */

					beta = max(VV/C,e);


					/* M*inverse(gamma+3*deltaTau/(2*deltaT)*M) */
					temp = K*T*(2*deltaT+3*deltaTau)*(2*deltaT+3*beta*deltaTau);
					temp2 = (K-1)*(beta-1)*deltaT*deltaT;
					temp3 = H-H*K-VV+K*(T+VV);

					d11 = (2*deltaT*(-2*(-(K-1)*(H-VV)-temp3*beta)*deltaT+3*K*T*beta*deltaTau))/temp;
					d12 = (-4*U*temp2)/temp;
					d13 = (-4*V*temp2)/temp;
					d14 = (-4*W*temp2)/temp;
					d15 = (4*temp2)/temp;

					d21 = (4*U*temp3*(beta-1)*deltaT*deltaT)/temp;
					d22 = (2*deltaT*(2*K*(T-U*U*(beta-1))*deltaT+2*U*U*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
					d23 = (-4*U*V*temp2)/temp;
					d24 = (-4*U*W*temp2)/temp;
					d25 = (4*U*temp2)/temp;

					d31 = (4*V*temp3*(beta-1)*deltaT*deltaT)/temp;
					d32 = (-4*U*V*temp2)/temp;
					d33 = (2*deltaT*(2*K*(T-V*V*(beta-1))*deltaT+2*V*V*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
					d34 = (-4*V*W*temp2)/temp;
					d35 = (4*V*temp2)/temp;

					d41 = (4*W*temp3*(beta-1)*deltaT*deltaT)/temp;
					d42 = (-4*U*W*temp2)/temp;
					d43 = (-4*V*W*temp2)/temp;
					d44 = (2*deltaT*(2*K*(T-W*W*(beta-1))*deltaT+2*W*W*(beta-1)*deltaT+3*K*T*beta*deltaTau))/temp;
					d45 = (4*W*temp2)/temp;

					d51 = (4*H*temp3*(beta-1)*deltaT*deltaT)/temp;
					d52 = (-4*H*U*temp2)/temp;
					d53 = (-4*H*V*temp2)/temp;
					d54 = (-4*H*W*temp2)/temp;
					d55 = (2*deltaT*(2*H*(K-1)*(beta-1)*deltaT+K*T*(2*deltaT+3*beta*deltaTau)))/temp;


					if (IF == IFLUID) {

						MR1 = deltaTau*(d11*Rk1+d12*Rk2+d13*Rk3+d14*Rk4+d15*Rk5);
						MR2 = deltaTau*(d21*Rk1+d22*Rk2+d23*Rk3+d24*Rk4+d25*Rk5);
						MR3 = deltaTau*(d31*Rk1+d32*Rk2+d33*Rk3+d34*Rk4+d35*Rk5);
						MR4 = deltaTau*(d41*Rk1+d42*Rk2+d43*Rk3+d44*Rk4+d45*Rk5);
						MR5 = deltaTau*(d51*Rk1+d52*Rk2+d53*Rk3+d54*Rk4+d55*Rk5);

					}
					else {

						MR1 = 0;
						MR2 = 0;
						MR3 = 0;
						MR4 = 0;
						MR5 = 0;




					}

					// --------------------------------------- Runge-Kutta --------------------------------------- //
					// ------------------------------------------------------------------------------------------- //

					Rku1[icube][i][j][k][0] = MR1;
					Rku1[icube][i][j][k][1] = MR2;
					Rku1[icube][i][j][k][2] = MR3;
					Rku1[icube][i][j][k][3] = MR4;
					Rku1[icube][i][j][k][4] = MR5;
					
					// if (csl[icube] == 0 | csl[icube] == 1 | csl[icube] == 2) {
					
						// Residual1[icube][i][j][k][0] = -(Rf2)*dx*dy*dz;
						// Residual1[icube][i][j][k][1] = -(vF2)*dx*dy*dz;
						
					// }
					
					
					// Residual1[icube][i][j][k][1] = -(Rp2+Rf2+vF2)*dx*dy*dz;
					
					// if (IF == IFLUID) {
					
						// Residual1[icube][i][j][k][0] = -(Rp2+Rf2+vF2)*dx*dy*dz;
						
					// }
					
					//if (IF < IFLUID) {
						
					if (IF < IFLUID) Residual1[icube][i][j][k][0] = -( -(dPxi-dPx)/dx + vF2 )*dx*dy*dz;
					

					if (IF < IFLUID) Residual1[icube][i][j][k][1] = -( -(dPzi-dPz)/dz + vF4 )*dx*dy*dz;
					

					//}
					
					
				}
				
				

			}
		}

	}





	if(RK == 1) {

#pragma omp parallel for private(i,j,k)
		for (icube = 1; icube < ncube; icube++) {  
			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {  

						if ( FWS[icube][i][j][k] == IFLUID ) {

							// U1_[icube][i][j][k][0] = U1_[icube][i][j][k][0]+Rku1[icube][i][j][k][0];
							// U1_[icube][i][j][k][1] = U1_[icube][i][j][k][1]+Rku1[icube][i][j][k][1];
							// U1_[icube][i][j][k][2] = U1_[icube][i][j][k][2]+Rku1[icube][i][j][k][2];
							// U1_[icube][i][j][k][3] = U1_[icube][i][j][k][3]+Rku1[icube][i][j][k][3];
							// U1_[icube][i][j][k][4] = U1_[icube][i][j][k][4]+Rku1[icube][i][j][k][4];
							
							U1p1[icube][i][j][k][0] = U1_[icube][i][j][k][0] = U1_[icube][i][j][k][0]+Rku1[icube][i][j][k][0];
							U1p1[icube][i][j][k][1] = U1_[icube][i][j][k][1] = U1_[icube][i][j][k][1]+Rku1[icube][i][j][k][1];
							U1p1[icube][i][j][k][2] = U1_[icube][i][j][k][2] = U1_[icube][i][j][k][2]+Rku1[icube][i][j][k][2];
							U1p1[icube][i][j][k][3] = U1_[icube][i][j][k][3] = U1_[icube][i][j][k][3]+Rku1[icube][i][j][k][3];
							U1p1[icube][i][j][k][4] = U1_[icube][i][j][k][4] = U1_[icube][i][j][k][4]+Rku1[icube][i][j][k][4];


						}
						// else {

						
							// U1_[icube][i][j][k][0] = rho0;
							// U1_[icube][i][j][k][1] = 0;
							// U1_[icube][i][j][k][2] = 0;
							// U1_[icube][i][j][k][3] = 0;
							// U1_[icube][i][j][k][4] = P0/(K-1);
						
						// }


						// U1p1[icube][i][j][k][0] = Rku1[icube][i][j][k][0];
						// U1p1[icube][i][j][k][1] = Rku1[icube][i][j][k][1];
						// U1p1[icube][i][j][k][2] = Rku1[icube][i][j][k][2];
						// U1p1[icube][i][j][k][3] = Rku1[icube][i][j][k][3];
						// U1p1[icube][i][j][k][4] = Rku1[icube][i][j][k][4];
						
					}
				}
			}
		}

	}


	if(RK == 2) {

#pragma omp parallel for private(i,j,k)
		for (icube = 1; icube < ncube; icube++) {  
			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {  

						if ( FWS[icube][i][j][k] == IFLUID ) {
						
							

							// U1_[icube][i][j][k][0] = U1_[icube][i][j][k][0]+A2*B2*U1p1[icube][i][j][k][0]+B2*Rku1[icube][i][j][k][0];
							// U1_[icube][i][j][k][1] = U1_[icube][i][j][k][1]+A2*B2*U1p1[icube][i][j][k][1]+B2*Rku1[icube][i][j][k][1];
							// U1_[icube][i][j][k][2] = U1_[icube][i][j][k][2]+A2*B2*U1p1[icube][i][j][k][2]+B2*Rku1[icube][i][j][k][2];
							// U1_[icube][i][j][k][3] = U1_[icube][i][j][k][3]+A2*B2*U1p1[icube][i][j][k][3]+B2*Rku1[icube][i][j][k][3];
							// U1_[icube][i][j][k][4] = U1_[icube][i][j][k][4]+A2*B2*U1p1[icube][i][j][k][4]+B2*Rku1[icube][i][j][k][4];
							
								U1_[icube][i][j][k][0] = U1p2[icube][i][j][k][0] = 0.75*U1_[icube][i][j][k][0]+0.25*U1p1[icube][i][j][k][0]+0.25*Rku1[icube][i][j][k][0];
								U1_[icube][i][j][k][1] = U1p2[icube][i][j][k][1] = 0.75*U1_[icube][i][j][k][1]+0.25*U1p1[icube][i][j][k][1]+0.25*Rku1[icube][i][j][k][1];
								U1_[icube][i][j][k][2] = U1p2[icube][i][j][k][2] = 0.75*U1_[icube][i][j][k][2]+0.25*U1p1[icube][i][j][k][2]+0.25*Rku1[icube][i][j][k][2];
								U1_[icube][i][j][k][3] = U1p2[icube][i][j][k][3] = 0.75*U1_[icube][i][j][k][3]+0.25*U1p1[icube][i][j][k][3]+0.25*Rku1[icube][i][j][k][3];
								U1_[icube][i][j][k][4] = U1p2[icube][i][j][k][4] = 0.75*U1_[icube][i][j][k][4]+0.25*U1p1[icube][i][j][k][4]+0.25*Rku1[icube][i][j][k][4];
			
							}
							
						// else {

							
							// U1_[icube][i][j][k][0] = rho0;
							// U1_[icube][i][j][k][1] = 0;
							// U1_[icube][i][j][k][2] = 0;
							// U1_[icube][i][j][k][3] = 0;
							// U1_[icube][i][j][k][4] = P0/(K-1);
						
						// }

						
						// U1p1[icube][i][j][k][0] = Rku1[icube][i][j][k][0];
						// U1p1[icube][i][j][k][1] = Rku1[icube][i][j][k][1];
						// U1p1[icube][i][j][k][2] = Rku1[icube][i][j][k][2];
						// U1p1[icube][i][j][k][3] = Rku1[icube][i][j][k][3];
						// U1p1[icube][i][j][k][4] = Rku1[icube][i][j][k][4];


					}
				}
			}
		}

	}


	if(RK == 3) {
	

#pragma omp parallel for private(j,k,rho,P,U,V,W,T,VV,rhoold,Uold,Vold,Wold,VVold,Pold,Told)//reduction(+:e1,e2,e3,e4,e5)

		for (icube = 1; icube < ncube; icube++) {  
			for (i = n_buffer; i < nxx; i++) {
				for (j = n_buffer; j < nyy; j++) {
					for (k = n_buffer; k < nzz; k++) {  

						rhoold = U1_[icube][i][j][k][0];
						Uold = U1_[icube][i][j][k][1]/rhoold;
						Vold = U1_[icube][i][j][k][2]/rhoold;
						Wold = U1_[icube][i][j][k][3]/rhoold;
						VVold = Uold*Uold+Vold*Vold+Wold*Wold;
						Pold = (U1_[icube][i][j][k][4]-0.5*rhoold*VVold)*(K-1);
						Told = Pold/rhoold;

						if ( FWS[icube][i][j][k] == IFLUID ) {
						
							// U1_[icube][i][j][k][0] = U1_[icube][i][j][k][0]+A3*B3*U1_[icube][i][j][k][0]+B3*Rku1[icube][i][j][k][0];
							// U1_[icube][i][j][k][1] = U1_[icube][i][j][k][1]+A3*B3*U1_[icube][i][j][k][1]+B3*Rku1[icube][i][j][k][1];
							// U1_[icube][i][j][k][2] = U1_[icube][i][j][k][2]+A3*B3*U1_[icube][i][j][k][2]+B3*Rku1[icube][i][j][k][2];
							// U1_[icube][i][j][k][3] = U1_[icube][i][j][k][3]+A3*B3*U1_[icube][i][j][k][3]+B3*Rku1[icube][i][j][k][3];
							// U1_[icube][i][j][k][4] = U1_[icube][i][j][k][4]+A3*B3*U1_[icube][i][j][k][4]+B3*Rku1[icube][i][j][k][4];
							
							U1_[icube][i][j][k][0] = 1.0/3.0*U1_[icube][i][j][k][0]+2.0/3.0*U1p2[icube][i][j][k][0]+2.0/3.0*Rku1[icube][i][j][k][0];
							U1_[icube][i][j][k][1] = 1.0/3.0*U1_[icube][i][j][k][1]+2.0/3.0*U1p2[icube][i][j][k][1]+2.0/3.0*Rku1[icube][i][j][k][1];
							U1_[icube][i][j][k][2] = 1.0/3.0*U1_[icube][i][j][k][2]+2.0/3.0*U1p2[icube][i][j][k][2]+2.0/3.0*Rku1[icube][i][j][k][2];
							U1_[icube][i][j][k][3] = 1.0/3.0*U1_[icube][i][j][k][3]+2.0/3.0*U1p2[icube][i][j][k][3]+2.0/3.0*Rku1[icube][i][j][k][3];
							U1_[icube][i][j][k][4] = 1.0/3.0*U1_[icube][i][j][k][4]+2.0/3.0*U1p2[icube][i][j][k][4]+2.0/3.0*Rku1[icube][i][j][k][4];
									
							
							rho = U1_[icube][i][j][k][0];
							U = U1_[icube][i][j][k][1]/rho;
							V = U1_[icube][i][j][k][2]/rho;
							W = U1_[icube][i][j][k][3]/rho;
							VV = U*U+V*V+W*W;
							P = (U1_[icube][i][j][k][4]-0.5*rho*VV)*(K-1);
							T = P/rho;

							
							e1 = e1+(P-Pold)*(P-Pold);
							e2 = e2+(U-Uold)*(U-Uold);
							e3 = e3+(V-Vold)*(V-Vold);	
							e4 = e4+(W-Wold)*(W-Wold);
							e5 = e5+(T-Told)*(T-Told);
							
							// Residual1[icube][i][j][k][2] = (P-Pold)*(P-Pold);
							// Residual1[icube][i][j][k][3] = (U-Uold)*(U-Uold);
							// Residual1[icube][i][j][k][4] = (V-Vold)*(V-Vold);	
							// Residual1[icube][i][j][k][5] = (W-Wold)*(W-Wold);
							// Residual1[icube][i][j][k][6] = (T-Told)*(T-Told);
							
							

						}
							
						else {
						
							// U1_[icube][i][j][k][0] = rho0;
							// U1_[icube][i][j][k][1] = 0;
							// U1_[icube][i][j][k][2] = 0;
							// U1_[icube][i][j][k][3] = 0;
							// U1_[icube][i][j][k][4] = P0/(K-1);
							
							// Residual1[icube][i][j][k][2] = 0;
							// Residual1[icube][i][j][k][3] = 0;
							// Residual1[icube][i][j][k][4] = 0;
							// Residual1[icube][i][j][k][5] = 0;
							// Residual1[icube][i][j][k][6] = 0;
							
						}



					}
				}
			}
		}

		
		
		for (icube = 1; icube < ncube; icube++) {  
			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) { 
					for (k = n_buffer; k <= nz; k++) { 


							e6 = e6+Residual1[icube][i][j][k][0];
							e7 = e7+Residual1[icube][i][j][k][1];
						
						//if ( FWS[icube][i][j][k] > ISOLID ) {
						
							// e1 = e1+Residual1[icube][i][j][k][2];
							// e2 = e2+Residual1[icube][i][j][k][3];
							// e3 = e3+Residual1[icube][i][j][k][4];
							// e4 = e4+Residual1[icube][i][j][k][5];
							// e5 = e5+Residual1[icube][i][j][k][6];
							
						//}
						
					}
				}
			}
		}

		
		double DN = 1./(MPI_Ncube*NcubeX*NcubeY*NcubeZ);


		e1 = sqrt(e1)*DN;
		e2 = sqrt(e2)*DN;
		e3 = sqrt(e3)*DN;
		e4 = sqrt(e4)*DN;
		e5 = sqrt(e5)*DN;
		
		e6 = e6/(0.5*rho0*U0*U0*2.602);
		e7 = e7/(0.5*rho0*U0*U0*2.602);


		MPI_Comm comm;
		comm=MPI_COMM_WORLD;

		MPI_Allreduce ((void*)&e1,(void*)&er[1],1,MPI_DOUBLE,MPI_SUM,comm );
		MPI_Allreduce ((void*)&e2,(void*)&er[2],1,MPI_DOUBLE,MPI_SUM,comm );
		MPI_Allreduce ((void*)&e3,(void*)&er[3],1,MPI_DOUBLE,MPI_SUM,comm );
		MPI_Allreduce ((void*)&e4,(void*)&er[4],1,MPI_DOUBLE,MPI_SUM,comm );
		MPI_Allreduce ((void*)&e5,(void*)&er[5],1,MPI_DOUBLE,MPI_SUM,comm );
		
		
		MPI_Allreduce ((void*)&e6,(void*)&er[6], 1, MPI_DOUBLE, MPI_SUM, comm );
		MPI_Allreduce ((void*)&e7,(void*)&er[7], 1, MPI_DOUBLE, MPI_SUM, comm );


	}


	
}