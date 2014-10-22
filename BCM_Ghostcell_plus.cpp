




#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "Resolution.h"

extern int Ncube;    
extern int NBC_plus;

void BCM_Ghostcell_plus
(
// =============================================================================== //
int myid,

int *NNBC,

double deltaT,
double deltaTau,

double (*weight) = new double[NBC_plus*8+1],
double (*N_dis) = new double[NBC_plus*4+1],
int (*GCindex) = new int[NBC_plus*4+1],
int (*IPsur) = new int[NBC_plus*4+1],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// =============================================================================== //
)

{

#include "BCM.h"
#include "prm.h"

int gicube,gi,gj,gk;
int iNBC, Ntemp;
double Nx, Ny, Nz;
double Ndis;

double rho, U, V, W, VV, P, T;
double rhoold,Uold,Vold,Wold,VVold,Pold,Told;
	

double rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111;
double U000, U100, U010, U001, U110, U101, U011, U111;
double V000, V100, V010, V001, V110, V101, V011, V111;
double W000, W100, W010, W001, W110, W101, W011, W111;
double P000, P100, P010, P001, P110, P101, P011, P111;
double VV000, VV100, VV010, VV001, VV110, VV101, VV011, VV111;
double T000, T100, T010, T001, T110, T101, T011, T111;

double w1, w2, w3, w4, w5, w6, w7, w8;


// ---- GCindex => IPsur => weight ---- //

#pragma omp parallel for private(Ntemp,icube,i,j,k,\
	rho000, rho100, rho010, rho001, rho110, rho101, rho011, rho111,\
	U000, U100, U010, U001, U110, U101, U011, U111,\
	V000, V100, V010, V001, V110, V101, V011, V111,\
	W000, W100, W010, W001, W110, W101, W011, W111,\
	P000, P100, P010, P001, P110, P101, P011, P111,\
	VV000, VV100, VV010, VV001, VV110, VV101, VV011, VV111,\
	T000, T100, T010, T001, T110, T101, T011, T111,\
	w1,w2,w3,w4,w5,w6,w7,w8,\
	gicube,gi,gj,gk,\
	Nx, Ny, Nz, Ndis,\
	rho,P,U,V,W,VV,\
	rhoold,Pold,Uold,Vold,Wold,VVold\
	)
	
	for (iNBC = 1; iNBC <= *NNBC; iNBC++) {

		Ntemp = (iNBC-1)*4;

		icube = IPsur[Ntemp+1];
		i = IPsur[Ntemp+2];
		j = IPsur[Ntemp+3];
		k = IPsur[Ntemp+4];
		
		rho000 = U1_[icube][i][j][k][0];
		U000 = U1_[icube][i][j][k][1]/rho000;
		V000 = U1_[icube][i][j][k][2]/rho000;
		W000 = U1_[icube][i][j][k][3]/rho000;
		VV000 = U000*U000+V000*V000+W000*W000;
		P000 = (U1_[icube][i][j][k][4]-0.5*rho000*VV000)*(K-1);
		T000 = P000/rho000/R;
	
		rho100 = U1_[icube][i+1][j][k][0];
		U100 = U1_[icube][i+1][j][k][1]/rho100;
		V100 = U1_[icube][i+1][j][k][2]/rho100;
		W100 = U1_[icube][i+1][j][k][3]/rho100;
		VV100 = U100*U100+V100*V100+W100*W100;
		P100 = (U1_[icube][i+1][j][k][4]-0.5*rho100*VV100)*(K-1);
		T100 = P100/rho100/R;

		rho010 = U1_[icube][i][j+1][k][0];
		U010 = U1_[icube][i][j+1][k][1]/rho010;
		V010 = U1_[icube][i][j+1][k][2]/rho010;
		W010 = U1_[icube][i][j+1][k][3]/rho010;
		VV010 = U010*U010+V010*V010+W010*W010;
		P010 = (U1_[icube][i][j+1][k][4]-0.5*rho010*VV010)*(K-1);
		T010 = P010/rho010/R;

		rho001 = U1_[icube][i][j][k+1][0];
		U001 = U1_[icube][i][j][k+1][1]/rho001;
		V001 = U1_[icube][i][j][k+1][2]/rho001;
		W001 = U1_[icube][i][j][k+1][3]/rho001;
		VV001 = U001*U001+V001*V001+W001*W001;
		P001 = (U1_[icube][i][j][k+1][4]-0.5*rho001*VV001)*(K-1);
		T001 = P001/rho001/R;

		rho110 = U1_[icube][i+1][j+1][k][0];
		U110 = U1_[icube][i+1][j+1][k][1]/rho110;
		V110 = U1_[icube][i+1][j+1][k][2]/rho110;
		W110 = U1_[icube][i+1][j+1][k][3]/rho110;
		VV110 = U110*U110+V110*V110+W110*W110;
		P110 = (U1_[icube][i+1][j+1][k][4]-0.5*rho110*VV110)*(K-1);
		T110 = P110/rho110/R;

		rho101 = U1_[icube][i+1][j][k+1][0];
		U101 = U1_[icube][i+1][j][k+1][1]/rho101;
		V101 = U1_[icube][i+1][j][k+1][2]/rho101;
		W101 = U1_[icube][i+1][j][k+1][3]/rho101;
		VV101 = U101*U101+V101*V101+W101*W101;
		P101 = (U1_[icube][i+1][j][k+1][4]-0.5*rho101*VV101)*(K-1);
		T101 = P101/rho101/R;

		rho011 = U1_[icube][i][j+1][k+1][0];
		U011 = U1_[icube][i][j+1][k+1][1]/rho011;
		V011 = U1_[icube][i][j+1][k+1][2]/rho011;
		W011 = U1_[icube][i][j+1][k+1][3]/rho011;
		VV011 = U011*U011+V011*V011+W011*W011;
		P011 = (U1_[icube][i][j+1][k+1][4]-0.5*rho011*VV011)*(K-1);
		T011 = P011/rho011/R;


		rho111 = U1_[icube][i+1][j+1][k+1][0];
		U111 = U1_[icube][i+1][j+1][k+1][1]/rho111;
		V111 = U1_[icube][i+1][j+1][k+1][2]/rho111;
		W111 = U1_[icube][i+1][j+1][k+1][3]/rho111;
		VV111 = U111*U111+V111*V111+W111*W111;
		P111 = (U1_[icube][i+1][j+1][k+1][4]-0.5*rho111*VV111)*(K-1);
		T111 = P111/rho111/R;


		Ntemp = (iNBC-1)*8;

		w1 = weight[Ntemp+1];
		w2 = weight[Ntemp+2];
		w3 = weight[Ntemp+3];
		w4 = weight[Ntemp+4];
		w5 = weight[Ntemp+5];
		w6 = weight[Ntemp+6];
		w7 = weight[Ntemp+7];
		w8 = weight[Ntemp+8];



// --------------------------------------------------------------- //
// ----------------- Pressure boundary condition ----------------- //
		
		Ntemp = (iNBC-1)*4;

		gicube = GCindex[Ntemp+1];
		gi = GCindex[Ntemp+2];
		gj = GCindex[Ntemp+3];
		gk = GCindex[Ntemp+4];
		
		Ntemp = (iNBC-1)*4;
		Nx = N_dis[Ntemp+1];
		Ny = N_dis[Ntemp+2];
		Nz = N_dis[Ntemp+3];
		Ndis = N_dis[Ntemp+4];

		rhoold = U1[icube][i][j][k][0];
		Uold = U1[icube][i][j][k][1]/rhoold;
		Vold = U1[icube][i][j][k][2]/rhoold;
		Wold = U1[icube][i][j][k][3]/rhoold;
		VVold = Uold*Uold+Vold*Vold+Wold*Wold;
		Pold = (U1[icube][i][j][k][4]-0.5*rhoold*VVold)*(K-1);
		
		VVold = Uold*Nx+Vold*Ny+Wold*Nz;    // ---- Velocity mag. in normal direction ---- //

// ----------------- Pressure boundary condition ----------------- //
// --------------------------------------------------------------- //		


		
		rho = w1*(rho000-rho0)+w2*(rho100-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0)+rho0;

		U = w1*U000+w2*U100+w3*U010+w4*U001+w5*U110+w6*U101+w7*U011+w8*U111;
		V = w1*V000+w2*V100+w3*V010+w4*V001+w5*V110+w6*V101+w7*V011+w8*V111;
		W = w1*W000+w2*W100+w3*W010+w4*W001+w5*W110+w6*W101+w7*W011+w8*W111;
		
		VV = U*U+V*V+W*W;

		P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;

		U1_[gicube][gi][gj][gk][0] = rho;
		U1_[gicube][gi][gj][gk][1] = 0.5*rho*U;
		U1_[gicube][gi][gj][gk][2] = 0.5*rho*V;
		U1_[gicube][gi][gj][gk][3] = 0.5*rho*W;
		U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*0.25*VV;
		
		
		if ( gi==(i  ) & gj==(j  ) & gk==(k  ) ) { 
		
			if (w1 > 0.9999) {

				rho = rho000;
			}
			else {

				rho = (w2*(rho100-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0))/(1-w1)+rho0;

			}
				
			U = (w2*U100+w3*U010+w4*U001+w5*U110+w6*U101+w7*U011+w8*U111)/(2-w1);
			V = (w2*V100+w3*V010+w4*V001+w5*V110+w6*V101+w7*V011+w8*V111)/(2-w1);
			W = (w2*W100+w3*W010+w4*W001+w5*W110+w6*W101+w7*W011+w8*W111)/(2-w1);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i+1) & gj==(j  ) & gk==(k  ) ) { 
		
			
			if (w2 > 0.9999) {

				rho = rho100;
			}
			else {

				rho = (w1*(rho000-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0))/(1-w2)+rho0;

			}
				
			
			U = (w1*U000+w3*U010+w4*U001+w5*U110+w6*U101+w7*U011+w8*U111)/(2-w2);
			V = (w1*V000+w3*V010+w4*V001+w5*V110+w6*V101+w7*V011+w8*V111)/(2-w2);
			W = (w1*W000+w3*W010+w4*W001+w5*W110+w6*W101+w7*W011+w8*W111)/(2-w2);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i  ) & gj==(j+1) & gk==(k  ) ) { 
		
			if (w3 > 0.9999) {

				rho = rho010;


			}
			else {

				rho = (w1*(rho000-rho0)+w2*(rho100-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0))/(1-w3)+rho0;

			}
				
				
			U = (w1*U000+w2*U100+w4*U001+w5*U110+w6*U101+w7*U011+w8*U111)/(2-w3);
			V = (w1*V000+w2*V100+w4*V001+w5*V110+w6*V101+w7*V011+w8*V111)/(2-w3);
			W = (w1*W000+w2*W100+w4*W001+w5*W110+w6*W101+w7*W011+w8*W111)/(2-w3);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;
				
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i  ) & gj==(j  ) & gk==(k+1) ) { 
		
			if (w4 > 0.9999) {

				rho = rho001;
			}
			else {

				rho = (w1*(rho000-rho0)+w2*(rho100-rho0)+w3*(rho010-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0))/(1-w4)+rho0;

			}
				
			
			U = (w1*U000+w2*U100+w3*U010+w5*U110+w6*U101+w7*U011+w8*U111)/(2-w4);
			V = (w1*V000+w2*V100+w3*V010+w5*V110+w6*V101+w7*V011+w8*V111)/(2-w4);
			W = (w1*W000+w2*W100+w3*W010+w5*W110+w6*W101+w7*W011+w8*W111)/(2-w4);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			
			}
			
		if ( gi==(i+1) & gj==(j+1) & gk==(k  ) ) { 
			
			
			if (w5 > 0.9999) {

				rho = rho110;
			}
			else {

				rho = (w1*(rho000-rho0)+w2*(rho100-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0))/(1-w5)+rho0;

			}
				
			
			
			U = (w1*U000+w2*U100+w3*U010+w4*U001+w6*U101+w7*U011+w8*U111)/(2-w5);
			V = (w1*V000+w2*V100+w3*V010+w4*V001+w6*V101+w7*V011+w8*V111)/(2-w5);
			W = (w1*W000+w2*W100+w3*W010+w4*W001+w6*W101+w7*W011+w8*W111)/(2-w5);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;

			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i+1) & gj==(j  ) & gk==(k+1) ) { 
			
			
			if (w6 > 0.9999) {

				rho = rho101;
			}
			else {

				rho = (w1*(rho000-rho0)+w2*(rho100-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w7*(rho011-rho0)+w8*(rho111-rho0))/(1-w6)+rho0;

			}
			
				
			
			U = (w1*U000+w2*U100+w3*U010+w4*U001+w5*U110+w7*U011+w8*U111)/(2-w6);
			V = (w1*V000+w2*V100+w3*V010+w4*V001+w5*V110+w7*V011+w8*V111)/(2-w6);
			W = (w1*W000+w2*W100+w3*W010+w4*W001+w5*W110+w7*W011+w8*W111)/(2-w6);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i  ) & gj==(j+1) & gk==(k+1) ) { 
			
			if (w7 > 0.9999) {

				rho = rho011;
			}
			else {

				rho = (w1*(rho000-rho0)+w2*(rho100-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w8*(rho111-rho0))/(1-w7)+rho0;

			}
				
				
			U = (w1*U000+w2*U100+w3*U010+w4*U001+w5*U110+w6*U101+w8*U111)/(2-w7);
			V = (w1*V000+w2*V100+w3*V010+w4*V001+w5*V110+w6*V101+w8*V111)/(2-w7);
			W = (w1*W000+w2*W100+w3*W010+w4*W001+w5*W110+w6*W101+w8*W111)/(2-w7);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;
			
			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}
			
		if ( gi==(i+1) & gj==(j+1) & gk==(k+1) ) {
			
			
			if (w8 > 0.9999) {

				rho = rho111;
			}
			else {

				rho = (w1*(rho000-rho0)+w2*(rho100-rho0)+w3*(rho010-rho0)+w4*(rho001-rho0)+w5*(rho110-rho0)+w6*(rho101-rho0)+w7*(rho011-rho0))/(1-w8)+rho0;
				
			}
				
			
			U = (w1*U000+w2*U100+w3*U010+w4*U001+w5*U110+w6*U101+w7*U011)/(2-w8);
			V = (w1*V000+w2*V100+w3*V010+w4*V001+w5*V110+w6*V101+w7*V011)/(2-w8);
			W = (w1*W000+w2*W100+w3*W010+w4*W001+w5*W110+w6*W101+w7*W011)/(2-w8);
			
			VV = U*U+V*V+W*W;

			P = Pold - (rho*U-rhoold*Uold)*Ndis/deltaT-rho*U*U;

			U1_[gicube][gi][gj][gk][0] = rho;
			U1_[gicube][gi][gj][gk][1] = rho*U;
			U1_[gicube][gi][gj][gk][2] = rho*V;
			U1_[gicube][gi][gj][gk][3] = rho*W;
			U1_[gicube][gi][gj][gk][4] = P/(K-1)+0.5*rho*VV;

			}

		
		
	}

}