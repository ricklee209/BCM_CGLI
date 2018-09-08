



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
extern int MPI_Nadj;

extern int Max_nei_eq;
extern int Max_nei_sb;
extern int Max_nei_bs;

extern int Ncpu_eq;
extern int Ncpu_sb;
extern int Ncpu_bs;

#include "implicit.h"



void BCM_Flux_XYZ_Viscous_DPLUSGS
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int mPI_Nadj,

int ncpu_bs,
int ncpu_eq, 
int ncpu_sb,

int max_nei_bs,
int max_nei_eq,
int max_nei_sb,

int nadjX_eq, 
int nadjY_eq, 
int nadjZ_eq,

int nadjX_bs_plus, 
int nadjX_sb_plus,
int nadjX_bs_minus,
int nadjX_sb_minus,

int nadjY_bs_plus, 
int nadjY_sb_plus, 
int nadjY_bs_minus,
int nadjY_sb_minus,

int nadjZ_bs_plus,
int nadjZ_sb_plus,
int nadjZ_bs_minus,
int nadjZ_sb_minus,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

int (*MPI_cpu) = new int[MPI_Nadj],
int (*MPI_cube) = new int[MPI_Nadj],

int (*MPI_cpu_adj) = new int[MPI_Nadj],
int (*MPI_cube_adj) = new int[MPI_Nadj],

int (*MPI_direction) = new int[MPI_Nadj],
int (*MPI_interface) = new int[MPI_Nadj],

int (*neighbor_cpu_eq) = new int[Ncpu_eq],
int (*Ncube_Ncpu_eq) = new int[Ncpu_eq], 

int (*neighbor_cpu_sb) = new int[Ncpu_sb],
int (*Ncube_Ncpu_sb) = new int[Ncpu_sb],

int (*neighbor_cpu_bs) = new int[Ncpu_bs],
int (*Ncube_Ncpu_bs) = new int[Ncpu_bs],

int (*Scube_Ncpu_eq) = new int[Max_nei_eq+1],
int (*Rcube_Ncpu_eq) = new int[Max_nei_eq+1],

double (*send_data_curr_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],
double (*recv_data_curr_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],

double (*send_data_neig_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],
double (*recv_data_neig_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],

int (*Sdir_eq) = new int[Max_nei_eq+1],
int (*Rdir_eq) = new int[Max_nei_eq+1],


int (*Scube_Ncpu_sb) = new int[Max_nei_sb+1],
int (*Rcube_Ncpu_sb) = new int[Max_nei_sb+1],

double (*send_data_curr_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1],
double (*recv_data_curr_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1],

double (*send_data_neig_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1],
double (*recv_data_neig_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1],

int (*Sdir_sb) = new int[Max_nei_sb+1],
int (*Rdir_sb) = new int[Max_nei_sb+1],


int (*Scube_Ncpu_bs) = new int[Max_nei_bs+1],
int (*Rcube_Ncpu_bs) = new int[Max_nei_bs+1],

double (*send_data_curr_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1],
double (*recv_data_curr_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1],

double (*send_data_neig_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1],
double (*recv_data_neig_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1],

int (*Sdir_bs) = new int[Max_nei_bs+1],
int (*Rdir_bs) = new int[Max_nei_bs+1],

int (*ist_eq) = new int[Ncpu_eq],
int (*ist_sb) = new int[Ncpu_sb],
int (*ist_bs) = new int[Ncpu_bs],

int (*adjN_sb) = new int[Max_nei_sb],

int (*RadjN_bs) = new int[Max_nei_bs],
int (*SadjN_bs) = new int[Max_nei_bs],


int (*adj_number)[5][7] = new int[Ncube][5][7],

int (*adjX_eq) = new int[Ncube],
int (*adjY_eq) = new int[Ncube],
int (*adjZ_eq) = new int[Ncube],

int (*adjX_bs_plus) = new int[Ncube],
int (*adjX_sb_plus) = new int[Ncube],
int (*adjX_bs_minus) = new int[Ncube],
int (*adjX_sb_minus) = new int[Ncube],

int (*adjY_bs_plus) = new int[Ncube],
int (*adjY_sb_plus) = new int[Ncube],
int (*adjY_bs_minus) = new int[Ncube],
int (*adjY_sb_minus) = new int[Ncube],

int (*adjZ_bs_plus) = new int[Ncube],
int (*adjZ_sb_plus) = new int[Ncube],
int (*adjZ_bs_minus) = new int[Ncube],
int (*adjZ_sb_minus) = new int[Ncube],


int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

int (*csl) = new int[Ncube],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1p1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1p2)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Roe_dis)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],

double (*CFL_tau)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],

double (*Rku1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Residual1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*er) = new double[10]

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


    double UN_rho,UN_U,UN_V,UN_W,UN_VV,UN_P,UN_T,UN_C,UN_H;
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

    double dTx,dTy,dTz,dTxi,dTyi,dTzi;

    double e1 = 0;
    double e2 = 0;
    double e3 = 0;
    double e4 = 0;
    double e5 = 0;
    
    double e6 = 0;
    double e7 = 0;
    double e8 = 0;
    double e9 = 0;

    double theda_p, U_p, C_p;

    double Cdiss;
    
    int nsweep = 4;
    int isweep;
    
    double Aplus11,Aplus12,Aplus13,Aplus14,
    Aplus21,Aplus22,
    Aplus31,        Aplus33,
    Aplus41,          	   Aplus44,
    Aplus51,Aplus52,Aplus53,Aplus54,Aplus55;

    double Bplus11,Bplus12,Bplus13,Bplus14,
    Bplus21,Bplus22,
    Bplus31,        Bplus33,
    Bplus41,          	   Bplus44,
    Bplus51,Bplus52,Bplus53,Bplus54,Bplus55;

    double Cplus11,Cplus12,Cplus13,Cplus14,
    Cplus21,Cplus22,
    Cplus31,        Cplus33,
    Cplus41,          	   Cplus44,
    Cplus51,Cplus52,Cplus53,Cplus54,Cplus55;

    double Amius11,Amius22,Amius33,Amius44,Amius55,
    Bmius11,Bmius22,Bmius33,Bmius44,Bmius55,
    Cmius11,Cmius22,Cmius33,Cmius44,Cmius55;
    
    double UU1,UU2,UU3,UU4,UU5;
    
    int ib=0;

    
    #pragma omp parallel for private(\
        i,j,k\
        )
    for (icube = 1; icube < ncube; icube++) {  
        
        for (i = 0; i <= nxxx; i++) {
            for (j = 0; j <= nxxx; j++) {
                for (k = 0; k <= nxxx; k++) {
                    
                    U1p1[icube][i][j][k][0] = 0;
                    U1p1[icube][i][j][k][1] = 0;
                    U1p1[icube][i][j][k][2] = 0;
                    U1p1[icube][i][j][k][3] = 0;
                    U1p1[icube][i][j][k][4] = 0;
                    
                    U1p2[icube][i][j][k][0] = 0;
                    U1p2[icube][i][j][k][1] = 0;
                    U1p2[icube][i][j][k][2] = 0;
                    U1p2[icube][i][j][k][3] = 0;
                    U1p2[icube][i][j][k][4] = 0;

                }
            }
        }
        
    }
    
    
    #pragma omp barrier
    

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
        UU1,UU2,UU3,UU4,UU5,\
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
        UN_rho,UN_U,UN_V,UN_W,UN_VV,UN_P,UN_T,UN_C,UN_H,\
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
        dPx,dPy,dPz,dPxi,dPyi,dPzi,\
        dTx,dTy,dTz,dTxi,dTyi,dTzi,\
        theda_p, U_p, C_p,\
        Cdiss,deltaTau,ib\
        ) 


    for (icube = 1; icube < ncube; icube++) {  

        dx = dy = dz = cube_size[icube]/NcubeX;
        invXI = invET = invZT = 1./dx;

        for (i = 2; i <= nx; i++) {
            for (j = 2; j <= ny; j++) {
                for (k = 2; k <= nz; k++) {



                    
                    IF = FWS[icube][i][j][k];
                    u1 = U1_[icube][i][j][k][0];
                    u2 = U1_[icube][i][j][k][1];
                    u3 = U1_[icube][i][j][k][2];
                    u4 = U1_[icube][i][j][k][3];
                    u5 = U1_[icube][i][j][k][4];



                    #ifdef NODT

                    Rp1 = 0.0;
                    Rp2 = 0.0;
                    Rp3 = 0.0;
                    Rp4 = 0.0;
                    Rp5 = 0.0;
                    
                    #else

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

                    
                    Rp1 = -(3*u1-4*uu1+u1q)/(2*deltaT);
                    Rp2 = -(3*u2-4*uu2+u2q)/(2*deltaT);
                    Rp3 = -(3*u3-4*uu3+u3q)/(2*deltaT);
                    Rp4 = -(3*u4-4*uu4+u4q)/(2*deltaT);
                    Rp5 = -(3*u5-4*uu5+u5q)/(2*deltaT);

                    #endif 



                    Cdiss = 1.0;

                    #ifdef ILES
                      Cdiss = Roe_dis[icube][i][j][k];
                    #endif 
                    
                    #ifdef NOUPD
                      Cdiss = 0.0;
                    #endif 



                    // -----------------------------------------------------------------//
                    // -------------------------- Z-direction --------------------------//

                    IF_k1 = FWS[icube][i][j][k-1];
                    IF_k2 = FWS[icube][i][j][k-2];
                    IFk1 = FWS[icube][i][j][k+1];
                    IFk2 = FWS[icube][i][j][k+2];


                    #ifdef Kcomputer
                    //if(k > 2) {
                    #else
                    if(k > 2) {
                        #endif
                        

                        IF_k3 = FWS[icube][i][j][k-3];
                        u1_k3 = U1_[icube][i][j][k-3][0];
                        u2_k3 = U1_[icube][i][j][k-3][1];
                        u3_k3 = U1_[icube][i][j][k-3][2];
                        u4_k3 = U1_[icube][i][j][k-3][3];
                        u5_k3 = U1_[icube][i][j][k-3][4];
                        #ifdef Kcomputer
                        //}
                        #else
                    }
                    #endif


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


                    #ifdef Kcomputer
                    //iif(k < nz) { 
                    #else
                    if(k < nz) { 
                        #endif

                        

                        IFk3 = FWS[icube][i][j][k+3];

                        u1k3 = U1_[icube][i][j][k+3][0];
                        u2k3 = U1_[icube][i][j][k+3][1];
                        u3k3 = U1_[icube][i][j][k+3][2];
                        u4k3 = U1_[icube][i][j][k+3][3];
                        u5k3 = U1_[icube][i][j][k+3][4];

                        #ifdef Kcomputer
                        //}
                        #else
                    }
                    #endif
                    
                    // -------------------------- Z-direction --------------------------//
                    // -----------------------------------------------------------------//




                    // ------------------------------------------------------------------- //
                    // ----------------------------- MUSCL-Z ----------------------------- //

                    #ifdef limiter

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
                            
                            ML1 = u1_k1 + max(0.0, min( min(u1_k1-u1_k2, u1-u1_k1), ML1-u1_k1 ));
                            ML2 = u2_k1 + max(0.0, min( min(u2_k1-u2_k2, u2-u1_k1), ML2-u2_k1 ));
                            ML3 = u3_k1 + max(0.0, min( min(u3_k1-u3_k2, u3-u1_k1), ML3-u3_k1 ));
                            ML4 = u4_k1 + max(0.0, min( min(u4_k1-u4_k2, u4-u1_k1), ML4-u4_k1 ));
                            ML5 = u5_k1 + max(0.0, min( min(u5_k1-u5_k2, u5-u1_k1), ML5-u5_k1 ));
                            
                            MR1 = u1 - max(0.0, min( min(u1k1-u1, u1-u1_k1), MR1-u1 ));
                            MR2 = u2 - max(0.0, min( min(u2k1-u2, u2-u2_k1), MR2-u2 ));
                            MR3 = u3 - max(0.0, min( min(u3k1-u3, u3-u3_k1), MR3-u3 ));
                            MR4 = u4 - max(0.0, min( min(u4k1-u4, u4-u4_k1), MR4-u4 ));
                            MR5 = u5 - max(0.0, min( min(u5k1-u5, u5-u5_k1), MR5-u5 ));
                            
                            
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
                            
                            
                            ML1 = u1_k1 + max(0.0, min( min(u1_k1-u1_k2, u1-u1_k1), ML1-u1_k1 ));
                            ML2 = u2_k1 + max(0.0, min( min(u2_k1-u2_k2, u2-u1_k1), ML2-u2_k1 ));
                            ML3 = u3_k1 + max(0.0, min( min(u3_k1-u3_k2, u3-u1_k1), ML3-u3_k1 ));
                            ML4 = u4_k1 + max(0.0, min( min(u4_k1-u4_k2, u4-u1_k1), ML4-u4_k1 ));
                            ML5 = u5_k1 + max(0.0, min( min(u5_k1-u5_k2, u5-u1_k1), ML5-u5_k1 ));
                            
                            MR1 = u1 - max(0.0, min( min(u1k1-u1, u1-u1_k1), MR1-u1 ));
                            MR2 = u2 - max(0.0, min( min(u2k1-u2, u2-u2_k1), MR2-u2 ));
                            MR3 = u3 - max(0.0, min( min(u3k1-u3, u3-u3_k1), MR3-u3 ));
                            MR4 = u4 - max(0.0, min( min(u4k1-u4, u4-u4_k1), MR4-u4 ));
                            MR5 = u5 - max(0.0, min( min(u5k1-u5, u5-u5_k1), MR5-u5 ));

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
                            
                            
                            ML1i = u1 + max(0.0, min( min(u1k1-u1, u1-u1_k1), MR1-u1 ));
                            ML2i = u2 + max(0.0, min( min(u2k1-u2, u2-u2_k1), MR2-u2 ));
                            ML3i = u3 + max(0.0, min( min(u3k1-u3, u3-u3_k1), MR3-u3 ));
                            ML4i = u4 + max(0.0, min( min(u4k1-u4, u4-u4_k1), MR4-u4 ));
                            ML5i = u5 + max(0.0, min( min(u5k1-u5, u5-u5_k1), MR5-u5 ));
                            
                            MR1i = u1k1 - max(0.0, min( min(u1k2-u1k1, u1k1-u1),MR1i-u1k1 ));
                            MR2i = u2k1 - max(0.0, min( min(u2k2-u2k1, u2k1-u2),MR2i-u2k1 ));
                            MR3i = u3k1 - max(0.0, min( min(u3k2-u3k1, u3k1-u3),MR3i-u3k1 ));
                            MR4i = u4k1 - max(0.0, min( min(u4k2-u4k1, u4k1-u4),MR4i-u4k1 ));
                            MR5i = u5k1 - max(0.0, min( min(u5k2-u5k1, u5k1-u5),MR5i-u5k1 ));
                            

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
                            
                            
                            
                            ML1i = u1 + max(0.0, min( min(u1k1-u1, u1-u1_k1), MR1-u1 ));
                            ML2i = u2 + max(0.0, min( min(u2k1-u2, u2-u2_k1), MR2-u2 ));
                            ML3i = u3 + max(0.0, min( min(u3k1-u3, u3-u3_k1), MR3-u3 ));
                            ML4i = u4 + max(0.0, min( min(u4k1-u4, u4-u4_k1), MR4-u4 ));
                            ML5i = u5 + max(0.0, min( min(u5k1-u5, u5-u5_k1), MR5-u5 ));
                            
                            MR1i = u1k1 - max(0.0, min( min(u1k2-u1k1, u1k1-u1),MR1i-u1k1 ));
                            MR2i = u2k1 - max(0.0, min( min(u2k2-u2k1, u2k1-u2),MR2i-u2k1 ));
                            MR3i = u3k1 - max(0.0, min( min(u3k2-u3k1, u3k1-u3),MR3i-u3k1 ));
                            MR4i = u4k1 - max(0.0, min( min(u4k2-u4k1, u4k1-u4),MR4i-u4k1 ));
                            MR5i = u5k1 - max(0.0, min( min(u5k2-u5k1, u5k1-u5),MR5i-u5k1 ));
                            

                        }
                        
                        
                        
                        
                        
                    }
                    else {

                        ML1 = u1_k1 + max( 0.0, min( u1_k1-u1_k2, u1-u1_k1 ) );
                        ML2 = u2_k1 + max( 0.0, min( u2_k1-u2_k2, u2-u2_k1 ) );
                        ML3 = u3_k1 + max( 0.0, min( u3_k1-u3_k2, u3-u3_k1 ) );
                        ML4 = u4_k1 + max( 0.0, min( u4_k1-u4_k2, u4-u4_k1 ) );
                        ML5 = u5_k1 + max( 0.0, min( u5_k1-u5_k2, u5-u5_k1 ) );
                        
                        MR1 = u1 - max( 0.0, min( u1k1-u1, u1-u1_k1 ) );
                        MR2 = u2 - max( 0.0, min( u2k1-u2, u2-u2_k1 ) );
                        MR3 = u3 - max( 0.0, min( u3k1-u3, u3-u3_k1 ) );
                        MR4 = u4 - max( 0.0, min( u4k1-u4, u4-u4_k1 ) );
                        MR5 = u5 - max( 0.0, min( u5k1-u5, u5-u5_k1 ) );
                        
                        ML1i = u1 + max( 0.0, min( u1-u1_k1, u1k1-u1 ) );
                        ML2i = u2 + max( 0.0, min( u2-u2_k1, u2k1-u2 ) );
                        ML3i = u3 + max( 0.0, min( u3-u3_k1, u3k1-u3 ) );
                        ML4i = u4 + max( 0.0, min( u4-u4_k1, u4k1-u4 ) );
                        ML5i = u5 + max( 0.0, min( u5-u5_k1, u5k1-u5 ) );
                        
                        MR1i = u1k1 - max( 0.0, min( u1k2-u1k1, u1k1-u1 ) );
                        MR2i = u2k1 - max( 0.0, min( u2k2-u2k1, u2k1-u2 ) );
                        MR3i = u3k1 - max( 0.0, min( u3k2-u3k1, u3k1-u3 ) );
                        MR4i = u4k1 - max( 0.0, min( u4k2-u4k1, u4k1-u4 ) );
                        MR5i = u5k1 - max( 0.0, min( u5k2-u5k1, u5k1-u5 ) );

                    }




                    #else


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

                    #endif

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
                    UN_rho = ML1;
                    UN_U = ML2/UN_rho;
                    UN_V = ML3/UN_rho;
                    UN_W = ML4/UN_rho;
                    UN_VV = UN_U*UN_U+UN_V*UN_V+UN_W*UN_W;
                    UN_P = (ML5-0.5*UN_rho*UN_VV)*(K-1);
                    UN_T = UN_P/UN_rho;
                    UN_C = K*UN_P/UN_rho;
                    UN_H = 0.5*UN_VV+UN_C/(K-1);

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

                    #if ROE != 4
                    
                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    #endif



                    #if ROE == 1 
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*W;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+W*W*(1-beta)*(1-beta));   // ---- C' ---- //

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = temp/S*(W_-UN_W);

                    temp3 = 0.5*(1-beta)*W*temp/S;

                    deltaU = (S-temp3-beta*fabs(W))*temp1+temp2;

                    deltaP = temp/S*(P_-UN_P)+(S-fabs(W)+temp3)*rho*(W_-UN_W);

                    #elif ROE == 2
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*W;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+W*W*(1-beta)*(1-beta));   // ---- C' ---- //


                    theda_p = VV/C;
                    U_p = 0.5*(1+theda_p)*W;
                    C_p = 0.5*sqrt(4*C*theda_p+W*W*(1-theda_p)*(1-theda_p));

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = U_p/S*(W_-UN_W);

                    temp3 = 0.5*(1-beta)*W*temp/S;

                    temp4 = 0.5*(1-theda_p)*W*U_p/S;

                    deltaU = (S-temp3-beta*fabs(W))*temp1+temp2;

                    deltaP = U_p/S*(P_-UN_P)+(C_p-fabs(W)+temp4)*rho*(W_-UN_W);

                    #elif ROE == 3

                    S = sqrt(C);

                    theda_p = (fabs(U)+fabs(V)+fabs(W))/S;

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);

                    beta = max(temp, 1e-8);
                    
                    C_p = beta*S;

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = W/S*beta*(W_-UN_W);

                    deltaU = (S-fabs(W))*temp1+temp2;

                    deltaP = W/S*(P_-UN_P)+(C_p-fabs(W))*rho*(W_-UN_W);

                    #elif ROE == 4

                    beta = sqrt(max(VV_/C_,UN_VV/UN_C));    // ---- theda ---- //

                    temp1 = 0.5*(UN_U+U_)+0.5*beta*(UN_U-U_);
                    temp2 = 0.5*(UN_V+V_)+0.5*beta*(UN_V-V_);
                    temp3 = 0.5*(UN_W+W_)+0.5*beta*(UN_W-W_);

                    U_ = 0.5*(U_+UN_U)+0.5*beta*(U_-UN_U);
                    V_ = 0.5*(V_+UN_V)+0.5*beta*(V_-UN_V);
                    W_ = 0.5*(W_+UN_W)+0.5*beta*(W_-UN_W);

                    UN_U = temp1;
                    UN_V = temp2;
                    UN_W = temp3;

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;


                    S = sqrt(C);

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = W/S*(W_-UN_W);

                    deltaU = (S-fabs(W))*temp1+temp2;

                    deltaP = W/S*(P_-UN_P)+(S-fabs(W))*rho*(W_-UN_W);

                    #elif ROE == 5
                    
                    S = sqrt(C);

                    theda_p = 0.5*(sqrt(VV_/C_)+sqrt(UN_VV/UN_C));

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);
                    
                    deltaP = W/S*(P_-UN_P)+temp*(S-fabs(W))*rho*(W_-UN_W);
                    
                    
                    // ==== deltaU_P + deltaU_RL ==== //            
                    temp1 = (P_-UN_P)/rho/C;

                    deltaU = ( 1.0-pow(temp,8) )*(S-fabs(W))*temp1;

                    temp2 = 0.5*W/S*(W_-UN_W);
                    
                    
                    // ==== KEASI ==== //
                    beta = sqrt(abs(W)/C);
                    
                    temp3 = beta*sqrt(4.0+(1.0-beta*beta)*(1.0-beta*beta))/(1.0+beta*beta);
                    
                    temp = max( 0.5*fabs(W_+UN_W)+(1.0-pow(temp3,8))*0.5*(W_-UN_W), 0.0 );

                    #endif






                    /* artificial viscosity */
                    #if ROE == 0
                      
                      S = sqrt(C);
                        
                      /* artificial viscosity */
                      Fav1 = 1.0/S*(P_-UN_P);
                      Fav2 =   U/S*(P_-UN_P);
                      Fav3 =   V/S*(P_-UN_P);
                      Fav4 =   W/S*(P_-UN_P);
                      Fav5 =   H/S*(P_-UN_P);
                      
                    #elif ROE == 5 
                      
                      Fav1 = Cdiss*temp*dU1+deltaU*rho  +temp2*(UN_rho+rho_);
                      Fav2 = Cdiss*temp*dU2+deltaU*rho*U+temp2*(UN_rho*UN_U+rho_*U_);
                      Fav3 = Cdiss*temp*dU3+deltaU*rho*V+temp2*(UN_rho*UN_V+rho_*V_);
                      Fav4 = Cdiss*temp*dU4+deltaU*rho*W+temp2*(UN_rho*UN_W+rho_*W_)+deltaP;
                      Fav5 = Cdiss*temp*dU5+deltaU*rho*H+temp2*(UN_rho*UN_H+rho_*H_)+deltaP*W;
                      
                    #else
                      
                      Fav1 = Cdiss*fabs(W)*dU1+deltaU*rho;
                      Fav2 = Cdiss*fabs(W)*dU2+deltaU*rho*U;
                      Fav3 = Cdiss*fabs(W)*dU3+deltaU*rho*V;
                      Fav4 = Cdiss*fabs(W)*dU4+deltaU*rho*W+deltaP;
                      Fav5 = Cdiss*fabs(W)*dU5+deltaU*rho*H+deltaP*W;

                    #endif

                    /* inviscid fluxes */

                    inFz1 = 0.5*((UN_rho*UN_W+rho_*W_-Ep*Fav1));
                    inFz2 = 0.5*((UN_rho*UN_U*UN_W+rho_*U_*W_)-Ep*Fav2);
                    inFz3 = 0.5*((UN_rho*UN_V*UN_W+rho_*V_*W_)-Ep*Fav3);
                    inFz4 = 0.5*((UN_rho*UN_W*UN_W+UN_P+rho_*W_*W_+P_)-Ep*Fav4);
                    inFz5 = 0.5*((UN_W*(3.5*UN_P+0.5*UN_rho*UN_VV)+W_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

                    dPz = 0.5*(UN_P+P_);
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
                    UN_rho = ML1i;
                    UN_U = ML2i/UN_rho;
                    UN_V = ML3i/UN_rho;
                    UN_W = ML4i/UN_rho;
                    UN_VV = UN_U*UN_U+UN_V*UN_V+UN_W*UN_W;
                    UN_P = (ML5i-0.5*UN_rho*UN_VV)*(K-1);
                    UN_T = UN_P/UN_rho;
                    UN_C = K*UN_P/UN_rho;
                    UN_H = 0.5*UN_VV+UN_C/(K-1);

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

                    
                    #if ROE != 4 

                    
                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    
                    #endif



                    #if ROE == 1

                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*W;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+W*W*(1-beta)*(1-beta));   // ---- C' ---- //


                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = temp/S*(W_-UN_W);

                    temp3 = 0.5*(1-beta)*W*temp/S;

                    deltaU = (S-temp3-beta*fabs(W))*temp1+temp2;

                    deltaP = temp/S*(P_-UN_P)+(S-fabs(W)+temp3)*rho*(W_-UN_W);

                    #elif ROE == 2

                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*W;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+W*W*(1-beta)*(1-beta));   // ---- C' ---- //


                    theda_p = VV/C;
                    U_p = 0.5*(1+theda_p)*W;
                    C_p = 0.5*sqrt(4*C*theda_p+W*W*(1-theda_p)*(1-theda_p));

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = U_p/S*(W_-UN_W);

                    temp3 = 0.5*(1-beta)*W*temp/S;

                    temp4 = 0.5*(1-theda_p)*W*U_p/S;

                    deltaU = (S-temp3-beta*fabs(W))*temp1+temp2;

                    deltaP = U_p/S*(P_-UN_P)+(C_p-fabs(W)+temp4)*rho*(W_-UN_W);
                    
                    #elif ROE == 3

                    S = sqrt(C);

                    theda_p = (fabs(U)+fabs(V)+fabs(W))/S;
                    
                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);

                    beta = max(temp, 1e-8);
                    
                    C_p = beta*S;

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = W/S*beta*(W_-UN_W);

                    deltaU = (S-fabs(W))*temp1+temp2;

                    deltaP = W/S*(P_-UN_P)+(C_p-fabs(W))*rho*(W_-UN_W);

                    #elif ROE == 4

                    beta = sqrt(max(VV_/C_,UN_VV/UN_C));    // ---- theda ---- //

                    temp1 = 0.5*(UN_U+U_)+0.5*beta*(UN_U-U_);
                    temp2 = 0.5*(UN_V+V_)+0.5*beta*(UN_V-V_);
                    temp3 = 0.5*(UN_W+W_)+0.5*beta*(UN_W-W_);

                    U_ = 0.5*(U_+UN_U)+0.5*beta*(U_-UN_U);
                    V_ = 0.5*(V_+UN_V)+0.5*beta*(V_-UN_V);
                    W_ = 0.5*(W_+UN_W)+0.5*beta*(W_-UN_W);

                    UN_U = temp1;
                    UN_V = temp2;
                    UN_W = temp3;

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;


                    S = sqrt(C);

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = W/S*(W_-UN_W);

                    deltaU = (S-fabs(W))*temp1+temp2;

                    deltaP = W/S*(P_-UN_P)+(S-fabs(W))*rho*(W_-UN_W);


                    #elif ROE == 5
                    
                    S = sqrt(C);

                    theda_p = 0.5*(sqrt(VV_/C_)+sqrt(UN_VV/UN_C));

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);
                    
                    deltaP = W/S*(P_-UN_P)+temp*(S-fabs(W))*rho*(W_-UN_W);
                    
                    
                    // ==== deltaU_P + deltaU_RL ==== //            
                    temp1 = (P_-UN_P)/rho/C;

                    deltaU = ( 1.0-pow(temp,8) )*(S-fabs(W))*temp1;

                    temp2 = 0.5*W/S*(W_-UN_W);
                    
                    
                    // ==== KEASI ==== //
                    beta = sqrt(abs(W)/C);
                    
                    temp3 = beta*sqrt(4.0+(1.0-beta*beta)*(1.0-beta*beta))/(1.0+beta*beta);
                    
                    temp = max( 0.5*fabs(W_+UN_W)+(1.0-pow(temp3,8))*0.5*(W_-UN_W), 0.0 );

                    #endif                                                                      
                    /* artificial viscosity */
                      
                    #if ROE == 0
                      
                      S = sqrt(C);
                        
                      /* artificial viscosity */
                      Fav1 = 1.0/S*(P_-UN_P);
                      Fav2 =   U/S*(P_-UN_P);
                      Fav3 =   V/S*(P_-UN_P);
                      Fav4 =   W/S*(P_-UN_P);
                      Fav5 =   H/S*(P_-UN_P);
                      
                    #elif ROE == 5 
                       
                      Fav1 = Cdiss*temp*dU1+deltaU*rho  +temp2*(UN_rho+rho_);
                      Fav2 = Cdiss*temp*dU2+deltaU*rho*U+temp2*(UN_rho*UN_U+rho_*U_);
                      Fav3 = Cdiss*temp*dU3+deltaU*rho*V+temp2*(UN_rho*UN_V+rho_*V_);
                      Fav4 = Cdiss*temp*dU4+deltaU*rho*W+temp2*(UN_rho*UN_W+rho_*W_)+deltaP;
                      Fav5 = Cdiss*temp*dU5+deltaU*rho*H+temp2*(UN_rho*UN_H+rho_*H_)+deltaP*W;
                      
                    #else
                      
                      Fav1 = Cdiss*fabs(W)*dU1+deltaU*rho;
                      Fav2 = Cdiss*fabs(W)*dU2+deltaU*rho*U;
                      Fav3 = Cdiss*fabs(W)*dU3+deltaU*rho*V;
                      Fav4 = Cdiss*fabs(W)*dU4+deltaU*rho*W+deltaP;
                      Fav5 = Cdiss*fabs(W)*dU5+deltaU*rho*H+deltaP*W;

                    #endif
                    /* inviscid fluxes */

                    inFz1i = 0.5*((UN_rho*UN_W+rho_*W_-Ep*Fav1));
                    inFz2i = 0.5*((UN_rho*UN_U*UN_W+rho_*U_*W_)-Ep*Fav2);
                    inFz3i = 0.5*((UN_rho*UN_V*UN_W+rho_*V_*W_)-Ep*Fav3);
                    inFz4i = 0.5*((UN_rho*UN_W*UN_W+UN_P+rho_*W_*W_+P_)-Ep*Fav4);
                    inFz5i = 0.5*((UN_W*(3.5*UN_P+0.5*UN_rho*UN_VV)+W_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);


                    dPzi = 0.5*(UN_P+P_);
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

                    #ifdef Kcomputer
                    //if(i > 2) {
                    #else
                    if(i > 2) {
                        #endif
                        

                        IF_i3 = FWS[icube][i-3][j][k];
                        u1_i3 = U1_[icube][i-3][j][k][0];
                        u2_i3 = U1_[icube][i-3][j][k][1];
                        u3_i3 = U1_[icube][i-3][j][k][2];
                        u4_i3 = U1_[icube][i-3][j][k][3];
                        u5_i3 = U1_[icube][i-3][j][k][4];

                        #ifdef Kcomputer
                        //}
                        #else
                    }
                    #endif

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

                    #ifdef Kcomputer
                    //if(i < nx)   {
                    #else
                    if(i < nx)   {
                        #endif

                        IFi3 = FWS[icube][i+3][j][k];
                        u1i3 = U1_[icube][i+3][j][k][0];
                        u2i3 = U1_[icube][i+3][j][k][1];
                        u3i3 = U1_[icube][i+3][j][k][2];
                        u4i3 = U1_[icube][i+3][j][k][3];
                        u5i3 = U1_[icube][i+3][j][k][4];

                        #ifdef Kcomputer
                        //}
                        #else
                    }
                    #endif


                    // -------------------------- X-direction --------------------------//
                    // -----------------------------------------------------------------//



                    // ------------------------------------------------------------------- //
                    // ----------------------------- MUSCL-X ----------------------------- //

                    #ifdef limiter
                    
                    
                    
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
                            
                            ML1 = u1_i1 + max(0.0, min( min(u1_i1-u1_i2, u1-u1_i1), ML1-u1_i1 ));
                            ML2 = u2_i1 + max(0.0, min( min(u2_i1-u2_i2, u2-u1_i1), ML2-u2_i1 ));
                            ML3 = u3_i1 + max(0.0, min( min(u3_i1-u3_i2, u3-u1_i1), ML3-u3_i1 ));
                            ML4 = u4_i1 + max(0.0, min( min(u4_i1-u4_i2, u4-u1_i1), ML4-u4_i1 ));
                            ML5 = u5_i1 + max(0.0, min( min(u5_i1-u5_i2, u5-u1_i1), ML5-u5_i1 ));
                            
                            MR1 = u1 - max(0.0, min( min(u1i1-u1, u1-u1_i1), MR1-u1 ));
                            MR2 = u2 - max(0.0, min( min(u2i1-u2, u2-u2_i1), MR2-u2 ));
                            MR3 = u3 - max(0.0, min( min(u3i1-u3, u3-u3_i1), MR3-u3 ));
                            MR4 = u4 - max(0.0, min( min(u4i1-u4, u4-u4_i1), MR4-u4 ));
                            MR5 = u5 - max(0.0, min( min(u5i1-u5, u5-u5_i1), MR5-u5 ));
                            
                            
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
                            
                            
                            ML1 = u1_i1 + max(0.0, min( min(u1_i1-u1_i2, u1-u1_i1), ML1-u1_i1 ));
                            ML2 = u2_i1 + max(0.0, min( min(u2_i1-u2_i2, u2-u1_i1), ML2-u2_i1 ));
                            ML3 = u3_i1 + max(0.0, min( min(u3_i1-u3_i2, u3-u1_i1), ML3-u3_i1 ));
                            ML4 = u4_i1 + max(0.0, min( min(u4_i1-u4_i2, u4-u1_i1), ML4-u4_i1 ));
                            ML5 = u5_i1 + max(0.0, min( min(u5_i1-u5_i2, u5-u1_i1), ML5-u5_i1 ));
                            
                            MR1 = u1 - max(0.0, min( min(u1i1-u1, u1-u1_i1), MR1-u1 ));
                            MR2 = u2 - max(0.0, min( min(u2i1-u2, u2-u2_i1), MR2-u2 ));
                            MR3 = u3 - max(0.0, min( min(u3i1-u3, u3-u3_i1), MR3-u3 ));
                            MR4 = u4 - max(0.0, min( min(u4i1-u4, u4-u4_i1), MR4-u4 ));
                            MR5 = u5 - max(0.0, min( min(u5i1-u5, u5-u5_i1), MR5-u5 ));

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
                            
                            
                            ML1i = u1 + max(0.0, min( min(u1i1-u1, u1-u1_i1), MR1-u1 ));
                            ML2i = u2 + max(0.0, min( min(u2i1-u2, u2-u2_i1), MR2-u2 ));
                            ML3i = u3 + max(0.0, min( min(u3i1-u3, u3-u3_i1), MR3-u3 ));
                            ML4i = u4 + max(0.0, min( min(u4i1-u4, u4-u4_i1), MR4-u4 ));
                            ML5i = u5 + max(0.0, min( min(u5i1-u5, u5-u5_i1), MR5-u5 ));
                            
                            MR1i = u1i1 - max(0.0, min( min(u1i2-u1i1, u1i1-u1),MR1i-u1i1 ));
                            MR2i = u2i1 - max(0.0, min( min(u2i2-u2i1, u2i1-u2),MR2i-u2i1 ));
                            MR3i = u3i1 - max(0.0, min( min(u3i2-u3i1, u3i1-u3),MR3i-u3i1 ));
                            MR4i = u4i1 - max(0.0, min( min(u4i2-u4i1, u4i1-u4),MR4i-u4i1 ));
                            MR5i = u5i1 - max(0.0, min( min(u5i2-u5i1, u5i1-u5),MR5i-u5i1 ));
                            

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
                            
                            
                            
                            ML1i = u1 + max(0.0, min( min(u1i1-u1, u1-u1_i1), MR1-u1 ));
                            ML2i = u2 + max(0.0, min( min(u2i1-u2, u2-u2_i1), MR2-u2 ));
                            ML3i = u3 + max(0.0, min( min(u3i1-u3, u3-u3_i1), MR3-u3 ));
                            ML4i = u4 + max(0.0, min( min(u4i1-u4, u4-u4_i1), MR4-u4 ));
                            ML5i = u5 + max(0.0, min( min(u5i1-u5, u5-u5_i1), MR5-u5 ));
                            
                            MR1i = u1i1 - max(0.0, min( min(u1i2-u1i1, u1i1-u1),MR1i-u1i1 ));
                            MR2i = u2i1 - max(0.0, min( min(u2i2-u2i1, u2i1-u2),MR2i-u2i1 ));
                            MR3i = u3i1 - max(0.0, min( min(u3i2-u3i1, u3i1-u3),MR3i-u3i1 ));
                            MR4i = u4i1 - max(0.0, min( min(u4i2-u4i1, u4i1-u4),MR4i-u4i1 ));
                            MR5i = u5i1 - max(0.0, min( min(u5i2-u5i1, u5i1-u5),MR5i-u5i1 ));
                            

                        }
                        
                        
                        
                        
                        
                    }
                    else {

                        ML1 = u1_i1 + max( 0.0, min( u1_i1-u1_i2, u1-u1_i1 ) );
                        ML2 = u2_i1 + max( 0.0, min( u2_i1-u2_i2, u2-u2_i1 ) );
                        ML3 = u3_i1 + max( 0.0, min( u3_i1-u3_i2, u3-u3_i1 ) );
                        ML4 = u4_i1 + max( 0.0, min( u4_i1-u4_i2, u4-u4_i1 ) );
                        ML5 = u5_i1 + max( 0.0, min( u5_i1-u5_i2, u5-u5_i1 ) );
                        
                        MR1 = u1 - max( 0.0, min( u1i1-u1, u1-u1_i1 ) );
                        MR2 = u2 - max( 0.0, min( u2i1-u2, u2-u2_i1 ) );
                        MR3 = u3 - max( 0.0, min( u3i1-u3, u3-u3_i1 ) );
                        MR4 = u4 - max( 0.0, min( u4i1-u4, u4-u4_i1 ) );
                        MR5 = u5 - max( 0.0, min( u5i1-u5, u5-u5_i1 ) );
                        
                        ML1i = u1 + max( 0.0, min( u1-u1_i1, u1i1-u1 ) );
                        ML2i = u2 + max( 0.0, min( u2-u2_i1, u2i1-u2 ) );
                        ML3i = u3 + max( 0.0, min( u3-u3_i1, u3i1-u3 ) );
                        ML4i = u4 + max( 0.0, min( u4-u4_i1, u4i1-u4 ) );
                        ML5i = u5 + max( 0.0, min( u5-u5_i1, u5i1-u5 ) );
                        
                        MR1i = u1i1 - max( 0.0, min( u1i2-u1i1, u1i1-u1 ) );
                        MR2i = u2i1 - max( 0.0, min( u2i2-u2i1, u2i1-u2 ) );
                        MR3i = u3i1 - max( 0.0, min( u3i2-u3i1, u3i1-u3 ) );
                        MR4i = u4i1 - max( 0.0, min( u4i2-u4i1, u4i1-u4 ) );
                        MR5i = u5i1 - max( 0.0, min( u5i2-u5i1, u5i1-u5 ) );

                    }
                    
                    

                    #else



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

                    #endif
                    
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
                    UN_rho = ML1;
                    UN_U = ML2/UN_rho;
                    UN_V = ML3/UN_rho;
                    UN_W = ML4/UN_rho;
                    UN_VV = UN_U*UN_U+UN_V*UN_V+UN_W*UN_W;
                    UN_P = (ML5-0.5*UN_rho*UN_VV)*(K-1);
                    UN_T = UN_P/UN_rho;
                    UN_C = K*UN_P/UN_rho;
                    UN_H = 0.5*UN_VV+UN_C/(K-1);

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

                    #if ROE != 4 

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    #endif




                    #if ROE == 1
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*U;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+U*U*(1-beta)*(1-beta));   // ---- C' ---- //

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = temp/S*(U_-UN_U);

                    temp3 = 0.5*(1-beta)*U*temp/S;



                    deltaU = (S-temp3-beta*fabs(U))*temp1+temp2;

                    deltaP = temp/S*(P_-UN_P)+(S-fabs(U)+temp3)*rho*(U_-UN_U);

                    
                    #elif ROE == 2
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*U;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+U*U*(1-beta)*(1-beta));   // ---- C' ---- //

                    theda_p = VV/C;

                    U_p = 0.5*(1+theda_p)*U;
                    C_p = 0.5*sqrt(4*C*theda_p+U*U*(1-theda_p)*(1-theda_p));

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = U_p/S*(U_-UN_U);

                    temp3 = 0.5*(1-beta)*U*temp/S;


                    
                    temp4 = 0.5*(1-theda_p)*U*U_p/S;

                    deltaU = (S-temp3-beta*fabs(U))*temp1+temp2;

                    deltaP = U_p/S*(P_-UN_P)+(C_p-fabs(U)+temp4)*rho*(U_-UN_U);

                    
                    #elif ROE == 3

                    S = sqrt(C);

                    theda_p = (fabs(U)+fabs(V)+fabs(W))/S;
                    
                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);

                    beta = max(temp, 1e-8);
                    
                    C_p = beta*S;

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = U/S*beta*(U_-UN_U);

                    deltaU = (S-fabs(U))*temp1+temp2;

                    deltaP = U/S*(P_-UN_P)+(C_p-fabs(U))*rho*(U_-UN_U);

                    #elif ROE == 4

                    
                    beta = sqrt(max(VV_/C_,UN_VV/UN_C));    // ---- theda ---- //

                    temp1 = 0.5*(UN_U+U_)+0.5*beta*(UN_U-U_);
                    temp2 = 0.5*(UN_V+V_)+0.5*beta*(UN_V-V_);
                    temp3 = 0.5*(UN_W+W_)+0.5*beta*(UN_W-W_);

                    U_ = 0.5*(U_+UN_U)+0.5*beta*(U_-UN_U);
                    V_ = 0.5*(V_+UN_V)+0.5*beta*(V_-UN_V);
                    W_ = 0.5*(W_+UN_W)+0.5*beta*(W_-UN_W);

                    UN_U = temp1;
                    UN_V = temp2;
                    UN_W = temp3;

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;


                    S = sqrt(C);

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = U/S*(U_-UN_U);

                    deltaU = (S-fabs(U))*temp1+temp2;

                    deltaP = U/S*(P_-UN_P)+(S-fabs(U))*rho*(U_-UN_U);

                    #elif ROE == 5
                    
                    S = sqrt(C);

                    theda_p = 0.5*(sqrt(VV_/C_)+sqrt(UN_VV/UN_C));

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);
                    
                    deltaP = U/S*(P_-UN_P)+temp*(S-fabs(U))*rho*(U_-UN_U);
                    
                    
                    // ==== deltaU_P + deltaU_RL ==== //            
                    temp1 = (P_-UN_P)/rho/C;

                    deltaU = ( 1.0-pow(temp,8) )*(S-fabs(U))*temp1;

                    temp2 = 0.5*U/S*(U_-UN_U);
                    
                    
                    // ==== KEASI ==== //
                    beta = sqrt(abs(U)/C);
                    
                    temp3 = beta*sqrt(4.0+(1.0-beta*beta)*(1.0-beta*beta))/(1.0+beta*beta);
                    
                    temp = max( 0.5*fabs(U_+UN_U)+(1.0-pow(temp3,8))*0.5*(U_-UN_U), 0.0 );

                    #endif
                    /* artificial viscosity */
                    
                   #if ROE == 0
                      
                      S = sqrt(C);
                        
                      /* artificial viscosity */
                      Fav1 = 1.0/S*(P_-UN_P);
                      Fav2 =   U/S*(P_-UN_P);
                      Fav3 =   V/S*(P_-UN_P);
                      Fav4 =   W/S*(P_-UN_P);
                      Fav5 =   H/S*(P_-UN_P);
                      
                    #elif ROE == 5 
                      
                    
                    Fav1 = Cdiss*temp*dU1+deltaU*rho  +temp2*(UN_rho+rho_);
                    Fav2 = Cdiss*temp*dU2+deltaU*rho*U+temp2*(UN_rho*UN_U+rho_*U_)+deltaP;
                    Fav3 = Cdiss*temp*dU3+deltaU*rho*V+temp2*(UN_rho*UN_V+rho_*V_);
                    Fav4 = Cdiss*temp*dU4+deltaU*rho*W+temp2*(UN_rho*UN_W+rho_*W_);
                    Fav5 = Cdiss*temp*dU5+deltaU*rho*H+temp2*(UN_rho*UN_H+rho_*H_)+deltaP*U;
                    
                    #else
                    
                    Fav1 = Cdiss*fabs(U)*dU1+deltaU*rho;
                    Fav2 = Cdiss*fabs(U)*dU2+deltaU*rho*U+deltaP;
                    Fav3 = Cdiss*fabs(U)*dU3+deltaU*rho*V;
                    Fav4 = Cdiss*fabs(U)*dU4+deltaU*rho*W;
                    Fav5 = Cdiss*fabs(U)*dU5+deltaU*rho*H+deltaP*U;
                    

                    #endif
                    /* inviscid fluxes */
                    inFx1 = 0.5*((UN_rho*UN_U+rho_*U_)-Ep*Fav1);
                    inFx2 = 0.5*((UN_rho*UN_U*UN_U+UN_P+rho_*U_*U_+P_)-Ep*Fav2);
                    inFx3 = 0.5*((UN_rho*UN_V*UN_U+rho_*V_*U_)-Ep*Fav3);
                    inFx4 = 0.5*((UN_rho*UN_W*UN_U+rho_*W_*U_)-Ep*Fav4);
                    inFx5 = 0.5*((UN_U*(3.5*UN_P+0.5*UN_rho*UN_VV)+U_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

                    
                    dPx = 0.5*(UN_P+P_);
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
                    UN_rho = ML1i;
                    UN_U = ML2i/UN_rho;
                    UN_V = ML3i/UN_rho;
                    UN_W = ML4i/UN_rho;
                    UN_VV = UN_U*UN_U+UN_V*UN_V+UN_W*UN_W;
                    UN_P = (ML5i-0.5*UN_rho*UN_VV)*(K-1);
                    UN_T = UN_P/UN_rho;
                    UN_C = K*UN_P/UN_rho;
                    UN_H = 0.5*UN_VV+UN_C/(K-1);

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
                    
                    #if ROE != 4

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    #endif




                    #if ROE == 1
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    temp = 0.5*(1+beta)*U;    // ---- U' ---- //
                    S = 0.5*sqrt(4*beta*C+U*U*(1-beta)*(1-beta));   // ---- C' ---- //

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = temp/S*(U_-UN_U);

                    temp3 = 0.5*(1-beta)*U*temp/S;

                    deltaU = (S-temp3-beta*fabs(U))*temp1+temp2;

                    deltaP = temp/S*(P_-UN_P)+(S-fabs(U)+temp3)*rho*(U_-UN_U);

                    #elif ROE == 2
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    
                    temp = 0.5*(1+beta)*U;    // ---- U' ---- //

                    S = 0.5*sqrt(4*beta*C+U*U*(1-beta)*(1-beta));   // ---- C' ---- //


                    theda_p = VV/C;



                    U_p = 0.5*(1+theda_p)*U;
                    C_p = 0.5*sqrt(4*C*theda_p+U*U*(1-theda_p)*(1-theda_p));

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = U_p/S*(U_-UN_U);

                    temp3 = 0.5*(1-beta)*U*temp/S;


                    
                    temp4 = 0.5*(1-theda_p)*U*U_p/S;

                    deltaU = (S-temp3-beta*fabs(U))*temp1+temp2;

                    deltaP = U_p/S*(P_-UN_P)+(C_p-fabs(U)+temp4)*rho*(U_-UN_U);

                    #elif ROE == 3

                    S = sqrt(C);

                    theda_p = (fabs(U)+fabs(V)+fabs(W))/S;
                    
                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);

                    beta = max(temp, 1e-8);
                    
                    C_p = beta*S;

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = U/S*beta*(U_-UN_U);

                    deltaU = (S-fabs(U))*temp1+temp2;

                    deltaP = U/S*(P_-UN_P)+(C_p-fabs(U))*rho*(U_-UN_U);
                    
                    #elif ROE == 4

                    
                    beta = sqrt(max(VV_/C_,UN_VV/UN_C));    // ---- theda ---- //

                    temp1 = 0.5*(UN_U+U_)+0.5*beta*(UN_U-U_);
                    temp2 = 0.5*(UN_V+V_)+0.5*beta*(UN_V-V_);
                    temp3 = 0.5*(UN_W+W_)+0.5*beta*(UN_W-W_);

                    U_ = 0.5*(U_+UN_U)+0.5*beta*(U_-UN_U);
                    V_ = 0.5*(V_+UN_V)+0.5*beta*(V_-UN_V);
                    W_ = 0.5*(W_+UN_W)+0.5*beta*(W_-UN_W);

                    UN_U = temp1;
                    UN_V = temp2;
                    UN_W = temp3;

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    S = sqrt(C);

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = U/S*(U_-UN_U);

                    deltaU = (S-fabs(U))*temp1+temp2;

                    deltaP = U/S*(P_-UN_P)+(S-fabs(U))*rho*(U_-UN_U);
                    
                    
                    #elif ROE == 5
                    
                    S = sqrt(C);

                    theda_p = 0.5*(sqrt(VV_/C_)+sqrt(UN_VV/UN_C));

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);
                    
                    deltaP = U/S*(P_-UN_P)+temp*(S-fabs(U))*rho*(U_-UN_U);
                    
                    
                    // ==== deltaU_P + deltaU_RL ==== //            
                    temp1 = (P_-UN_P)/rho/C;

                    deltaU = ( 1.0-pow(temp,8) )*(S-fabs(U))*temp1;

                    temp2 = 0.5*U/S*(U_-UN_U);
                    
                    
                    // ==== KEASI ==== //
                    beta = sqrt(abs(U)/C);
                    
                    temp3 = beta*sqrt(4.0+(1.0-beta*beta)*(1.0-beta*beta))/(1.0+beta*beta);
                    
                    temp = max( 0.5*fabs(U_+UN_U)+(1.0-pow(temp3,8))*0.5*(U_-UN_U), 0.0 );

                    #endif

                    /* artificial viscosity */
                    #if ROE == 0
                      
                      S = sqrt(C);
                        
                      /* artificial viscosity */
                      Fav1 = 1.0/S*(P_-UN_P);
                      Fav2 =   U/S*(P_-UN_P);
                      Fav3 =   V/S*(P_-UN_P);
                      Fav4 =   W/S*(P_-UN_P);
                      Fav5 =   H/S*(P_-UN_P);
                      
                    #elif ROE == 5 
                      
                    
                    Fav1 = Cdiss*temp*dU1+deltaU*rho  +temp2*(UN_rho+rho_);
                    Fav2 = Cdiss*temp*dU2+deltaU*rho*U+temp2*(UN_rho*UN_U+rho_*U_)+deltaP;
                    Fav3 = Cdiss*temp*dU3+deltaU*rho*V+temp2*(UN_rho*UN_V+rho_*V_);
                    Fav4 = Cdiss*temp*dU4+deltaU*rho*W+temp2*(UN_rho*UN_W+rho_*W_);
                    Fav5 = Cdiss*temp*dU5+deltaU*rho*H+temp2*(UN_rho*UN_H+rho_*H_)+deltaP*U;
                    
                    #else
                    
                    Fav1 = Cdiss*fabs(U)*dU1+deltaU*rho;
                    Fav2 = Cdiss*fabs(U)*dU2+deltaU*rho*U+deltaP;
                    Fav3 = Cdiss*fabs(U)*dU3+deltaU*rho*V;
                    Fav4 = Cdiss*fabs(U)*dU4+deltaU*rho*W;
                    Fav5 = Cdiss*fabs(U)*dU5+deltaU*rho*H+deltaP*U;
                    

                    #endif

                    /* inviscid fluxes */
                    inFx1i = 0.5*((UN_rho*UN_U+rho_*U_)-Ep*Fav1);
                    inFx2i = 0.5*((UN_rho*UN_U*UN_U+UN_P+rho_*U_*U_+P_)-Ep*Fav2);
                    inFx3i = 0.5*((UN_rho*UN_V*UN_U+rho_*V_*U_)-Ep*Fav3);
                    inFx4i = 0.5*((UN_rho*UN_W*UN_U+rho_*W_*U_)-Ep*Fav4);
                    inFx5i = 0.5*((UN_U*(3.5*UN_P+0.5*UN_rho*UN_VV)+U_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

                    
                    dPxi = 0.5*(UN_P+P_);
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

                    #ifdef Kcomputer
                    //if(j > 2) {
                    #else
                    if(j > 2) {
                        #endif
                        

                        IF_j3 = FWS[icube][i][j-3][k];
                        u1_j3 = U1_[icube][i][j-3][k][0];
                        u2_j3 = U1_[icube][i][j-3][k][1];
                        u3_j3 = U1_[icube][i][j-3][k][2];
                        u4_j3 = U1_[icube][i][j-3][k][3];
                        u5_j3 = U1_[icube][i][j-3][k][4];

                        #ifdef Kcomputer
                        //}
                        #else
                    }
                    #endif

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

                    #ifdef Kcomputer
                    //if(j < ny)   {
                    #else
                    if(j < ny)   {
                        #endif
                        

                        IFj3 = FWS[icube][i][j+3][k];
                        u1j3 = U1_[icube][i][j+3][k][0];
                        u2j3 = U1_[icube][i][j+3][k][1];
                        u3j3 = U1_[icube][i][j+3][k][2];
                        u4j3 = U1_[icube][i][j+3][k][3];
                        u5j3 = U1_[icube][i][j+3][k][4];

                        #ifdef Kcomputer
                        //}
                        #else
                    }
                    #endif
                    

                    // -------------------------- Y-direction --------------------------//
                    // -----------------------------------------------------------------//


                    // ------------------------------------------------------------------- //
                    // ----------------------------- MUSCL-Y ----------------------------- //

                    #ifdef limiter
                    
                    
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
                            
                            ML1 = u1_j1 + max(0.0, min( min(u1_j1-u1_j2, u1-u1_j1), ML1-u1_j1 ));
                            ML2 = u2_j1 + max(0.0, min( min(u2_j1-u2_j2, u2-u1_j1), ML2-u2_j1 ));
                            ML3 = u3_j1 + max(0.0, min( min(u3_j1-u3_j2, u3-u1_j1), ML3-u3_j1 ));
                            ML4 = u4_j1 + max(0.0, min( min(u4_j1-u4_j2, u4-u1_j1), ML4-u4_j1 ));
                            ML5 = u5_j1 + max(0.0, min( min(u5_j1-u5_j2, u5-u1_j1), ML5-u5_j1 ));
                            
                            MR1 = u1 - max(0.0, min( min(u1j1-u1, u1-u1_j1), MR1-u1 ));
                            MR2 = u2 - max(0.0, min( min(u2j1-u2, u2-u2_j1), MR2-u2 ));
                            MR3 = u3 - max(0.0, min( min(u3j1-u3, u3-u3_j1), MR3-u3 ));
                            MR4 = u4 - max(0.0, min( min(u4j1-u4, u4-u4_j1), MR4-u4 ));
                            MR5 = u5 - max(0.0, min( min(u5j1-u5, u5-u5_j1), MR5-u5 ));
                            
                            
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
                            
                            
                            ML1 = u1_j1 + max(0.0, min( min(u1_j1-u1_j2, u1-u1_j1), ML1-u1_j1 ));
                            ML2 = u2_j1 + max(0.0, min( min(u2_j1-u2_j2, u2-u1_j1), ML2-u2_j1 ));
                            ML3 = u3_j1 + max(0.0, min( min(u3_j1-u3_j2, u3-u1_j1), ML3-u3_j1 ));
                            ML4 = u4_j1 + max(0.0, min( min(u4_j1-u4_j2, u4-u1_j1), ML4-u4_j1 ));
                            ML5 = u5_j1 + max(0.0, min( min(u5_j1-u5_j2, u5-u1_j1), ML5-u5_j1 ));
                            
                            MR1 = u1 - max(0.0, min( min(u1j1-u1, u1-u1_j1), MR1-u1 ));
                            MR2 = u2 - max(0.0, min( min(u2j1-u2, u2-u2_j1), MR2-u2 ));
                            MR3 = u3 - max(0.0, min( min(u3j1-u3, u3-u3_j1), MR3-u3 ));
                            MR4 = u4 - max(0.0, min( min(u4j1-u4, u4-u4_j1), MR4-u4 ));
                            MR5 = u5 - max(0.0, min( min(u5j1-u5, u5-u5_j1), MR5-u5 ));

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
                            
                            
                            ML1i = u1 + max(0.0, min( min(u1j1-u1, u1-u1_j1), MR1-u1 ));
                            ML2i = u2 + max(0.0, min( min(u2j1-u2, u2-u2_j1), MR2-u2 ));
                            ML3i = u3 + max(0.0, min( min(u3j1-u3, u3-u3_j1), MR3-u3 ));
                            ML4i = u4 + max(0.0, min( min(u4j1-u4, u4-u4_j1), MR4-u4 ));
                            ML5i = u5 + max(0.0, min( min(u5j1-u5, u5-u5_j1), MR5-u5 ));
                            
                            MR1i = u1j1 - max(0.0, min( min(u1j2-u1j1, u1j1-u1),MR1i-u1j1 ));
                            MR2i = u2j1 - max(0.0, min( min(u2j2-u2j1, u2j1-u2),MR2i-u2j1 ));
                            MR3i = u3j1 - max(0.0, min( min(u3j2-u3j1, u3j1-u3),MR3i-u3j1 ));
                            MR4i = u4j1 - max(0.0, min( min(u4j2-u4j1, u4j1-u4),MR4i-u4j1 ));
                            MR5i = u5j1 - max(0.0, min( min(u5j2-u5j1, u5j1-u5),MR5i-u5j1 ));
                            

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
                            
                            
                            
                            ML1i = u1 + max(0.0, min( min(u1j1-u1, u1-u1_j1), MR1-u1 ));
                            ML2i = u2 + max(0.0, min( min(u2j1-u2, u2-u2_j1), MR2-u2 ));
                            ML3i = u3 + max(0.0, min( min(u3j1-u3, u3-u3_j1), MR3-u3 ));
                            ML4i = u4 + max(0.0, min( min(u4j1-u4, u4-u4_j1), MR4-u4 ));
                            ML5i = u5 + max(0.0, min( min(u5j1-u5, u5-u5_j1), MR5-u5 ));
                            
                            MR1i = u1j1 - max(0.0, min( min(u1j2-u1j1, u1j1-u1),MR1i-u1j1 ));
                            MR2i = u2j1 - max(0.0, min( min(u2j2-u2j1, u2j1-u2),MR2i-u2j1 ));
                            MR3i = u3j1 - max(0.0, min( min(u3j2-u3j1, u3j1-u3),MR3i-u3j1 ));
                            MR4i = u4j1 - max(0.0, min( min(u4j2-u4j1, u4j1-u4),MR4i-u4j1 ));
                            MR5i = u5j1 - max(0.0, min( min(u5j2-u5j1, u5j1-u5),MR5i-u5j1 ));
                            

                        }
                        
                        
                        
                        
                        
                    }
                    else {

                        ML1 = u1_j1 + max( 0.0, min( u1_j1-u1_j2, u1-u1_j1 ) );
                        ML2 = u2_j1 + max( 0.0, min( u2_j1-u2_j2, u2-u2_j1 ) );
                        ML3 = u3_j1 + max( 0.0, min( u3_j1-u3_j2, u3-u3_j1 ) );
                        ML4 = u4_j1 + max( 0.0, min( u4_j1-u4_j2, u4-u4_j1 ) );
                        ML5 = u5_j1 + max( 0.0, min( u5_j1-u5_j2, u5-u5_j1 ) );
                        
                        MR1 = u1 - max( 0.0, min( u1j1-u1, u1-u1_j1 ) );
                        MR2 = u2 - max( 0.0, min( u2j1-u2, u2-u2_j1 ) );
                        MR3 = u3 - max( 0.0, min( u3j1-u3, u3-u3_j1 ) );
                        MR4 = u4 - max( 0.0, min( u4j1-u4, u4-u4_j1 ) );
                        MR5 = u5 - max( 0.0, min( u5j1-u5, u5-u5_j1 ) );
                        
                        ML1i = u1 + max( 0.0, min( u1-u1_j1, u1j1-u1 ) );
                        ML2i = u2 + max( 0.0, min( u2-u2_j1, u2j1-u2 ) );
                        ML3i = u3 + max( 0.0, min( u3-u3_j1, u3j1-u3 ) );
                        ML4i = u4 + max( 0.0, min( u4-u4_j1, u4j1-u4 ) );
                        ML5i = u5 + max( 0.0, min( u5-u5_j1, u5j1-u5 ) );
                        
                        MR1i = u1j1 - max( 0.0, min( u1j2-u1j1, u1j1-u1 ) );
                        MR2i = u2j1 - max( 0.0, min( u2j2-u2j1, u2j1-u2 ) );
                        MR3i = u3j1 - max( 0.0, min( u3j2-u3j1, u3j1-u3 ) );
                        MR4i = u4j1 - max( 0.0, min( u4j2-u4j1, u4j1-u4 ) );
                        MR5i = u5j1 - max( 0.0, min( u5j2-u5j1, u5j1-u5 ) );

                    }
                    
                    
                    
                    #else





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
                    
                    #endif


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
                    UN_rho = ML1;
                    UN_U = ML2/UN_rho;
                    UN_V = ML3/UN_rho;
                    UN_W = ML4/UN_rho;
                    UN_VV = UN_U*UN_U+UN_V*UN_V+UN_W*UN_W;
                    UN_P = (ML5-0.5*UN_rho*UN_VV)*(K-1);
                    UN_T = UN_P/UN_rho;
                    UN_C = K*UN_P/UN_rho;
                    UN_H = 0.5*UN_VV+UN_C/(K-1);

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

                    #if ROE != 4

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    #endif


                    #if ROE == 1 

                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    temp = 0.5*(1+beta)*V;    // ---- U' ---- //
                    S = 0.5*sqrt(4*beta*C+V*V*(1-beta)*(1-beta));   // ---- C' ---- //

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = temp/S*(V_-UN_V);

                    temp3 = 0.5*(1-beta)*V*temp/S;

                    deltaU = (S-temp3-beta*fabs(V))*temp1+temp2;

                    deltaP = temp/S*(P_-UN_P)+(S-fabs(V)+temp3)*rho*(V_-UN_V);

                    #elif ROE == 2
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    temp = 0.5*(1+beta)*V;    // ---- U' ---- //
                    S = 0.5*sqrt(4*beta*C+V*V*(1-beta)*(1-beta));   // ---- C' ---- //


                    theda_p = VV/C;
                    U_p = 0.5*(1+theda_p)*V;
                    C_p = 0.5*sqrt(4*C*theda_p+V*V*(1-theda_p)*(1-theda_p));

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = U_p/S*(V_-UN_V);

                    temp3 = 0.5*(1-beta)*V*temp/S;


                    
                    temp4 = 0.5*(1-theda_p)*V*U_p/S;

                    deltaU = (S-temp3-beta*fabs(V))*temp1+temp2;


                    deltaP = U_p/S*(P_-UN_P)+(C_p-fabs(V)+temp4)*rho*(V_-UN_V);

                    #elif ROE == 3

                    S = sqrt(C);

                    theda_p = (fabs(U)+fabs(V)+fabs(W))/S;
                    
                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);

                    beta = max(temp, 1e-8);
                    
                    C_p = beta*S;

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = V/S*beta*(V_-UN_V);

                    deltaU = (S-fabs(V))*temp1+temp2;

                    deltaP = V/S*(P_-UN_P)+(C_p-fabs(V))*rho*(V_-UN_V);

                    #elif ROE == 4
                    
                    
                    beta = sqrt(max(VV_/C_,UN_VV/UN_C));    // ---- theda ---- //

                    temp1 = 0.5*(UN_U+U_)+0.5*beta*(UN_U-U_);
                    temp2 = 0.5*(UN_V+V_)+0.5*beta*(UN_V-V_);
                    temp3 = 0.5*(UN_W+W_)+0.5*beta*(UN_W-W_);

                    U_ = 0.5*(U_+UN_U)+0.5*beta*(U_-UN_U);
                    V_ = 0.5*(V_+UN_V)+0.5*beta*(V_-UN_V);
                    W_ = 0.5*(W_+UN_W)+0.5*beta*(W_-UN_W);

                    UN_U = temp1;
                    UN_V = temp2;
                    UN_W = temp3;

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    S = sqrt(C);

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = V/S*(V_-UN_V);

                    deltaU = (S-fabs(V))*temp1+temp2;

                    deltaP = V/S*(P_-UN_P)+(S-fabs(V))*rho*(V_-UN_V);

                    #elif ROE == 5
                    
                    S = sqrt(C);

                    theda_p = 0.5*(sqrt(VV_/C_)+sqrt(UN_VV/UN_C));

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);
                    
                    deltaP = V/S*(P_-UN_P)+temp*(S-fabs(V))*rho*(V_-UN_V);
                    
                    
                    // ==== deltaU_P + deltaU_RL ==== //            
                    temp1 = (P_-UN_P)/rho/C;

                    deltaU = ( 1.0-pow(temp,8) )*(S-fabs(V))*temp1;

                    temp2 = 0.5*V/S*(V_-UN_V);
                    
                    
                    // ==== KEASI ==== //
                    beta = sqrt(abs(V)/C);
                    
                    temp3 = beta*sqrt(4.0+(1.0-beta*beta)*(1.0-beta*beta))/(1.0+beta*beta);
                    
                    temp = max( 0.5*fabs(V_+UN_V)+(1.0-pow(temp3,8))*0.5*(V_-UN_V), 0.0 );


                    #endif

                    /* artificial viscosity */
                    
                    #if ROE == 0
                      
                      S = sqrt(C);
                        
                      /* artificial viscosity */
                      Fav1 = 1.0/S*(P_-UN_P);
                      Fav2 =   U/S*(P_-UN_P);
                      Fav3 =   V/S*(P_-UN_P);
                      Fav4 =   W/S*(P_-UN_P);
                      Fav5 =   H/S*(P_-UN_P);
                      
                    #elif ROE == 5 
                      
                    
                    Fav1 = Cdiss*temp*dU1+deltaU*rho  +temp2*(UN_rho+rho_);
                    Fav2 = Cdiss*temp*dU2+deltaU*rho*U+temp2*(UN_rho*UN_U+rho_*U_);
                    Fav3 = Cdiss*temp*dU3+deltaU*rho*V+temp2*(UN_rho*UN_V+rho_*V_)+deltaP;
                    Fav4 = Cdiss*temp*dU4+deltaU*rho*W+temp2*(UN_rho*UN_W+rho_*W_);
                    Fav5 = Cdiss*temp*dU5+deltaU*rho*H+temp2*(UN_rho*UN_H+rho_*H_)+deltaP*V;
                    
                    #else
                    
                    Fav1 = Cdiss*fabs(V)*dU1+deltaU*rho;
                    Fav2 = Cdiss*fabs(V)*dU2+deltaU*rho*U;
                    Fav3 = Cdiss*fabs(V)*dU3+deltaU*rho*V+deltaP;
                    Fav4 = Cdiss*fabs(V)*dU4+deltaU*rho*W;
                    Fav5 = Cdiss*fabs(V)*dU5+deltaU*rho*H+deltaP*V;
                    

                    #endif
                    /* inviscid fluxes */
                    inFy1 = 0.5*((UN_rho*UN_V+rho_*V_-Ep*Fav1));
                    inFy2 = 0.5*((UN_rho*UN_U*UN_V+rho_*U_*V_)-Ep*Fav2);
                    inFy3 = 0.5*((UN_rho*UN_V*UN_V+UN_P+rho_*V_*V_+P_)-Ep*Fav3);
                    inFy4 = 0.5*((UN_rho*UN_W*UN_V+rho_*W_*V_)-Ep*Fav4);
                    inFy5 = 0.5*((UN_V*(3.5*UN_P+0.5*UN_rho*UN_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);

                    dPy = 0.5*(UN_P+P_);


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
                    UN_rho = ML1i;
                    UN_U = ML2i/UN_rho;
                    UN_V = ML3i/UN_rho;
                    UN_W = ML4i/UN_rho;
                    UN_VV = UN_U*UN_U+UN_V*UN_V+UN_W*UN_W;
                    UN_P = (ML5i-0.5*UN_rho*UN_VV)*(K-1);
                    UN_T = UN_P/UN_rho;
                    UN_C = K*UN_P/UN_rho;
                    UN_H = 0.5*UN_VV+UN_C/(K-1);

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

                    #if ROE != 4

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    #endif




                    #if ROE == 1 
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    temp = 0.5*(1+beta)*V;    // ---- U' ---- //
                    S = 0.5*sqrt(4*beta*C+V*V*(1-beta)*(1-beta));   // ---- C' ---- //

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = temp/S*(V_-UN_V);

                    temp3 = 0.5*(1-beta)*V*temp/S;



                    deltaU = (S-temp3-beta*fabs(V))*temp1+temp2;

                    deltaP = temp/S*(P_-UN_P)+(S-fabs(V)+temp3)*rho*(V_-UN_V);

                    #elif ROE == 2
                    
                    beta = max(VV/C,e);    // ---- theda ---- //
                    temp = 0.5*(1+beta)*V;    // ---- U' ---- //
                    S = 0.5*sqrt(4*beta*C+V*V*(1-beta)*(1-beta));   // ---- C' ---- //

                    theda_p = VV/C;
                    U_p = 0.5*(1+theda_p)*V;
                    C_p = 0.5*sqrt(4*C*theda_p+V*V*(1-theda_p)*(1-theda_p));

                    temp1 = (P_-UN_P)/rho/beta/C;

                    temp2 = U_p/S*(V_-UN_V);

                    temp3 = 0.5*(1-beta)*V*temp/S;


                    
                    temp4 = 0.5*(1-theda_p)*V*U_p/S;

                    deltaU = (S-temp3-beta*fabs(V))*temp1+temp2;


                    deltaP = U_p/S*(P_-UN_P)+(C_p-fabs(V)+temp4)*rho*(V_-UN_V);

                    #elif ROE == 3

                    S = sqrt(C);

                    theda_p = (fabs(U)+fabs(V)+fabs(W))/S;
                    
                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);

                    beta = max(temp, 1e-8);
                    
                    C_p = beta*S;

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = V/S*beta*(V_-UN_V);

                    deltaU = (S-fabs(V))*temp1+temp2;

                    deltaP = V/S*(P_-UN_P)+(C_p-fabs(V))*rho*(V_-UN_V);

                    #elif ROE == 4
                    
                    beta = sqrt(max(VV_/C_,UN_VV/UN_C));    // ---- theda ---- //

                    temp1 = 0.5*(UN_U+U_)+0.5*beta*(UN_U-U_);
                    temp2 = 0.5*(UN_V+V_)+0.5*beta*(UN_V-V_);
                    temp3 = 0.5*(UN_W+W_)+0.5*beta*(UN_W-W_);

                    U_ = 0.5*(U_+UN_U)+0.5*beta*(U_-UN_U);
                    V_ = 0.5*(V_+UN_V)+0.5*beta*(V_-UN_V);
                    W_ = 0.5*(W_+UN_W)+0.5*beta*(W_-UN_W);

                    UN_U = temp1;
                    UN_V = temp2;
                    UN_W = temp3;

                    temp5 = sqrt(UN_rho);
                    temp6 = sqrt(rho_);
                    
                    temp4 = temp5+temp6;

                    rho = sqrt(UN_rho*rho_);
                    U = (temp5*UN_U+temp6*U_)/temp4;
                    V = (temp5*UN_V+temp6*V_)/temp4;
                    W = (temp5*UN_W+temp6*W_)/temp4;
                    VV = U*U+V*V+W*W;
                    H = (temp5*UN_H+temp6*H_)/temp4;
                    C = (H-0.5*VV)*(K-1);    // ---- C*C ---- //
                    P = rho*C/K;

                    S = sqrt(C);

                    temp1 = (P_-UN_P)/rho/C;

                    temp2 = V/S*(V_-UN_V);

                    deltaU = (S-fabs(V))*temp1+temp2;

                    deltaP = V/S*(P_-UN_P)+(S-fabs(V))*rho*(V_-UN_V);

                    #elif ROE == 5
                    
                    S = sqrt(C);

                    theda_p = 0.5*(sqrt(VV_/C_)+sqrt(UN_VV/UN_C));

                    temp = theda_p*sqrt(4.0+(1.0-theda_p*theda_p)*(1.0-theda_p*theda_p))/(1.0+theda_p*theda_p);
                    
                    deltaP = V/S*(P_-UN_P)+temp*(S-fabs(V))*rho*(V_-UN_V);
                    
                    
                    // ==== deltaU_P + deltaU_RL ==== //            
                    temp1 = (P_-UN_P)/rho/C;

                    deltaU = ( 1.0-pow(temp,8) )*(S-fabs(V))*temp1;

                    temp2 = 0.5*V/S*(V_-UN_V);
                    
                    
                    // ==== KEASI ==== //
                    beta = sqrt(abs(V)/C);
                    
                    temp3 = beta*sqrt(4.0+(1.0-beta*beta)*(1.0-beta*beta))/(1.0+beta*beta);
                    
                    temp = max( 0.5*fabs(V_+UN_V)+(1.0-pow(temp3,8))*0.5*(V_-UN_V), 0.0 );


                    #endif

                    /* artificial viscosity */
                    
                    #if ROE == 0
                      
                      S = sqrt(C);
                        
                      /* artificial viscosity */
                      Fav1 = 1.0/S*(P_-UN_P);
                      Fav2 =   U/S*(P_-UN_P);
                      Fav3 =   V/S*(P_-UN_P);
                      Fav4 =   W/S*(P_-UN_P);
                      Fav5 =   H/S*(P_-UN_P);
                      
                    #elif ROE == 5 
                      
                    
                    Fav1 = Cdiss*temp*dU1+deltaU*rho  +temp2*(UN_rho+rho_);
                    Fav2 = Cdiss*temp*dU2+deltaU*rho*U+temp2*(UN_rho*UN_U+rho_*U_);
                    Fav3 = Cdiss*temp*dU3+deltaU*rho*V+temp2*(UN_rho*UN_V+rho_*V_)+deltaP;
                    Fav4 = Cdiss*temp*dU4+deltaU*rho*W+temp2*(UN_rho*UN_W+rho_*W_);
                    Fav5 = Cdiss*temp*dU5+deltaU*rho*H+temp2*(UN_rho*UN_H+rho_*H_)+deltaP*V;
                    
                    #else
                    
                    Fav1 = Cdiss*fabs(V)*dU1+deltaU*rho;
                    Fav2 = Cdiss*fabs(V)*dU2+deltaU*rho*U;
                    Fav3 = Cdiss*fabs(V)*dU3+deltaU*rho*V+deltaP;
                    Fav4 = Cdiss*fabs(V)*dU4+deltaU*rho*W;
                    Fav5 = Cdiss*fabs(V)*dU5+deltaU*rho*H+deltaP*V;
                    

                    #endif
                    
                    
                    /* inviscid fluxes */
                    inFy1i = 0.5*((UN_rho*UN_V+rho_*V_-Ep*Fav1));
                    inFy2i = 0.5*((UN_rho*UN_U*UN_V+rho_*U_*V_)-Ep*Fav2);
                    inFy3i = 0.5*((UN_rho*UN_V*UN_V+UN_P+rho_*V_*V_+P_)-Ep*Fav3);
                    inFy4i = 0.5*((UN_rho*UN_W*UN_V+rho_*W_*V_)-Ep*Fav4);
                    inFy5i = 0.5*((UN_V*(3.5*UN_P+0.5*UN_rho*UN_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Ep*Fav5);


                    dPyi = 0.5*(UN_P+P_);
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


                    dTx = dT_dx;
                    


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


                    dTxi = dT_dx;

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

                    dTy = dT_dy;

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

                    dTyi = dT_dy;

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

                    dTz = dT_dz;

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

                    dTzi = dT_dz;

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

                    
                    
                    rho = u1;
                    U = u2/u1;
                    V = u3/u1;
                    W = u4/u1;
                    VV = U*U+V*V+W*W;
                    P = (u5-0.5*rho*VV)*(K-1);
                    C = K*P/rho;
                    T = P/rho;
                    H = 0.5*VV+C/(K-1);
                    
                    
                    Rk1 = Fabs[icube][i][j][k][0]+Rp1+Rf1;
                    Rk2 = Fabs[icube][i][j][k][1]+Rp2+Rf2+vF2;
                    Rk3 = Fabs[icube][i][j][k][2]+Rp3+Rf3+vF3;
                    Rk4 = Fabs[icube][i][j][k][3]+Rp4+Rf4+vF4;
                    Rk5 = Fabs[icube][i][j][k][4]+Rp5+Rf5+vF5;
                    

                    /* preconditioning */

                    beta = max(VV/C,e);


                    temp = rho*K;

                    d11 = -( H*(K-1) + VV - K*(T+VV) )*beta;
                    d12 = -(K-1)*U*beta;
                    d13 = -(K-1)*V*beta;
                    d14 = -(K-1)*W*beta;
                    d15 = (K-1)*beta;

                    d21 = -U/rho;
                    d22 = 1/rho;
                    d23 = 0;
                    d24 = 0;
                    d25 = 0;

                    d31 = -V/rho;
                    d32 = 0;
                    d33 = 1/rho;
                    d34 = 0;
                    d35 = 0;

                    d41 = -W/rho;
                    d42 = 0;
                    d43 = 0;
                    d44 = 1/rho;
                    d45 = 0;

                    d51 = (K-1)*( VV + K*T*beta + (K-1)*VV*beta + H*(-1+beta-K*beta) )/temp;

                    d52 = -(K-1)*U*(1+(K-1)*beta)/temp;
                    d53 = -(K-1)*V*(1+(K-1)*beta)/temp;
                    d54 = -(K-1)*W*(1+(K-1)*beta)/temp;
                    d55 = (K-1)*(1+(K-1)*beta)/temp;
                    

                    if (IF == IFLUID) {

                        MR1 = d11*Rk1+d12*Rk2+d13*Rk3+d14*Rk4+d15*Rk5;
                        MR2 = d21*Rk1+d22*Rk2+d23*Rk3+d24*Rk4+d25*Rk5;
                        MR3 = d31*Rk1+d32*Rk2+d33*Rk3+d34*Rk4+d35*Rk5;
                        MR4 = d41*Rk1+d42*Rk2+d43*Rk3+d44*Rk4+d45*Rk5;
                        MR5 = d51*Rk1+d52*Rk2+d53*Rk3+d54*Rk4+d55*Rk5;

                    }
                    else {

                        MR1 = 0;
                        MR2 = 0;
                        MR3 = 0;
                        MR4 = 0;
                        MR5 = 0;

                    }

                    Rku1[icube][i][j][k][0] = MR1;
                    Rku1[icube][i][j][k][1] = MR2;
                    Rku1[icube][i][j][k][2] = MR3;
                    Rku1[icube][i][j][k][3] = MR4;
                    Rku1[icube][i][j][k][4] = MR5;
                    

                    Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);
                    Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);
                    Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);

                    Ux = ((beta+1)*fabs(U)+Sx)/2;
                    Uy = ((beta+1)*fabs(V)+Sy)/2;
                    Uz = ((beta+1)*fabs(W)+Sz)/2;



                    U_p = Ux*invXI+Uy*invET+Uz*invZT;

                    
                    #if defined(DTau)
                    
                    if(CFL_tau[icube][i][j][k] > -minimum) {
                        
                        deltaTau = DTau_CFL/max(U_p,1.0e-8);
                        
                        Rk1 = deltaTau*Rk1;
                        Rk2 = deltaTau*Rk2;
                        Rk3 = deltaTau*Rk3;
                        Rk4 = deltaTau*Rk4;
                        Rk5 = deltaTau*Rk5;
                        
                        UU1 = 0.1*max(rho*VV, 1.0e-9*P0);
                        UU2 = 2.0*max(2.0*sqrt(VV),1.0e-8);
                        UU3 = 2.0*max(2.0*sqrt(VV),1.0e-8);
                        UU4 = 2.0*max(2.0*sqrt(VV),1.0e-8);
                        UU5 = 0.1*T;
                        
                        LL1 = DTau_CFL;
                        LL2 = DTau_CFL;
                        LL3 = DTau_CFL;
                        LL4 = DTau_CFL;
                        LL5 = DTau_CFL;
                        
                        if( fabs(Rk1) > UU1 | fabs(Rk2) > UU2 | fabs(Rk3) > UU3 | fabs(Rk4) > UU4 | fabs(Rk5) > UU5 ) {
                            
                            LL1 = DTau_CFL*UU1/fabs(Rk1);
                            LL2 = DTau_CFL*UU2/fabs(Rk2);
                            LL3 = DTau_CFL*UU3/fabs(Rk3);
                            LL4 = DTau_CFL*UU4/fabs(Rk4);
                            LL5 = DTau_CFL*UU5/fabs(Rk5);
                            
                            CFL_tau[icube][i][j][k] = min(LL5,min(LL4,min(LL3,min(LL2,LL1))));
                            
                        }
                        else CFL_tau[icube][i][j][k] = DTau_CFL;
                        
                    }
                    
                    else CFL_tau[icube][i][j][k] = 0.3;
                    
                    deltaTau = CFL_tau[icube][i][j][k]/max(U_p,1.0e-8);
                    
                    temp = 1.0/deltaTau + U_p;
                    
                    d11 = 1.0 / (3.0*beta/(2.0*deltaT) + temp);

                    d22 = d33 = d44 = d55 = 1.0 / ( 3.0/(2.0*deltaT) +temp);

                    temp1 = 6*(K-1)*(beta-1)*deltaT*deltaTau*deltaTau;
                    temp2 = K*( 3*deltaTau+2*(1.0+U_p*deltaTau)*deltaT )*(3*beta*deltaTau+2*(1.0+U_p*deltaTau)*deltaT)*rho;
                    temp = -temp1/temp2;

                    d51 = temp;
                    
                    #elif defined(DTau_fix)
                    
                    temp = 1.0/deltaTau + U_p;
                    
                    d11 = 1.0 / (3.0*beta/(2.0*deltaT) + temp);

                    d22 = d33 = d44 = d55 = 1.0 / ( 3.0/(2.0*deltaT) +temp);

                    temp1 = 6*(K-1)*(beta-1)*deltaT*deltaTau*deltaTau;
                    temp2 = K*( 3*deltaTau+2*(1.0+U_p*deltaTau)*deltaT )*(3*beta*deltaTau+2*(1.0+U_p*deltaTau)*deltaT)*rho;
                    temp = -temp1/temp2;

                    d51 = temp;
                    
                    #elif defined(NODTau)
                    
                    d11 = 2*deltaT/(3*beta+2*U_p*deltaT);

                    d22 = d33 = d44 = d55 = 2*deltaT/(3+2*U_p*deltaT);

                    temp1 = 6*(K-1)*(beta-1)*deltaT;
                    temp2 = K*(3+2*U_p*deltaT)*(3*beta+2*U_p*deltaT)*rho;
                    temp = -temp1/temp2;

                    d51 = temp;
                    
                    #elif defined(DTauCAA)

                    deltaTau = DTau_CFL/max(U_p,1.0e-8);
                    
                    Rk1 = deltaTau*Rk1;
                    Rk2 = deltaTau*Rk2;
                    Rk3 = deltaTau*Rk3;
                    Rk4 = deltaTau*Rk4;
                    Rk5 = deltaTau*Rk5;
                    
                    UU1 = 0.1*max(rho*VV, 1.0e-9*P0);
                    UU2 = 2.0*max(2.0*sqrt(VV),1.0e-8);
                    UU3 = 2.0*max(2.0*sqrt(VV),1.0e-8);
                    UU4 = 2.0*max(2.0*sqrt(VV),1.0e-8);
                    UU5 = 0.1*T;
                    
                    LL1 = DTau_CFL;
                    LL2 = DTau_CFL;
                    LL3 = DTau_CFL;
                    LL4 = DTau_CFL;
                    LL5 = DTau_CFL;

                    ib = 0;
                    
                    if( fabs(Rk1) > UU1 | fabs(Rk2) > UU2 | fabs(Rk3) > UU3 | fabs(Rk4) > UU4 | fabs(Rk5) > UU5 ) {
                        
                        LL1 = DTau_CFL*UU1/fabs(Rk1);
                        LL2 = DTau_CFL*UU2/fabs(Rk2);
                        LL3 = DTau_CFL*UU3/fabs(Rk3);
                        LL4 = DTau_CFL*UU4/fabs(Rk4);
                        LL5 = DTau_CFL*UU5/fabs(Rk5);

                        ib = 1;
                        
                        CFL_tau[icube][i][j][k] = 2.0*ib;
                        
                    }
                    
                    if (ib == 0) {

                        d11 = 2*deltaT/(3*beta+2*U_p*deltaT);

                        d22 = d33 = d44 = d55 = 2*deltaT/(3+2*U_p*deltaT);

                        temp1 = 6*(K-1)*(beta-1)*deltaT;
                        temp2 = K*(3+2*U_p*deltaT)*(3*beta+2*U_p*deltaT)*rho;
                        temp = -temp1/temp2;

                        d51 = temp;

                    }
                    else {

                        if (CFL_tau[icube][i][j][k]>0)	CFL_tau[icube][i][j][k] = min(LL5,min(LL4,min(LL3,min(LL2,LL1))));
                        else   CFL_tau[icube][i][j][k]=0.3;
                        
                        temp = 1.0/deltaTau + U_p;
                        
                        d11 = 1.0 / (3.0*beta/(2.0*deltaT) + temp);

                        d22 = d33 = d44 = d55 = 1.0 / ( 3.0/(2.0*deltaT) +temp);

                        temp1 = 6*(K-1)*(beta-1)*deltaT*deltaTau*deltaTau;
                        temp2 = K*( 3*deltaTau+2*(1.0+U_p*deltaTau)*deltaT )*(3*beta*deltaTau+2*(1.0+U_p*deltaTau)*deltaT)*rho;
                        temp = -temp1/temp2;

                        d51 = temp;
                        
                    }

                    
                    Residual1[icube][i][j][k][4] = ib;
                    
                    #endif

                    
                    U1p1[icube][i][j][k][0] = d11*MR1;
                    U1p1[icube][i][j][k][1] = d22*MR2;
                    U1p1[icube][i][j][k][2] = d33*MR3;
                    U1p1[icube][i][j][k][3] = d44*MR4;
                    U1p1[icube][i][j][k][4] = d51*MR1+d55*MR5;

                    
                    if (IF < IFLUID) {
                        
                        Residual1[icube][i][j][k][0] = ( dTxi-dTx+dTyi-dTy+dTzi-dTz )/dx*dx*dy*dz;
                        
                        Residual1[icube][i][j][k][1] = -( -(dPxi-dPx)/dx )*dx*dy*dz-( vF2 )*dx*dy*dz;

                        Residual1[icube][i][j][k][2] = -( -(dPyi-dPy)/dy )*dx*dy*dz-( vF3 )*dx*dy*dz;
                        
                        Residual1[icube][i][j][k][3] = -( -(dPzi-dPz)/dz )*dx*dy*dz-( vF4 )*dx*dy*dz;

                    }

                    
                }    // ---- for (k = 2; k <= nz; k++) { ---- //
            }
        }

    }
    
    
    #pragma omp barrier
    


    
    for (isweep = 1; isweep < nsweep+1; isweep++) {
        
        
        
        
        BCM_Interface(myid,ncube,

        mPI_Nadj,

        ncpu_bs, ncpu_eq, ncpu_sb,
        max_nei_bs,max_nei_eq,max_nei_bs,

        nadjX_eq, nadjY_eq, nadjZ_eq,
        nadjX_bs_plus, nadjX_sb_plus, nadjX_bs_minus, nadjX_sb_minus,
        nadjY_bs_plus, nadjY_sb_plus, nadjY_bs_minus, nadjY_sb_minus,
        nadjZ_bs_plus, nadjZ_sb_plus, nadjZ_bs_minus, nadjZ_sb_minus,

        rank_map,

        MPI_cpu, MPI_cube, MPI_cpu_adj, MPI_cube_adj, MPI_direction, MPI_interface,
        neighbor_cpu_eq, Ncube_Ncpu_eq, neighbor_cpu_sb, Ncube_Ncpu_sb, neighbor_cpu_bs, Ncube_Ncpu_bs,
        Scube_Ncpu_eq, Rcube_Ncpu_eq, send_data_curr_eq, recv_data_curr_eq, send_data_neig_eq, recv_data_neig_eq, Sdir_eq, Rdir_eq,  
        Scube_Ncpu_sb, Rcube_Ncpu_sb, send_data_curr_sb, recv_data_curr_sb, send_data_neig_sb, recv_data_neig_sb, Sdir_sb, Rdir_sb,
        Scube_Ncpu_bs, Rcube_Ncpu_bs, send_data_curr_bs, recv_data_curr_bs, send_data_neig_bs, recv_data_neig_bs, Sdir_bs, Rdir_bs,
        ist_eq,ist_sb,ist_bs, adjN_sb, RadjN_bs, SadjN_bs,

        csl, 
        adj_number, 
        adjX_eq, adjY_eq, adjZ_eq,
        adjX_bs_plus, adjX_sb_plus, adjX_bs_minus, adjX_sb_minus,
        adjY_bs_plus, adjY_sb_plus, adjY_bs_minus, adjY_sb_minus,
        adjZ_bs_plus, adjZ_sb_plus, adjZ_bs_minus, adjZ_sb_minus,
        U1p1);


        
        // ---- data transfer ---- //
        
        
        #pragma omp parallel for private(\
            IF,\
            dx,dy,dz,invXI,invET,invZT,\
            i,j,k,\
            rho,U,V,W,VV,P,C,T,h,H,\
            Ux,Uy,Uz,\
            beta,Sx,Sy,Sz,U_p,\
            temp,temp1,temp2,\
            d11,d12,d13,d14,d15,\
            d21,d22,d23,d24,d25,\
            d31,d32,d33,d34,d35,\
            d41,d42,d43,d44,d45,\
            d51,d52,d53,d54,d55,\
            LL1,LL2,LL3,LL4,LL5,\
            UU1,UU2,UU3,UU4,UU5,\
            Aplus11,Aplus12,\
            Aplus21,Aplus22,\
            Aplus33,\
            Aplus44,\
            Aplus51,Aplus52,Aplus55,\
            Bplus11,Bplus13,\
            Bplus22,\
            Bplus31,Bplus33,\
            Bplus41,Bplus44,\
            Bplus51,Bplus53,Bplus55,\
            Cplus11,Cplus14,\
            Cplus22,\
            Cplus33,\
            Cplus41,Cplus44,\
            Cplus51,Cplus54,Cplus55,\
            Amius11,Amius22,Amius33,Amius44,Amius55,\
            Bmius11,Bmius22,Bmius33,Bmius44,Bmius55,\
            Cmius11,Cmius22,Cmius33,Cmius44,Cmius55,\
            deltaTau\
            )

        for (icube = 1; icube < ncube; icube++) {  

            dx = dy = dz = cube_size[icube]/NcubeX;
            invXI = invET = invZT = 1./dx;

            for (i = n_buffer; i <= nx; i++) {
                for (j = n_buffer; j <= ny; j++) {
                    for (k = n_buffer; k <= nz; k++) {

                        IF = FWS[icube][i][j][k];
                        
                        rho = U1_[icube][i-1][j][k][0];
                        U = U1_[icube][i-1][j][k][1]/rho;
                        V = U1_[icube][i-1][j][k][2]/rho;
                        W = U1_[icube][i-1][j][k][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i-1][j][k][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;
                        H = 0.5*VV+C/(K-1.0);
                        
                        beta = max(VV/C,e);
                        
                        Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);
                        Ux = 0.5*((beta+1)*fabs(U)+Sx) + 2*K*mu_L/Pr_L/rho/dx;
                        
                        d11 = fabs(Ux);
                        d22 = d11;
                        d33 = d11;
                        d44 = d11;
                        d55 = d11;

                        Aplus11 = 0.5*(beta*U+d11);
                        Aplus12 = 0.5*P*K*beta;

                        Aplus21 = 0.5/rho;
                        Aplus22 = 0.5*(U+d22);

                        Aplus33 = 0.5*(U+d33);

                        Aplus44 = 0.5*(U+d44);

                        Aplus51 = 0.5*U*(K-1)*(beta-1)/rho/K;
                        Aplus52 = 0.5*beta*T*(K-1);
                        Aplus55 = 0.5*(U+d55);
                        
                        
                        
                        
                        rho = U1_[icube][i][j-1][k][0];
                        U = U1_[icube][i][j-1][k][1]/rho;
                        V = U1_[icube][i][j-1][k][2]/rho;
                        W = U1_[icube][i][j-1][k][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i][j-1][k][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;
                        H = 0.5*VV+C/(K-1.0);
                        
                        beta = max(VV/C,e);
                        
                        Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);
                        Uy = 0.5*((beta+1)*fabs(V)+Sy) + 2*K*mu_L/Pr_L/rho/dy;
                        
                        d11 = fabs(Uy);
                        d22 = d11;
                        d33 = d11;
                        d44 = d11;
                        d55 = d11;

                        Bplus11 = 0.5*(beta*V+d11);
                        Bplus13 = 0.5*P*K*beta;

                        Bplus22 = 0.5*(V+d22);

                        Bplus31 = 0.5/rho;
                        Bplus33 = 0.5*(V+d33);

                        Bplus44 = 0.5*(V+d44);

                        Bplus51 = 0.5*V*(K-1)*(beta-1)/rho/K;
                        Bplus53 = 0.5*beta*T*(K-1);
                        Bplus55 = 0.5*(V+d55);
                        
                        
                        
                        
                        rho = U1_[icube][i][j][k-1][0];
                        U = U1_[icube][i][j][k-1][1]/rho;
                        V = U1_[icube][i][j][k-1][2]/rho;
                        W = U1_[icube][i][j][k-1][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i][j][k-1][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;
                        H = 0.5*VV+C/(K-1.0);
                        
                        beta = max(VV/C,e);
                        
                        Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);
                        Uz = 0.5*((beta+1)*fabs(W)+Sz) + 2*K*mu_L/Pr_L/rho/dz;
                        
                        d11 = fabs(Uz);
                        d22 = d11;
                        d33 = d11;
                        d44 = d11;
                        d55 = d11;

                        Cplus11 = 0.5*(beta*W+d11);
                        Cplus14 = 0.5*P*K*beta;

                        Cplus22 = 0.5*(W+d22);

                        Cplus33 = 0.5*(W+d33);

                        Cplus41 = 0.5/rho;
                        Cplus44 = 0.5*(W+d44);

                        Cplus51 = 0.5*W*(K-1)*(beta-1)/rho/K;
                        Cplus54 = 0.5*beta*T*(K-1);
                        Cplus55 = 0.5*(W+d55);
                        
                        
                        LL1 = (-Aplus11*U1p1[icube][i-1][j][k][0]       \
                        -Aplus12*U1p1[icube][i-1][j][k][1])/dx + \
                        (-Bplus11*U1p1[icube][i][j-1][k][0]       \
                        -Bplus13*U1p1[icube][i][j-1][k][2])/dy + \
                        (-Cplus11*U1p1[icube][i][j][k-1][0]       \
                        -Cplus14*U1p1[icube][i][j][k-1][3])/dz;
                        
                        LL2 = (-Aplus21*U1p1[icube][i-1][j][k][0]       \
                        -Aplus22*U1p1[icube][i-1][j][k][1])/dx + \
                        (-Bplus22*U1p1[icube][i][j-1][k][1])/dy + \
                        (-Cplus22*U1p1[icube][i][j][k-1][1])/dz;
                        
                        LL3 = (-Aplus33*U1p1[icube][i-1][j][k][2])/dx + \
                        (-Bplus31*U1p1[icube][i][j-1][k][0]       \
                        -Bplus33*U1p1[icube][i][j-1][k][2])/dy + \
                        (-Cplus33*U1p1[icube][i][j][k-1][2])/dz;
                        
                        LL4 = (-Aplus44*U1p1[icube][i-1][j][k][3])/dx + \
                        (-Bplus44*U1p1[icube][i][j-1][k][3])/dy + \
                        (-Cplus41*U1p1[icube][i][j][k-1][0]       \
                        -Cplus44*U1p1[icube][i][j][k-1][3])/dz;     
                        
                        LL5 = (-Aplus51*U1p1[icube][i-1][j][k][0]       \
                        -Aplus52*U1p1[icube][i-1][j][k][1]       \
                        -Aplus55*U1p1[icube][i-1][j][k][4])/dx + \
                        (-Bplus51*U1p1[icube][i][j-1][k][0]       \
                        -Bplus53*U1p1[icube][i][j-1][k][2]       \
                        -Bplus55*U1p1[icube][i][j-1][k][4])/dy + \
                        (-Cplus51*U1p1[icube][i][j][k-1][0]       \
                        -Cplus54*U1p1[icube][i][j][k-1][3]       \
                        -Cplus55*U1p1[icube][i][j][k-1][4])/dz;
                        
                        // =================================  Lower part end  ================================= //



                        rho = U1_[icube][i+1][j][k][0];
                        U = U1_[icube][i+1][j][k][1]/rho;
                        V = U1_[icube][i+1][j][k][2]/rho;
                        W = U1_[icube][i+1][j][k][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i+1][j][k][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;
                        H = 0.5*VV+C/(K-1.0);
                        
                        beta = max(VV/C,e);
                        
                        Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);
                        Ux = 0.5*((beta+1)*fabs(U)+Sx) + 2*K*mu_L/Pr_L/rho/dx;
                        
                        d11 = fabs(Ux);
                        d22 = d11;
                        d33 = d11;
                        d44 = d11;
                        d55 = d11;

                        Aplus11 = 0.5*(beta*U+d11);
                        Aplus12 = 0.5*P*K*beta;

                        Aplus21 = 0.5/rho;
                        Aplus22 = 0.5*(U+d22);

                        Aplus33 = 0.5*(U+d33);

                        Aplus44 = 0.5*(U+d44);

                        Aplus51 = 0.5*U*(K-1)*(beta-1)/rho/K;
                        Aplus52 = 0.5*beta*T*(K-1);
                        Aplus55 = 0.5*(U+d55);
                        
                        Amius11 = 0.5*(beta*U-d11);
                        Amius22 = 0.5*(U-d22);
                        Amius33 = 0.5*(U-d33);
                        Amius44 = 0.5*(U-d44);
                        Amius55 = 0.5*(U-d55);
                        
                        
                        
                        rho = U1_[icube][i][j+1][k][0];
                        U = U1_[icube][i][j+1][k][1]/rho;
                        V = U1_[icube][i][j+1][k][2]/rho;
                        W = U1_[icube][i][j+1][k][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i][j+1][k][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;
                        H = 0.5*VV+C/(K-1.0);
                        
                        beta = max(VV/C,e);
                        
                        Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);
                        Uy = 0.5*((beta+1)*fabs(V)+Sy) + 2*K*mu_L/Pr_L/rho/dy;
                        
                        d11 = fabs(Uy);
                        d22 = d11;
                        d33 = d11;
                        d44 = d11;
                        d55 = d11;

                        Bplus11 = 0.5*(beta*V+d11);
                        Bplus13 = 0.5*P*K*beta;

                        Bplus22 = 0.5*(V+d22);

                        Bplus31 = 0.5/rho;
                        Bplus33 = 0.5*(V+d33);

                        Bplus44 = 0.5*(V+d44);

                        Bplus51 = 0.5*V*(K-1)*(beta-1)/rho/K;
                        Bplus53 = 0.5*beta*T*(K-1);
                        Bplus55 = 0.5*(V+d55);
                        
                        Bmius11 = 0.5*(beta*V-d11);
                        Bmius22 = 0.5*(V-d22);
                        Bmius33 = 0.5*(V-d33);
                        Bmius44 = 0.5*(V-d44);
                        Bmius55 = 0.5*(V-d55);
                        
                        
                        
                        
                        rho = U1_[icube][i][j][k+1][0];
                        U = U1_[icube][i][j][k+1][1]/rho;
                        V = U1_[icube][i][j][k+1][2]/rho;
                        W = U1_[icube][i][j][k+1][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i][j][k+1][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;
                        H = 0.5*VV+C/(K-1.0);
                        
                        beta = max(VV/C,e);
                        
                        Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);
                        Uz = 0.5*((beta+1)*fabs(W)+Sz) + 2*K*mu_L/Pr_L/rho/dz;
                        
                        d11 = fabs(Uz);
                        d22 = d11;
                        d33 = d11;
                        d44 = d11;
                        d55 = d11;

                        Cplus11 = 0.5*(beta*W+d11);
                        Cplus14 = 0.5*P*K*beta;

                        Cplus22 = 0.5*(W+d22);

                        Cplus33 = 0.5*(W+d33);

                        Cplus41 = 0.5/rho;
                        Cplus44 = 0.5*(W+d22);

                        Cplus51 = 0.5*W*(K-1)*(beta-1)/rho/K;
                        Cplus54 = 0.5*beta*T*(K-1);
                        Cplus55 = 0.5*(W+d55);
                        
                        Cmius11 = 0.5*(beta*W-d11);
                        Cmius22 = 0.5*(W-d22);
                        Cmius33 = 0.5*(W-d33);
                        Cmius44 = 0.5*(W-d44);
                        Cmius55 = 0.5*(W-d55);
                        
                        
                        UU1 = (+Amius11*U1p1[icube][i+1][j][k][0]       \
                        +Aplus12*U1p1[icube][i+1][j][k][1])/dx + \
                        (+Bmius11*U1p1[icube][i][j+1][k][0]       \
                        +Bplus13*U1p1[icube][i][j+1][k][2])/dy + \
                        (+Cmius11*U1p1[icube][i][j][k+1][0]       \
                        +Cplus14*U1p1[icube][i][j][k+1][3])/dz;
                        
                        UU2 = (+Aplus21*U1p1[icube][i+1][j][k][0]       \
                        +Amius22*U1p1[icube][i+1][j][k][1])/dx + \
                        (+Bmius22*U1p1[icube][i][j+1][k][1])/dy + \
                        (+Cmius22*U1p1[icube][i][j][k+1][1])/dz;
                        
                        UU3 = (+Amius33*U1p1[icube][i+1][j][k][2])/dx + \
                        (+Bplus31*U1p1[icube][i][j+1][k][0]       \
                        +Bmius33*U1p1[icube][i][j+1][k][2])/dy + \
                        (+Cmius33*U1p1[icube][i][j][k+1][2])/dz;
                        
                        UU4 = (+Amius44*U1p1[icube][i+1][j][k][3])/dx + \
                        (+Bmius44*U1p1[icube][i][j+1][k][3])/dy + \
                        (+Cplus41*U1p1[icube][i][j][k+1][0]       \
                        +Cmius44*U1p1[icube][i][j][k+1][3])/dz;     
                        
                        UU5 = (+Aplus51*U1p1[icube][i+1][j][k][0]       \
                        +Aplus52*U1p1[icube][i+1][j][k][1]       \
                        +Amius55*U1p1[icube][i+1][j][k][4])/dx + \
                        (+Bplus51*U1p1[icube][i][j+1][k][0]       \
                        +Bplus53*U1p1[icube][i][j+1][k][2]       \
                        +Bmius55*U1p1[icube][i][j+1][k][4])/dy + \
                        (+Cplus51*U1p1[icube][i][j][k+1][0]       \
                        +Cplus54*U1p1[icube][i][j][k+1][3]       \
                        +Cmius55*U1p1[icube][i][j][k+1][4])/dz;
                        
                        
                        rho = U1_[icube][i][j][k][0];
                        U = U1_[icube][i][j][k][1]/rho;
                        V = U1_[icube][i][j][k][2]/rho;
                        W = U1_[icube][i][j][k][3]/rho;
                        VV = U*U+V*V+W*W;
                        P = (U1_[icube][i][j][k][4]-0.5*rho*VV)*(K-1);
                        T = P/rho;
                        C = K*P/rho;

                        /* preconditioning */

                        beta = max(VV/C,e);
                        
                        Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C);
                        Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C);
                        Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C);

                        Ux = ((beta+1)*fabs(U)+Sx)/2;
                        Uy = ((beta+1)*fabs(V)+Sy)/2;
                        Uz = ((beta+1)*fabs(W)+Sz)/2;

                        

                        U_p = Ux*invXI+Uy*invET+Uz*invZT;

                        #if defined(DTau)
                        
                        deltaTau = CFL_tau[icube][i][j][k]/max(U_p,1.0e-8);
                        
                        temp = 1.0/deltaTau + U_p;
                        
                        d11 = 1.0 / (3.0*beta/(2.0*deltaT) + temp);

                        d22 = d33 = d44 = d55 = 1.0 / ( 3.0/(2.0*deltaT) +temp);

                        temp1 = 6*(K-1)*(beta-1)*deltaT*deltaTau*deltaTau;
                        temp2 = K*( 3*deltaTau+2*(1.0+U_p*deltaTau)*deltaT )*(3*beta*deltaTau+2*(1.0+U_p*deltaTau)*deltaT)*rho;
                        temp = -temp1/temp2;

                        d51 = temp;
                        
                        #elif defined(DTau_fix)
                        
                        temp = 1.0/deltaTau + U_p;
                        
                        d11 = 1.0 / (3.0*beta/(2.0*deltaT) + temp);

                        d22 = d33 = d44 = d55 = 1.0 / ( 3.0/(2.0*deltaT) +temp);

                        temp1 = 6*(K-1)*(beta-1)*deltaT*deltaTau*deltaTau;
                        temp2 = K*( 3*deltaTau+2*(1.0+U_p*deltaTau)*deltaT )*(3*beta*deltaTau+2*(1.0+U_p*deltaTau)*deltaT)*rho;
                        temp = -temp1/temp2;

                        d51 = temp;
                        
                        #elif defined(NODTau)
                        
                        d11 = 2*deltaT/(3*beta+2*U_p*deltaT);

                        d22 = d33 = d44 = d55 = 2*deltaT/(3+2*U_p*deltaT);

                        temp1 = 6*(K-1)*(beta-1)*deltaT;
                        temp2 = K*(3+2*U_p*deltaT)*(3*beta+2*U_p*deltaT)*rho;
                        temp = -temp1/temp2;

                        d51 = temp;
                        
                        #elif defined(DTauCAA)
                        
                        if(Residual1[icube][i][j][k][4]<0.5) {
                            
                            d11 = 2*deltaT/(3*beta+2*U_p*deltaT);

                            d22 = d33 = d44 = d55 = 2*deltaT/(3+2*U_p*deltaT);

                            temp1 = 6*(K-1)*(beta-1)*deltaT;
                            temp2 = K*(3+2*U_p*deltaT)*(3*beta+2*U_p*deltaT)*rho;
                            temp = -temp1/temp2;

                            d51 = temp;
                            
                        }
                        else {
                            
                            deltaTau = CFL_tau[icube][i][j][k]/max(U_p,1.0e-8);
                            
                            temp = 1.0/deltaTau + U_p;
                            
                            d11 = 1.0 / (3.0*beta/(2.0*deltaT) + temp);

                            d22 = d33 = d44 = d55 = 1.0 / ( 3.0/(2.0*deltaT) +temp);

                            temp1 = 6*(K-1)*(beta-1)*deltaT*deltaTau*deltaTau;
                            temp2 = K*( 3*deltaTau+2*(1.0+U_p*deltaTau)*deltaT )*(3*beta*deltaTau+2*(1.0+U_p*deltaTau)*deltaT)*rho;
                            temp = -temp1/temp2;

                            d51 = temp;
                            
                        }
                        
                        #endif

                        
                        
                        if (IF == IFLUID) {

                            U1p2[icube][i][j][k][0] = d11*(Rku1[icube][i][j][k][0] - LL1 - UU1);
                            U1p2[icube][i][j][k][1] = d22*(Rku1[icube][i][j][k][1] - LL2 - UU2);
                            U1p2[icube][i][j][k][2] = d33*(Rku1[icube][i][j][k][2] - LL3 - UU3);
                            U1p2[icube][i][j][k][3] = d44*(Rku1[icube][i][j][k][3] - LL4 - UU4);
                            U1p2[icube][i][j][k][4] = d51*(Rku1[icube][i][j][k][0] - LL1 - UU1) + \
                            d55*(Rku1[icube][i][j][k][4] - LL5 - UU5);
                            
                        }
                        else {

                            U1p2[icube][i][j][k][0] = 0;
                            U1p2[icube][i][j][k][1] = 0;
                            U1p2[icube][i][j][k][2] = 0;
                            U1p2[icube][i][j][k][3] = 0;
                            U1p2[icube][i][j][k][4] = 0;

                        }
                        
                        
                    }
                }
            }
            
        }    // ---- for (icube = 1; icube < ncube; icube++) ---- //

        #pragma omp barrier

        #pragma omp parallel for private(\
            i,j,k\
            )
        for (icube = 1; icube < ncube; icube++) {  
            
            for (i = n_buffer ; i < nxx; i++) {
                for (j = n_buffer; j < nyy; j++) {
                    for (k = n_buffer; k < nzz; k++) {

                        U1p1[icube][i][j][k][0] = U1p2[icube][i][j][k][0];
                        U1p1[icube][i][j][k][1] = U1p2[icube][i][j][k][1];
                        U1p1[icube][i][j][k][2] = U1p2[icube][i][j][k][2];
                        U1p1[icube][i][j][k][3] = U1p2[icube][i][j][k][3];
                        U1p1[icube][i][j][k][4] = U1p2[icube][i][j][k][4];
                        
                    }
                }
            }
            
        }

        
        #pragma omp barrier
        

        
    }    // ---- for (isweep = 1; isweep < nsweep+1; isweep++) ---- //


#pragma omp parallel for private(IF,i,j,k,rho,P,U,V,W,T,rhoold,Uold,Vold,Wold,VVold,Pold,Told)reduction(+:e1,e2,e3,e4,e5)
    
    for (icube = 1; icube < ncube; icube++) {  

        for (i = n_buffer ; i < nxx; i++) {
            for (j = n_buffer; j < nyy; j++) {
                for (k = n_buffer; k < nzz; k++) {

                    IF = FWS[icube][i][j][k];
                    
                    /* flux parameter */
                    rhoold = U1_[icube][i][j][k][0];
                    Uold = U1_[icube][i][j][k][1]/rhoold;
                    Vold = U1_[icube][i][j][k][2]/rhoold;
                    Wold = U1_[icube][i][j][k][3]/rhoold;
                    VVold = Uold*Uold+Vold*Vold+Wold*Wold;
                    Pold = (U1_[icube][i][j][k][4]-0.5*rhoold*VVold)*(K-1);
                    Told = Pold/rhoold;

                    P = Pold + U1p1[icube][i][j][k][0];
                    U = Uold + U1p1[icube][i][j][k][1];
                    V = Vold + U1p1[icube][i][j][k][2];
                    W = Wold + U1p1[icube][i][j][k][3];
                    T = Told + U1p1[icube][i][j][k][4];
                    rho = P/T;

                    if (IF == IFLUID) {

                        U1_[icube][i][j][k][0] = rho;
                        U1_[icube][i][j][k][1] = rho*U;
                        U1_[icube][i][j][k][2] = rho*V;
                        U1_[icube][i][j][k][3] = rho*W;
                        U1_[icube][i][j][k][4] = P/(K-1.0)+0.5*rho*(U*U+V*V+W*W);


                        e1 = e1+(P-Pold)*(P-Pold);
                        e2 = e2+(U-Uold)*(U-Uold);
                        e3 = e3+(V-Vold)*(V-Vold);	
                        e4 = e4+(W-Wold)*(W-Wold);
                        e5 = e5+(T-Told)*(T-Told);
                        
                    }
                    
                    else {
                        
                        e1 = 0.0;
                        e2 = 0.0;
                        e3 = 0.0;
                        e4 = 0.0;
                        e5 = 0.0;
                        
                    }


                }
            }
        }
        
        
    }
    
    
#pragma omp barrier
    


    
    for (icube = 1; icube < ncube; icube++) {  
        for (i = n_buffer; i <= nx; i++) {
            for (j = n_buffer; j <= ny; j++) { 
                for (k = n_buffer; k <= nz; k++) { 

                    e6 = e6+Residual1[icube][i][j][k][0];  // Averaged Nusselt number //
                    e7 = e7+Residual1[icube][i][j][k][1];  // Cd //
                    e8 = e8+Residual1[icube][i][j][k][2];  // Cl //
                    
                }
            }
        }
    }

    double DN = 1./(ncube*NcubeX*NcubeY*NcubeZ);

    e1 = sqrt(e1)*DN;
    e2 = sqrt(e2)*DN;
    e3 = sqrt(e3)*DN;
    e4 = sqrt(e4)*DN;
    e5 = sqrt(e5)*DN;
    
    MPI_Comm comm;
    comm=MPI_COMM_WORLD;

    MPI_Allreduce ((void*)&e1,(void*)&er[1],1,MPI_DOUBLE,MPI_SUM,comm );
    MPI_Allreduce ((void*)&e2,(void*)&er[2],1,MPI_DOUBLE,MPI_SUM,comm );
    MPI_Allreduce ((void*)&e3,(void*)&er[3],1,MPI_DOUBLE,MPI_SUM,comm );
    MPI_Allreduce ((void*)&e4,(void*)&er[4],1,MPI_DOUBLE,MPI_SUM,comm );
    MPI_Allreduce ((void*)&e5,(void*)&er[5],1,MPI_DOUBLE,MPI_SUM,comm );
    
    er[1] = er[1]/np;
    er[2] = er[2]/np;
    er[3] = er[3]/np;
    er[4] = er[4]/np;
    er[5] = er[5]/np;
    
    
    MPI_Allreduce ((void*)&e6,(void*)&er[6], 1, MPI_DOUBLE, MPI_SUM, comm );
    MPI_Allreduce ((void*)&e7,(void*)&er[7], 1, MPI_DOUBLE, MPI_SUM, comm );
    MPI_Allreduce ((void*)&e8,(void*)&er[8], 1, MPI_DOUBLE, MPI_SUM, comm );


    
    
    
}