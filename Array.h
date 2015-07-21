



#include "Resolution.h"

//extern int Ncube;
//extern int N_wallcube;


double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*U1p1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*U1p2)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*Residual1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*Rku1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];

double (*cube_size) = new double[Ncube];
int (*csl) = new int[Ncube];

double (*Xcube) = new double[Ncube];
double (*Ycube) = new double[Ncube];
double (*Zcube) = new double[Ncube];

double (*Xcnt)[X_size] = new double[Ncube][X_size];
double (*Ycnt)[Y_size] = new double[Ncube][Y_size];
double (*Zcnt)[Z_size] = new double[Ncube][Z_size];

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size];
// ---- F:fluid -1 ---- W:wall 0 ---- S:Sloid 1 ---- //

int (*adj_number)[5][7] = new int[Ncube][5][7];
int (*wallcube) = new int[N_wallcube];

int (*adjX_eq) = new int[Ncube];
int (*adjY_eq) = new int[Ncube];
int (*adjZ_eq) = new int[Ncube];

int (*adjX_bs_plus) = new int[Ncube];
int (*adjX_sb_plus) = new int[Ncube];

int (*adjX_bs_minus) = new int[Ncube];
int (*adjX_sb_minus) = new int[Ncube];

int (*adjY_bs_plus) = new int[Ncube];
int (*adjY_sb_plus) = new int[Ncube];

int (*adjY_bs_minus) = new int[Ncube];
int (*adjY_sb_minus) = new int[Ncube];

int (*adjZ_bs_plus) = new int[Ncube];
int (*adjZ_sb_plus) = new int[Ncube];

int (*adjZ_bs_minus) = new int[Ncube];
int (*adjZ_sb_minus) = new int[Ncube];


/**** MPI_data_communication ****/
double (*send_data_cur) = new double[2*X_size*X_size*5];
double (*recv_data_cur) = new double[2*X_size*X_size*5];

double (*send_data_adj) = new double[2*X_size*X_size*5];
double (*recv_data_adj) = new double[2*X_size*X_size*5];


/**** Statistical data ****/
double (*Pall)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size];

double (*VVall)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size];


/**** Absorbing force term ****/
double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim];


/**** artificial dissipation coefficient ****/
double (*Wig)[3][X_size][Y_size][Z_size] = new double[Ncube][3][X_size][Y_size][Z_size];

