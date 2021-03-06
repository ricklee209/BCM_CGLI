



#include "Resolution.h"

extern int Ncube;    
extern int N_wallcube; 

extern int MPI_Nadj;

extern int Max_nei_eq;
extern int Max_nei_sb;
extern int Max_nei_bs;

extern int Ncpu_eq;
extern int Ncpu_sb;
extern int Ncpu_bs;

void BCM_Interface
(
// ================================================================================ //
int myid,
int ncube,

int mPI_Nadj,

int ncpu_bs,
int ncpu_eq, 
int ncpu_sb,

int max_nei_bs,
int max_nei_eq,
int max_nei_sb,

int nadjX_bs, 
int nadjX_eq, 
int nadjX_sb,

int nadjY_bs, 
int nadjY_eq, 
int nadjY_sb,

int nadjZ_bs, 
int nadjZ_eq, 
int nadjZ_sb,

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
int (*adjN_bs) = new int[Max_nei_bs],


int (*csl) = new int[Ncube],

int (*adj_number)[5][7] = new int[Ncube][5][7],

int (*adjX_bs) = new int[Ncube],
int (*adjX_eq) = new int[Ncube],
int (*adjX_sb) = new int[Ncube],

int (*adjY_bs) = new int[Ncube],
int (*adjY_eq) = new int[Ncube],
int (*adjY_sb) = new int[Ncube],

int (*adjZ_bs) = new int[Ncube],
int (*adjZ_eq) = new int[Ncube],
int (*adjZ_sb) = new int[Ncube],


double (*U1_)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U2_)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U3_)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U4_)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U5_)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]

// ================================================================================ //
);

