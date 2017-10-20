



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


extern int NBC_minus;
extern int NBC_plus;


void BCM_Grid 
(
// =================================================== //
int myid,

int ncube,

int n_wallcube,

double (*cube_size) = new double[Ncube],
int (*csl) = new int[Ncube],

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

int (*adj_number)[5][7] = new int[Ncube][5][7],

int (*wallcube) = new int[Ncube]
// ================================================== //
);



 void BCM_FWS_Interface
 (
 // ================================================================================ //
 int myid,
 int ncube,

 int mPI_Nadj,

 int ncpu_eq,

 int max_nei_eq,

 int nadjX_eq, 
 int nadjY_eq, 
 int nadjZ_eq,

 int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

 int (*MPI_cpu) = new int[MPI_Nadj],
 int (*MPI_cube) = new int[MPI_Nadj],

 int (*MPI_cpu_adj) = new int[MPI_Nadj],
 int (*MPI_cube_adj) = new int[MPI_Nadj],

 int (*MPI_direction) = new int[MPI_Nadj],
 int (*MPI_interface) = new int[MPI_Nadj],

 int (*neighbor_cpu_eq) = new int[Ncpu_eq],
 int (*Ncube_Ncpu_eq) = new int[Ncpu_eq], 

 int (*Scube_Ncpu_eq) = new int[Max_nei_eq+1],
 int (*Rcube_Ncpu_eq) = new int[Max_nei_eq+1],

 double (*send_data_curr_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],
 double (*recv_data_curr_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],

 double (*send_data_neig_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],
 double (*recv_data_neig_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1],

 int (*Sdir_eq) = new int[Max_nei_eq+1],
 int (*Rdir_eq) = new int[Max_nei_eq+1],

 int (*ist_eq) = new int[Ncpu_eq],

 int (*csl) = new int[Ncube],

 int (*adj_number)[5][7] = new int[Ncube][5][7],

 int (*adjX_eq) = new int[Ncube],
 int (*adjY_eq) = new int[Ncube],
 int (*adjZ_eq) = new int[Ncube],

 int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size]

 // ================================================================================ //
 );


void BCM_Immersed_boundary
(
// ================================================================================ //
int myid,
int ncube,
int n_wallcube,

int *NBC,

double Xmin, 
double Xmax,
double Ymin, 
double Ymax, 
double Zmin, 
double Zmax, 

double (*cube_size) = new double[Ncube],
int (*csl) = new int[Ncube],

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

int (*wallcube) = new int[N_wallcube]
// ================================================================================ //
);



void BCM_Initial_condition
(
// ================================================================================ //
int myid,
int ncube,
int n_wallcube,
int switch_initial,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);




void BCM_Interface_table
(
// ================================================================================ //
int myid,
int ncube,
int n_wallcube,

int *nadjX_eq, 
int *nadjY_eq, 
int *nadjZ_eq,

int *nadjX_bs_plus, 
int *nadjX_sb_plus,
int *nadjX_bs_minus,
int *nadjX_sb_minus,

int *nadjY_bs_plus, 
int *nadjY_sb_plus, 
int *nadjY_bs_minus,
int *nadjY_sb_minus,

int *nadjZ_bs_plus,
int *nadjZ_sb_plus,
int *nadjZ_bs_minus,
int *nadjZ_sb_minus,

int (*csl) = new int[Ncube],

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
int (*adjZ_sb_minus) = new int[Ncube]
// ================================================================================ //
);




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


int (*csl) = new int[Ncube],

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


double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// ================================================================================ //
);


extern int NBC;

void BCM_Reading_IBM
(
// =================================================== //
int myid,

int *NNBC,

int mp_switch,

double (*weight) = new double[NBC*8+1],

int (*GCindex) = new int[NBC*4+1],

int (*IPsur) = new int[NBC*4+1],

double (*Nor_D) = new double[NBC+1],

double (*Nvec) = new double[NBC*3+1]
// =================================================== //
);


 void BCM_Ghostcell_minus
 (
 // =============================================================================== //
 int myid,

 int *NNBC,

 double Th,

 double (*weight) = new double[NBC_minus*8+1],
 int (*GCindex) = new int[NBC_minus*4+1],
 int (*IPsur) = new int[NBC_minus*4+1],
 double (*Nor_D) = new double[NBC_minus+1],
 double (*Nvec) = new double[NBC_minus*3+1],

 int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

 double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

 // =============================================================================== //
 );


 void BCM_Ghostcell_plus
 (
 // =============================================================================== //
 int myid,

 int *NNBC,

 double Th,

 double (*weight) = new double[NBC_plus*8+1],
 int (*GCindex) = new int[NBC_plus*4+1],
 int (*IPsur) = new int[NBC_plus*4+1],
 double (*Nor_D) = new double[NBC_plus+1],
 double (*Nvec) = new double[NBC_plus*3+1],

 int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

 double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

 // =============================================================================== //
 );



 

 void BCM_Ghostcell_minus_Tem
 (
 // =============================================================================== //
 int myid,

 int *NNBC,

 double Th,

 double (*weight) = new double[NBC_minus*8+1],
 int (*GCindex) = new int[NBC_minus*4+1],
 int (*IPsur) = new int[NBC_minus*4+1],
  double (*Nor_D) = new double[NBC_minus+1],
  double (*Nvec) = new double[NBC_minus*3+1],

 int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

 double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

 // =============================================================================== //
 );


 void BCM_Ghostcell_plus_Tem
 (
 // =============================================================================== //
 int myid,

 int *NNBC,

 double Th,

 double (*weight) = new double[NBC_plus*8+1],
 int (*GCindex) = new int[NBC_plus*4+1],
 int (*IPsur) = new int[NBC_plus*4+1],
 double (*Nor_D) = new double[NBC_plus+1],
 double (*Nvec) = new double[NBC_plus*3+1],

 int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

 double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

 // =============================================================================== //
 );




extern int NXbc_l;
extern int NXbc_u;

void BCM_X_boundary_condition
(
// ================================================================================ //
int myid,

int nXbc_l,
int nXbc_u,

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);



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
);




extern int NZbc_l;
extern int NZbc_u;


void BCM_Z_boundary_condition
(
// ================================================================================ //
int nZbc_l,
int nZbc_u,

int (*Zbc_l) = new int[NZbc_l+1],
int (*Zbc_u) = new int[NZbc_u+1],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);




// ================================================================================================================//
// ========================================= Absorbing boundary conidtion =========================================//
void BCM_Abs_XYZ_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nXbc_l,
int nXbc_u,
int nYbc_l,
int nYbc_u,
int nZbc_l,
int nZbc_u,

double gXmax,
double gXmin,
double gYmax,
double gYmin,
double gZmax,
double gZmin,


double gdXmax,
double gdYmax,
double gdZmax,

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],
int (*Ybc_l) = new int[NYbc_l+1],
int (*Ybc_u) = new int[NYbc_u+1],
int (*Zbc_l) = new int[NYbc_l+1],
int (*Zbc_u) = new int[NYbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*CFL_tau)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]
// ================================================================================ //
);

void BCM_Abs_XY_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nXbc_l,
int nXbc_u,
int nYbc_l,
int nYbc_u,

double gXmax,
double gXmin,
double gYmax,
double gYmin,

double gdXmax,
double gdYmax,

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],
int (*Ybc_l) = new int[NYbc_l+1],
int (*Ybc_u) = new int[NYbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);


void BCM_Abs_X_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nXbc_l,
int nXbc_u,

int (*Xbc_l) = new int[NXbc_l+1],
int (*Xbc_u) = new int[NXbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);



void BCM_Abs_Y_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nYbc_l,
int nYbc_u,

int (*Ybc_l) = new int[NYbc_l+1],
int (*Ybc_u) = new int[NYbc_u+1],

double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);


void BCM_Abs_Z_boundary_condition
(
// ================================================================================ //
int myid,

int ncube,

double deltaT,

double deltaTau,

double e,

int nZbc_l,
int nZbc_u,

int (*Zbc_l) = new int[NZbc_l+1],
int (*Zbc_u) = new int[NZbc_u+1],


double (*cube_size) = new double[Ncube],

double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]
// ================================================================================ //
);



// ========================================= Absorbing boundary conidtion =========================================//
// ================================================================================================================//




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

double (*Fabs)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Roe_dis)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],

double (*Rku1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*Residual1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],


double (*er) = new double[10]

// ================================================================================ //
);



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
);




void BCM_Flux_XYZ_Viscous_LUSGS
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
);



 void BCM_Output
 (
 // ============================================================================ //
 int myid,
int ncube,

int step,

int switch_output,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1q)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*cube_size) = new double[Ncube],


double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size]
 // ============================================================================ //
 );


void BCM_Statistic
(
// =================================================== //
int myid,
int ncube,

int step,
int start_step,
int statistic_step,
int dp_step,

int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*U1_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U2_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U3_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size],
double (*U4_sum)[X_size][Y_size][Z_size] = new double[Ncube][X_size][Y_size][Z_size]
// =================================================== //
);



void BCM_Nusselt_Sphere 
(
// ================================================================================ //
int myid,
int ncube,

double Th,

int (*csl) = new int[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],


double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// ================================================================================ //
);


void BCM_Mean_pressure_coefficient_Sphere 
(
// ================================================================================ //
int myid,
int ncube,

int (*csl) = new int[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],


double (*U1_)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim]

// ================================================================================ //
);


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
);




void BCM_Point_probe
(
// =================================================== //
int myid,
int ncube,
int n_wallcube,

double Xp,
double Yp,
double Zp,

char probe_name[100],

double (*U1)[X_size][Y_size][Z_size][Ndim] = new double[Ncube][X_size][Y_size][Z_size][Ndim],

double (*cube_size) = new double[Ncube],
int (*csl) = new int[Ncube],

double (*Xcube) = new double[Ncube],
double (*Ycube) = new double[Ncube],
double (*Zcube) = new double[Ncube],

double (*Xcnt)[X_size] = new double[Ncube][X_size],
double (*Ycnt)[Y_size] = new double[Ncube][Y_size],
double (*Zcnt)[Z_size] = new double[Ncube][Z_size],

int (*FWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size],

int (*adj_number)[5][7] = new int[Ncube][5][7],

int (*wallcube) = new int[Ncube]
// ================================================== //
);



