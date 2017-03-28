



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;    
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
)

{

#include "MPI_prm.h"
#include "BCM.h"
#include "prm.h"


	int count_index = 0;

	int acube;

	int adjX,adjY,adjZ;
	int ic0,ic1,ic2,ic3,ic4,ic5,ic6;
	int i0,i1,i2,i3,i4,i5,i6;
	int j0,j1,j2,j3,j4,j5,j6;
	int k0,k1,k2,k3,k4,k5,k6;

	int icpu_neig_eq,icpu_neig_sb,icpu_neig_bs;

	int icube_send, icube_recv;

	int iL, iadj;

	int ii,jj,kk,_ii,_jj,_kk,ii_,jj_,kk_;

	int L1, L2, L3, L4, L5;

	int zone = NcubeX*NcubeY*n_buffer;

	
	int index;

	MPI_Request (*Sreq_eq) = new MPI_Request[ncpu_eq];
	MPI_Request (*Rreq_eq) = new MPI_Request[ncpu_eq];

	MPI_Request (*Sreq_sb) = new MPI_Request[ncpu_sb];
	MPI_Request (*Rreq_sb) = new MPI_Request[ncpu_sb];

	MPI_Request (*Sreq_bs) = new MPI_Request[ncpu_bs];
	MPI_Request (*Rreq_bs) = new MPI_Request[ncpu_bs];

	
	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Status istat[8];
	MPI_Request requ1, reqps1;

	
	


	


// ================================================================================================== //
// ============================== # Package for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	
	//start_collection("region1");
	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		
		ii = ii+Ncube_Ncpu_eq[icpu_neig_eq];
		count_index = (ii-Ncube_Ncpu_eq[icpu_neig_eq])*zone;
		
		L1 = 0*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L2 = 1*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L3 = 2*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L4 = 3*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L5 = 4*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;

		
#pragma omp parallel for private(iicube,icube,i,j,k,iL,k2,j2,i2)

		for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_eq[icpu_neig_eq]+icube_send;
			icube = Scube_Ncpu_eq[iicube];
			
			

// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Sdir_eq[iicube] == 6) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							//count_index = count_index+1;
							//iL = iL + 1;
							iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k2 = NcubeZ+k;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							
						}
					}
				}

			}    // ---- if (Sdir_eq[iicube] == 6)  ---- //

// ------------------------------------- [Z+] direction ------------------------------------- //

// ------------------------------------- [Z-] direction ------------------------------------- //

			else if (Sdir_eq[iicube] == 5) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k2 = n_buffer+k;
							
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							

						}
					}
				}

			}    // ---- else if (Sdir_eq[iicube] == 5)  ---- //

// ------------------------------------- [Z-] direction ------------------------------------- //


// ------------------------------------- [Y+] direction ------------------------------------- //

			else if (Sdir_eq[iicube] == 4) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							j2 = NcubeY+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}

			}    // ---- else if (Sdir_eq[iicube] == 4) ---- //


// ------------------------------------- [Y+] direction ------------------------------------- //


// ------------------------------------- [Y-] direction ------------------------------------- //

			else if (Sdir_eq[iicube] == 3) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							j2 = n_buffer+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}

			}    // ---- else if (Sdir_eq[iicube] == 3) ---- //


// ------------------------------------- [Y-] direction ------------------------------------- //


// ------------------------------------- [X+] direction ------------------------------------- //

			else if (Sdir_eq[iicube] == 2) {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
							i2 = NcubeX+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}

			}    // ---- else if (Sdir_eq[iicube] == 2) ---- //


// ------------------------------------- [X+] direction ------------------------------------- //


// ------------------------------------- [X-] direction ------------------------------------- //

			else {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
							i2 = n_buffer+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) ---- //

#pragma omp barrier
		
	}    // ---- for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) ---- //

// ============================== # Package for sending in Current CPU ============================== //
// ================================================================================================== //

	//stop_collection("region1");
	
	
	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		istart = ist_eq[icpu_neig_eq]*5*NcubeX*NcubeY*n_buffer;

		idest = neighbor_cpu_eq[icpu_neig_eq];

		icount = Ncube_Ncpu_eq[icpu_neig_eq]*5*NcubeX*NcubeY*n_buffer;

		itag = 0;

		MPI_Isend((void *)&send_data_curr_eq[istart], icount, MPI_DOUBLE, idest, itag, comm, &Sreq_eq[icpu_neig_eq]);

	}
	
	
	

 //start_collection("region2");

// ================================================================================================== //
// ============================== # Package for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	

	for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) {
	
		ii = ii+Ncube_Ncpu_sb[icpu_neig_sb];
		count_index = (ii-Ncube_Ncpu_sb[icpu_neig_sb])*zone/4;
		
		L1 = 0*Ncube_Ncpu_sb[icpu_neig_sb]*zone/4+5*count_index;
		L2 = 1*Ncube_Ncpu_sb[icpu_neig_sb]*zone/4+5*count_index;
		L3 = 2*Ncube_Ncpu_sb[icpu_neig_sb]*zone/4+5*count_index;
		L4 = 3*Ncube_Ncpu_sb[icpu_neig_sb]*zone/4+5*count_index;
		L5 = 4*Ncube_Ncpu_sb[icpu_neig_sb]*zone/4+5*count_index;


#pragma omp parallel for private(iicube,icube,i,j,k,iL,k2,j2,i2)

		for (icube_send = 1; icube_send <= Ncube_Ncpu_sb[icpu_neig_sb]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_sb[icpu_neig_sb]+icube_send;
			icube = Scube_Ncpu_sb[iicube];
			
			


// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Sdir_sb[iicube] == 6) {

				for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
					for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;
							
							i2 = (i-n_buffer+1)*2;
							j2 = (j-n_buffer+1)*2;
							k2 = nz-2*n_buffer+n_buffer*(k+1)-1;
							
							send_data_curr_sb[L1+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][0]+
														      U1_[icube][i2+1][j2  ][k2  ][0]+
															  U1_[icube][i2  ][j2+1][k2  ][0]+
															  U1_[icube][i2  ][j2  ][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2  ][0]+
															  U1_[icube][i2+1][j2  ][k2+1][0]+
															  U1_[icube][i2  ][j2+1][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2+1][0]);
															  
							send_data_curr_sb[L2+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][1]+
														      U1_[icube][i2+1][j2  ][k2  ][1]+
															  U1_[icube][i2  ][j2+1][k2  ][1]+
															  U1_[icube][i2  ][j2  ][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2  ][1]+
															  U1_[icube][i2+1][j2  ][k2+1][1]+
															  U1_[icube][i2  ][j2+1][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2+1][1]);
															  
							send_data_curr_sb[L3+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][2]+
														      U1_[icube][i2+1][j2  ][k2  ][2]+
															  U1_[icube][i2  ][j2+1][k2  ][2]+
															  U1_[icube][i2  ][j2  ][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2  ][2]+
															  U1_[icube][i2+1][j2  ][k2+1][2]+
															  U1_[icube][i2  ][j2+1][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2+1][2]);
															  
							send_data_curr_sb[L4+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][3]+
														      U1_[icube][i2+1][j2  ][k2  ][3]+
															  U1_[icube][i2  ][j2+1][k2  ][3]+
															  U1_[icube][i2  ][j2  ][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2  ][3]+
															  U1_[icube][i2+1][j2  ][k2+1][3]+
															  U1_[icube][i2  ][j2+1][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2+1][3]);
															  
							send_data_curr_sb[L5+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][4]+
														      U1_[icube][i2+1][j2  ][k2  ][4]+
															  U1_[icube][i2  ][j2+1][k2  ][4]+
															  U1_[icube][i2  ][j2  ][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2  ][4]+
															  U1_[icube][i2+1][j2  ][k2+1][4]+
															  U1_[icube][i2  ][j2+1][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2+1][4]);
							
							//if(myid == 0 & icpu_neig_sb == 2) printf("%f\t",U1_[icube][i][j][k2]);
							
						}
					}
				}

			}    // ---- if (Sdir_sb[iicube] == 6)  ---- //

// ------------------------------------- [Z+] direction ------------------------------------- //

// ------------------------------------- [Z-] direction ------------------------------------- //

			else if (Sdir_sb[iicube] == 5) {

				for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
					for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

							i2 = (i-n_buffer+1)*2;
							j2 = (j-n_buffer+1)*2;
							k2 = n_buffer*(k+1);
							
							send_data_curr_sb[L1+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][0]+
														      U1_[icube][i2+1][j2  ][k2  ][0]+
															  U1_[icube][i2  ][j2+1][k2  ][0]+
															  U1_[icube][i2  ][j2  ][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2  ][0]+
															  U1_[icube][i2+1][j2  ][k2+1][0]+
															  U1_[icube][i2  ][j2+1][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2+1][0]);
															  
							send_data_curr_sb[L2+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][1]+
														      U1_[icube][i2+1][j2  ][k2  ][1]+
															  U1_[icube][i2  ][j2+1][k2  ][1]+
															  U1_[icube][i2  ][j2  ][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2  ][1]+
															  U1_[icube][i2+1][j2  ][k2+1][1]+
															  U1_[icube][i2  ][j2+1][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2+1][1]);
															  
							send_data_curr_sb[L3+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][2]+
														      U1_[icube][i2+1][j2  ][k2  ][2]+
															  U1_[icube][i2  ][j2+1][k2  ][2]+
															  U1_[icube][i2  ][j2  ][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2  ][2]+
															  U1_[icube][i2+1][j2  ][k2+1][2]+
															  U1_[icube][i2  ][j2+1][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2+1][2]);
															  
							send_data_curr_sb[L4+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][3]+
														      U1_[icube][i2+1][j2  ][k2  ][3]+
															  U1_[icube][i2  ][j2+1][k2  ][3]+
															  U1_[icube][i2  ][j2  ][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2  ][3]+
															  U1_[icube][i2+1][j2  ][k2+1][3]+
															  U1_[icube][i2  ][j2+1][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2+1][3]);
															  
							send_data_curr_sb[L5+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][4]+
														      U1_[icube][i2+1][j2  ][k2  ][4]+
															  U1_[icube][i2  ][j2+1][k2  ][4]+
															  U1_[icube][i2  ][j2  ][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2  ][4]+
															  U1_[icube][i2+1][j2  ][k2+1][4]+
															  U1_[icube][i2  ][j2+1][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2+1][4]);
							
						}
					}
				}

			}    // ---- else if (Sdir_sb[iicube] == 5)  ---- //

// ------------------------------------- [Z-] direction ------------------------------------- //


// ------------------------------------- [Y+] direction ------------------------------------- //

			else if (Sdir_sb[iicube] == 4) {

				for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

							iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);
							
							i2 = (i-n_buffer+1)*2;
							j2 = ny-2*n_buffer+n_buffer*(j+1)-1;
							k2 = (k-n_buffer+1)*2;

							send_data_curr_sb[L1+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][0]+
														      U1_[icube][i2+1][j2  ][k2  ][0]+
															  U1_[icube][i2  ][j2+1][k2  ][0]+
															  U1_[icube][i2  ][j2  ][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2  ][0]+
															  U1_[icube][i2+1][j2  ][k2+1][0]+
															  U1_[icube][i2  ][j2+1][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2+1][0]);
															  
							send_data_curr_sb[L2+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][1]+
														      U1_[icube][i2+1][j2  ][k2  ][1]+
															  U1_[icube][i2  ][j2+1][k2  ][1]+
															  U1_[icube][i2  ][j2  ][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2  ][1]+
															  U1_[icube][i2+1][j2  ][k2+1][1]+
															  U1_[icube][i2  ][j2+1][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2+1][1]);
															  
							send_data_curr_sb[L3+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][2]+
														      U1_[icube][i2+1][j2  ][k2  ][2]+
															  U1_[icube][i2  ][j2+1][k2  ][2]+
															  U1_[icube][i2  ][j2  ][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2  ][2]+
															  U1_[icube][i2+1][j2  ][k2+1][2]+
															  U1_[icube][i2  ][j2+1][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2+1][2]);
															  
							send_data_curr_sb[L4+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][3]+
														      U1_[icube][i2+1][j2  ][k2  ][3]+
															  U1_[icube][i2  ][j2+1][k2  ][3]+
															  U1_[icube][i2  ][j2  ][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2  ][3]+
															  U1_[icube][i2+1][j2  ][k2+1][3]+
															  U1_[icube][i2  ][j2+1][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2+1][3]);
															  
							send_data_curr_sb[L5+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][4]+
														      U1_[icube][i2+1][j2  ][k2  ][4]+
															  U1_[icube][i2  ][j2+1][k2  ][4]+
															  U1_[icube][i2  ][j2  ][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2  ][4]+
															  U1_[icube][i2+1][j2  ][k2+1][4]+
															  U1_[icube][i2  ][j2+1][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2+1][4]);
							
						}
					}
				}

			}    // ---- else if (Sdir_sb[iicube] == 4) ---- //


// ------------------------------------- [Y+] direction ------------------------------------- //


// ------------------------------------- [Y-] direction ------------------------------------- //

			else if (Sdir_sb[iicube] == 3) {

				for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  
							
							iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);
							
							i2 = (i-n_buffer+1)*2;
							j2 = n_buffer*(j+1);
							k2 = (k-n_buffer+1)*2;

							
							send_data_curr_sb[L1+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][0]+
														      U1_[icube][i2+1][j2  ][k2  ][0]+
															  U1_[icube][i2  ][j2+1][k2  ][0]+
															  U1_[icube][i2  ][j2  ][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2  ][0]+
															  U1_[icube][i2+1][j2  ][k2+1][0]+
															  U1_[icube][i2  ][j2+1][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2+1][0]);
															  
							send_data_curr_sb[L2+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][1]+
														      U1_[icube][i2+1][j2  ][k2  ][1]+
															  U1_[icube][i2  ][j2+1][k2  ][1]+
															  U1_[icube][i2  ][j2  ][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2  ][1]+
															  U1_[icube][i2+1][j2  ][k2+1][1]+
															  U1_[icube][i2  ][j2+1][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2+1][1]);
															  
							send_data_curr_sb[L3+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][2]+
														      U1_[icube][i2+1][j2  ][k2  ][2]+
															  U1_[icube][i2  ][j2+1][k2  ][2]+
															  U1_[icube][i2  ][j2  ][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2  ][2]+
															  U1_[icube][i2+1][j2  ][k2+1][2]+
															  U1_[icube][i2  ][j2+1][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2+1][2]);
															  
							send_data_curr_sb[L4+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][3]+
														      U1_[icube][i2+1][j2  ][k2  ][3]+
															  U1_[icube][i2  ][j2+1][k2  ][3]+
															  U1_[icube][i2  ][j2  ][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2  ][3]+
															  U1_[icube][i2+1][j2  ][k2+1][3]+
															  U1_[icube][i2  ][j2+1][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2+1][3]);
															  
							send_data_curr_sb[L5+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][4]+
														      U1_[icube][i2+1][j2  ][k2  ][4]+
															  U1_[icube][i2  ][j2+1][k2  ][4]+
															  U1_[icube][i2  ][j2  ][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2  ][4]+
															  U1_[icube][i2+1][j2  ][k2+1][4]+
															  U1_[icube][i2  ][j2+1][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2+1][4]);
							
						}
					}
				}

			}    // ---- else if (Sdir_sb[iicube] == 3) ---- //


// ------------------------------------- [Y-] direction ------------------------------------- //


// ------------------------------------- [X+] direction ------------------------------------- //

			else if (Sdir_sb[iicube] == 2) {

				
				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
						for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

							iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);
							
							i2 = nx-2*n_buffer+n_buffer*(i+1)-1;
							j2 = (j-n_buffer+1)*2;
							k2 = (k-n_buffer+1)*2;

							
							send_data_curr_sb[L1+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][0]+
														      U1_[icube][i2+1][j2  ][k2  ][0]+
															  U1_[icube][i2  ][j2+1][k2  ][0]+
															  U1_[icube][i2  ][j2  ][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2  ][0]+
															  U1_[icube][i2+1][j2  ][k2+1][0]+
															  U1_[icube][i2  ][j2+1][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2+1][0]);
															  
							send_data_curr_sb[L2+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][1]+
														      U1_[icube][i2+1][j2  ][k2  ][1]+
															  U1_[icube][i2  ][j2+1][k2  ][1]+
															  U1_[icube][i2  ][j2  ][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2  ][1]+
															  U1_[icube][i2+1][j2  ][k2+1][1]+
															  U1_[icube][i2  ][j2+1][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2+1][1]);
															  
							send_data_curr_sb[L3+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][2]+
														      U1_[icube][i2+1][j2  ][k2  ][2]+
															  U1_[icube][i2  ][j2+1][k2  ][2]+
															  U1_[icube][i2  ][j2  ][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2  ][2]+
															  U1_[icube][i2+1][j2  ][k2+1][2]+
															  U1_[icube][i2  ][j2+1][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2+1][2]);
															  
							send_data_curr_sb[L4+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][3]+
														      U1_[icube][i2+1][j2  ][k2  ][3]+
															  U1_[icube][i2  ][j2+1][k2  ][3]+
															  U1_[icube][i2  ][j2  ][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2  ][3]+
															  U1_[icube][i2+1][j2  ][k2+1][3]+
															  U1_[icube][i2  ][j2+1][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2+1][3]);
															  
							send_data_curr_sb[L5+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][4]+
														      U1_[icube][i2+1][j2  ][k2  ][4]+
															  U1_[icube][i2  ][j2+1][k2  ][4]+
															  U1_[icube][i2  ][j2  ][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2  ][4]+
															  U1_[icube][i2+1][j2  ][k2+1][4]+
															  U1_[icube][i2  ][j2+1][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2+1][4]);
							
						}
					}
				}

			}    // ---- else if (Sdir_sb[iicube] == 2) ---- //


// ------------------------------------- [X+] direction ------------------------------------- //


// ------------------------------------- [X-] direction ------------------------------------- //

			else {

				
				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
						for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  
							
							iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);
							
							i2 = n_buffer*(i+1);
							j2 = (j-n_buffer+1)*2;
							k2 = (k-n_buffer+1)*2;

							send_data_curr_sb[L1+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][0]+
														      U1_[icube][i2+1][j2  ][k2  ][0]+
															  U1_[icube][i2  ][j2+1][k2  ][0]+
															  U1_[icube][i2  ][j2  ][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2  ][0]+
															  U1_[icube][i2+1][j2  ][k2+1][0]+
															  U1_[icube][i2  ][j2+1][k2+1][0]+
															  U1_[icube][i2+1][j2+1][k2+1][0]);
															  
							send_data_curr_sb[L2+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][1]+
														      U1_[icube][i2+1][j2  ][k2  ][1]+
															  U1_[icube][i2  ][j2+1][k2  ][1]+
															  U1_[icube][i2  ][j2  ][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2  ][1]+
															  U1_[icube][i2+1][j2  ][k2+1][1]+
															  U1_[icube][i2  ][j2+1][k2+1][1]+
															  U1_[icube][i2+1][j2+1][k2+1][1]);
															  
							send_data_curr_sb[L3+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][2]+
														      U1_[icube][i2+1][j2  ][k2  ][2]+
															  U1_[icube][i2  ][j2+1][k2  ][2]+
															  U1_[icube][i2  ][j2  ][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2  ][2]+
															  U1_[icube][i2+1][j2  ][k2+1][2]+
															  U1_[icube][i2  ][j2+1][k2+1][2]+
															  U1_[icube][i2+1][j2+1][k2+1][2]);
															  
							send_data_curr_sb[L4+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][3]+
														      U1_[icube][i2+1][j2  ][k2  ][3]+
															  U1_[icube][i2  ][j2+1][k2  ][3]+
															  U1_[icube][i2  ][j2  ][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2  ][3]+
															  U1_[icube][i2+1][j2  ][k2+1][3]+
															  U1_[icube][i2  ][j2+1][k2+1][3]+
															  U1_[icube][i2+1][j2+1][k2+1][3]);
															  
							send_data_curr_sb[L5+iL] = 0.125*(U1_[icube][i2  ][j2  ][k2  ][4]+
														      U1_[icube][i2+1][j2  ][k2  ][4]+
															  U1_[icube][i2  ][j2+1][k2  ][4]+
															  U1_[icube][i2  ][j2  ][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2  ][4]+
															  U1_[icube][i2+1][j2  ][k2+1][4]+
															  U1_[icube][i2  ][j2+1][k2+1][4]+
															  U1_[icube][i2+1][j2+1][k2+1][4]);
							
							
						}
					}
				}

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_sb[icpu_neig_sb]; icube_send++) ---- //

#pragma omp barrier
		
	}    // ---- for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) ---- //

// ============================== # Package for sending in Current CPU ============================== //
// ================================================================================================== //

//stop_collection("region2");


	
//start_collection("region3");

// ================================================================================================== //
// ============================== # Package for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	

	for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) {

	
		ii = ii+Ncube_Ncpu_bs[icpu_neig_bs];
		count_index = (ii-Ncube_Ncpu_bs[icpu_neig_bs])*zone;
		
		L1 = 0*Ncube_Ncpu_bs[icpu_neig_bs]*zone+5*count_index;
		L2 = 1*Ncube_Ncpu_bs[icpu_neig_bs]*zone+5*count_index;
		L3 = 2*Ncube_Ncpu_bs[icpu_neig_bs]*zone+5*count_index;
		L4 = 3*Ncube_Ncpu_bs[icpu_neig_bs]*zone+5*count_index;
		L5 = 4*Ncube_Ncpu_bs[icpu_neig_bs]*zone+5*count_index;


#pragma omp parallel for private(iicube,icube,i,j,k,iL,k2,j2,i2)

		for (icube_send = 1; icube_send <= Ncube_Ncpu_bs[icpu_neig_bs]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_bs[icpu_neig_bs]+icube_send;
			icube = Scube_Ncpu_bs[iicube];

// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Sdir_bs[iicube] == 6) {

// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 1) {

					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1;
								j2 = j/n_buffer+1;
								k2 = nz+k/n_buffer;

								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 2) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = j/n_buffer+1;
								k2 = nz+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 3) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = nz+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 4) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = nz+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//
				
			}    // ---- if (Sdir_bs[iicube] == 6)  ---- //

// ------------------------------------- [Z+] direction ------------------------------------- //


// ------------------------------------- [Z-] direction ------------------------------------- //

			else if (Sdir_bs[iicube] == 5) {

				
// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 1) {

					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1;
								j2 = j/n_buffer+1;
								k2 = n_buffer+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 2) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
								
								i2 = i/n_buffer+1+NcubeX/2;
								j2 = j/n_buffer+1;
								k2 = n_buffer+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 3) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = n_buffer+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 4) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = n_buffer+k/n_buffer;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- else if (Sdir_bs[iicube] == 5)  ---- //

// ------------------------------------- [Z-] direction ------------------------------------- //


// ------------------------------------- [Y+] direction ------------------------------------- //

			else if (Sdir_bs[iicube] == 4) {

					
// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 1) {

					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1;
								j2 = ny+j/n_buffer;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 2) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = ny+j/n_buffer;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 3) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1;
								j2 = ny+j/n_buffer;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 4) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = ny+j/n_buffer;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


			}    // ---- else if (Sdir_bs[iicube] == 4) ---- //


// ------------------------------------- [Y+] direction ------------------------------------- //


// ------------------------------------- [Y-] direction ------------------------------------- //

			else if (Sdir_bs[iicube] == 3) {

				
// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 1) {

					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1;
								j2 = n_buffer+j/n_buffer;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 2) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = n_buffer+j/n_buffer;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 3) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1;
								j2 = n_buffer+j/n_buffer;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 4) {
				
					for (i = n_buffer; i <= nx; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);

								i2 = i/n_buffer+1+NcubeX/2;
								j2 = n_buffer+j/n_buffer;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];

							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- else if (Sdir_bs[iicube] == 3) ---- //


// ------------------------------------- [Y-] direction ------------------------------------- //


// ------------------------------------- [X+] direction ------------------------------------- //

			else if (Sdir_bs[iicube] == 2) {

					
// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 1) {

					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);

								i2 = nx+i/n_buffer;
								j2 = j/n_buffer+1;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 2) {
				
					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
								
								i2 = nx+i/n_buffer;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 3) {
				
					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);

								i2 = nx+i/n_buffer;
								j2 = j/n_buffer+1;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 4) {
				
					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);

								i2 = nx+i/n_buffer;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- else if (Sdir_bs[iicube] == 2) ---- //


// ------------------------------------- [X+] direction ------------------------------------- //


// ------------------------------------- [X-] direction ------------------------------------- //

			else {

								
// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 1) {

					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);

								i2 = n_buffer+i/n_buffer;
								j2 = j/n_buffer+1;
								k2 = k/n_buffer+1;

								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 2) {
				
					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);

								i2 = n_buffer+i/n_buffer;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = k/n_buffer+1;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 3) {
				
					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);

								i2 = n_buffer+i/n_buffer;
								j2 = j/n_buffer+1;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (SadjN_bs[iicube] == 4) {
				
					for (i = 0; i < n_buffer; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
								
								i2 = n_buffer+i/n_buffer;
								j2 = j/n_buffer+1+NcubeY/2;
								k2 = k/n_buffer+1+NcubeZ/2;
								
								send_data_curr_bs[L1+iL] = U1_[icube][i2][j2][k2][0];
								send_data_curr_bs[L2+iL] = U1_[icube][i2][j2][k2][1];
								send_data_curr_bs[L3+iL] = U1_[icube][i2][j2][k2][2];
								send_data_curr_bs[L4+iL] = U1_[icube][i2][j2][k2][3];
								send_data_curr_bs[L5+iL] = U1_[icube][i2][j2][k2][4];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_bs[icpu_neig_bs]; icube_send++) ---- //
		
#pragma omp barrier

	}    // ---- for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) ---- //

// ============================== # Package for sending in Current CPU ============================== //
// ================================================================================================== //

//stop_collection("region3");

	
	/*
	MPI_Datatype Type_recv_eq, Type_recv_bs, Type_recv_sb, Type_send_eq, Type_send_bs, Type_send_sb;

	MPI_Type_vector(1, 5*NcubeX*NcubeY*n_buffer*max_nei_bs+1, 5*NcubeX*NcubeY*n_buffer*max_nei_bs+1, MPI_DOUBLE, &Type_recv_eq);


	MPI_Type_vector(1, 5*NcubeX*NcubeY*n_buffer*max_nei_bs+1, 5*NcubeX*NcubeY*n_buffer*max_nei_bs+1, MPI_DOUBLE, &Type_recv_bs);
	MPI_Type_vector(1, 5*NcubeX*NcubeY*n_buffer*max_nei_sb+1, 5*NcubeX*NcubeY*n_buffer*max_nei_sb+1, MPI_DOUBLE, &Type_recv_sb);

	MPI_Type_vector(1, 5*NcubeX*NcubeY*n_buffer*max_nei_eq+1, 5*NcubeX*NcubeY*n_buffer*max_nei_eq+1, MPI_DOUBLE, &Type_send_eq);
	MPI_Type_vector(1, 5*NcubeX*NcubeY*n_buffer*max_nei_bs+1, 5*NcubeX*NcubeY*n_buffer*max_nei_bs+1, MPI_DOUBLE, &Type_send_bs);
	MPI_Type_vector(1, 5*NcubeX*NcubeY*n_buffer*max_nei_sb+1, 5*NcubeX*NcubeY*n_buffer*max_nei_sb+1, MPI_DOUBLE, &Type_send_sb);

	MPI_Type_commit(&Type_recv_eq);
	MPI_Type_commit(&Type_recv_bs);
	MPI_Type_commit(&Type_recv_sb);

	MPI_Type_commit(&Type_send_eq);
	MPI_Type_commit(&Type_send_bs);
	MPI_Type_commit(&Type_send_sb);
	*/


	
	for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) {

		istart = ist_sb[icpu_neig_sb]*5*NcubeX*NcubeY*n_buffer/4;

		idest = neighbor_cpu_sb[icpu_neig_sb];

		icount = Ncube_Ncpu_sb[icpu_neig_sb]*5*NcubeX*NcubeY*n_buffer/4;

		itag = 0;

		MPI_Isend((void *)&send_data_curr_sb[istart], icount, MPI_DOUBLE, idest, itag, comm, &Sreq_sb[icpu_neig_sb]);

	}
	
	
	for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) {

		istart = ist_bs[icpu_neig_bs]*5*NcubeX*NcubeY*n_buffer;

		idest = neighbor_cpu_bs[icpu_neig_bs];

		icount = Ncube_Ncpu_bs[icpu_neig_bs]*5*NcubeX*NcubeY*n_buffer;

		itag = 0;

		MPI_Isend((void *)&send_data_curr_bs[istart], icount, MPI_DOUBLE, idest, itag, comm, &Sreq_bs[icpu_neig_bs]);

	}
	
	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		istart = ist_eq[icpu_neig_eq]*5*NcubeX*NcubeY*n_buffer;

		isrc = neighbor_cpu_eq[icpu_neig_eq];

		icount = Ncube_Ncpu_eq[icpu_neig_eq]*5*NcubeX*NcubeY*n_buffer;

		itag = 0;

		MPI_Irecv((void *)&recv_data_curr_eq[istart], icount, MPI_DOUBLE, isrc, itag, comm, &Rreq_eq[icpu_neig_eq]);

	}
	
	
	for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) {

		istart = ist_bs[icpu_neig_bs]*5*NcubeX*NcubeY*n_buffer/4;

		isrc = neighbor_cpu_bs[icpu_neig_bs];

		icount = Ncube_Ncpu_bs[icpu_neig_bs]*5*NcubeX*NcubeY*n_buffer/4;

		itag = 0;

		MPI_Irecv((void *)&recv_data_curr_bs[istart], icount, MPI_DOUBLE, isrc, itag, comm, &Rreq_bs[icpu_neig_bs]);

	}
	
	
	for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) {

		istart = ist_sb[icpu_neig_sb]*5*NcubeX*NcubeY*n_buffer;

		isrc = neighbor_cpu_sb[icpu_neig_sb];

		icount = Ncube_Ncpu_sb[icpu_neig_sb]*5*NcubeX*NcubeY*n_buffer;

		itag = 0;

		MPI_Irecv((void *)&recv_data_curr_sb[istart], icount, MPI_DOUBLE, isrc, itag, comm, &Rreq_sb[icpu_neig_sb]);

	}
	
	
	
	//start_collection("region4");

// ---- the same size connection in X direction ---- //

#pragma omp parallel for private(ic0, ic2, i,j,k, i0,i2) schedule(dynamic)
// ============================================== //
	for (adjX = 1; adjX <= nadjX_eq; adjX++) {    //  
// ============================================== //

		ic0 = adjX_eq[adjX];
		ic2 = adj_number[ic0][1][2];
		
		for (i = 0; i < n_buffer; i++) {
//// #pragma omp parallel for private(k,i0,i2)
			for (j = n_buffer; j <= ny; j++) {
				for (k = n_buffer; k <= nz; k++) {  

					i0 = NcubeX+i;
					i2 = i;
					
					U1_[ic2][i2][j][k][0] = U1_[ic0][i0][j][k][0];
					U1_[ic2][i2][j][k][1] = U1_[ic0][i0][j][k][1];
					U1_[ic2][i2][j][k][2] = U1_[ic0][i0][j][k][2];
					U1_[ic2][i2][j][k][3] = U1_[ic0][i0][j][k][3];
					U1_[ic2][i2][j][k][4] = U1_[ic0][i0][j][k][4];
					
					i0 = n_buffer+NcubeX+i;
					i2 = n_buffer+i;

					U1_[ic0][i0][j][k][0] = U1_[ic2][i2][j][k][0];
					U1_[ic0][i0][j][k][1] = U1_[ic2][i2][j][k][1];
					U1_[ic0][i0][j][k][2] = U1_[ic2][i2][j][k][2];
					U1_[ic0][i0][j][k][3] = U1_[ic2][i2][j][k][3];
					U1_[ic0][i0][j][k][4] = U1_[ic2][i2][j][k][4];
					
				}
			}
		}
// ============================================== //
	}                                             //
// ============================================== //





	
// ---- big cube to small cube in [X+] direction for small cube ---- //


#pragma omp parallel for private(ic0, ic2, i,j,k, i0,j0,k0) schedule(dynamic)
// ============================================== //
	for (adjX = 1; adjX <= nadjX_bs_plus; adjX++) {    //  
// ============================================== //
		
		ic0 = adjX_bs_plus[adjX];
		ic2 = adj_number[ic0][1][2];


		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i/n_buffer;
						j0 = j/n_buffer+1;
						k0 = k/n_buffer+1;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];

					}
				}
			}

		}


		ic2 = adj_number[ic0][2][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i/n_buffer;
						j0 = j/n_buffer+1+NcubeY/2;
						k0 = k/n_buffer+1;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];
					}
				}
			}

		}


		ic2 = adj_number[ic0][3][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i/n_buffer;
						j0 = j/n_buffer+1;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];
					}
				}
			}

		}


		ic2 = adj_number[ic0][4][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i/n_buffer;
						j0 = j/n_buffer+1+NcubeY/2;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];
					}
				}
			}

		}


// ============================================== //
	}                                             //
// ============================================== //
	


// ---- small cube to big cube in [X+] direction for small cube ---- //


#pragma omp parallel for private(ic0, ic2,i,j,k,i0,i2,j0,k0) schedule(dynamic)
// ============================================== //
	for (adjX = 1; adjX <= nadjX_bs_minus; adjX++) {    //  
// ============================================== //

		ic2 = adjX_bs_minus[adjX];

		ic0 = adj_number[ic2][1][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
// #pragma omp parallel for private(k,i0,i2,j0,k0)
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i+1;
						i2 = n_buffer+i/n_buffer;
						j0 = j/n_buffer+1;
						k0 = k/n_buffer+1;

						U1_[ic0][i0][j][k][0] = U1_[ic2][i2][j0][k0][0];
						U1_[ic0][i0][j][k][1] = U1_[ic2][i2][j0][k0][1];
						U1_[ic0][i0][j][k][2] = U1_[ic2][i2][j0][k0][2];
						U1_[ic0][i0][j][k][3] = U1_[ic2][i2][j0][k0][3];
						U1_[ic0][i0][j][k][4] = U1_[ic2][i2][j0][k0][4];

					}
				}
			}

		}


		ic0 = adj_number[ic2][2][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
// #pragma omp parallel for private(k,i0,i2,j0,k0)
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i+1;
						i2 = n_buffer+i/n_buffer;
						j0 = j/n_buffer+1+NcubeY/2;
						k0 = k/n_buffer+1;

						U1_[ic0][i0][j][k][0] = U1_[ic2][i2][j0][k0][0];
						U1_[ic0][i0][j][k][1] = U1_[ic2][i2][j0][k0][1];
						U1_[ic0][i0][j][k][2] = U1_[ic2][i2][j0][k0][2];
						U1_[ic0][i0][j][k][3] = U1_[ic2][i2][j0][k0][3];
						U1_[ic0][i0][j][k][4] = U1_[ic2][i2][j0][k0][4];

					}
				}
			}

		}


		ic0 = adj_number[ic2][3][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
// #pragma omp parallel for private(k,i0,i2,j0,k0)
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i+1;
						i2 = n_buffer+i/n_buffer;
						j0 = j/n_buffer+1;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic0][i0][j][k][0] = U1_[ic2][i2][j0][k0][0];
						U1_[ic0][i0][j][k][1] = U1_[ic2][i2][j0][k0][1];
						U1_[ic0][i0][j][k][2] = U1_[ic2][i2][j0][k0][2];
						U1_[ic0][i0][j][k][3] = U1_[ic2][i2][j0][k0][3];
						U1_[ic0][i0][j][k][4] = U1_[ic2][i2][j0][k0][4];

					}
				}
			}

		}


		ic0 = adj_number[ic2][4][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
// #pragma omp parallel for private(k,i0,i2,j0,k0)
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = nx+i+1;
						i2 = n_buffer+i/n_buffer;
						j0 = j/n_buffer+1+NcubeY/2;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic0][i0][j][k][0] = U1_[ic2][i2][j0][k0][0];
						U1_[ic0][i0][j][k][1] = U1_[ic2][i2][j0][k0][1];
						U1_[ic0][i0][j][k][2] = U1_[ic2][i2][j0][k0][2];
						U1_[ic0][i0][j][k][3] = U1_[ic2][i2][j0][k0][3];
						U1_[ic0][i0][j][k][4] = U1_[ic2][i2][j0][k0][4];

					}
				}
			}

		}

// ============================================== //
	}                                             //
// ============================================== //



// ---- big cube to small cube in [X+] direction for big cube ---- //

#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,i2,j2,k2) schedule(dynamic)
// ============================================== //
	for (adjX = 1; adjX <= nadjX_bs_plus; adjX++) {    //  
// ============================================== //
		
		ic0 = adjX_bs_plus[adjX];

		ic2 = adj_number[ic0][1][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						//ic2 = adj_number[ic0][1][2];

						i0 = n_buffer+NcubeX+i;
						j0 = j;
						k0 = k;

						i2 = n_buffer*(i+1);
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);
					}
				}
			}

		}



		ic2 = adj_number[ic0][2][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = n_buffer+NcubeX+i;
						j0 = j+NcubeY/2;
						k0 = k;

						i2 = n_buffer*(i+1);
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);
					}
				}
			}

		}


		ic2 = adj_number[ic0][3][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = n_buffer+NcubeX+i;
						j0 = j;
						k0 = k+NcubeZ/2;;

						i2 = n_buffer*(i+1);
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);
					}
				}
			}

		}


		ic2 = adj_number[ic0][4][2];

		if (ic2 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = n_buffer+NcubeX+i;
						j0 = j+NcubeY/2;
						k0 = k+NcubeZ/2;

						i2 = n_buffer*(i+1);
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


// ============================================== //
	}                                             //
// ============================================== //




// ---- small cube to big cube in [X+] direction for big cube ---- //

#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,i2,j2,k2) schedule(dynamic)
// ============================================== //
	for (adjX = 1; adjX <= nadjX_bs_minus; adjX++) {    //  
// ============================================== //

		ic2 = adjX_bs_minus[adjX];

		ic0 = adj_number[ic2][1][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i;
						j0 = j;
						k0 = k;

						i2 = nx-2*n_buffer+n_buffer*(i+1)-1;
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic0 = adj_number[ic2][2][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  


						i0 = i;
						j0 = j+NcubeY/2;
						k0 = k;

						i2 = nx-2*n_buffer+n_buffer*(i+1)-1;
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic0 = adj_number[ic2][3][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  


						i0 = i;
						j0 = j;
						k0 = k+NcubeZ/2;

						i2 = nx-2*n_buffer+n_buffer*(i+1)-1;
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic0 = adj_number[ic2][4][1];

		if (ic0 > 0) {

			for (i = 0; i < n_buffer; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  


						i0 = i;
						j0 = j+NcubeY/2;
						k0 = k+NcubeZ/2;

						i2 = nx-2*n_buffer+n_buffer*(i+1)-1;
						j2 = (j-n_buffer+1)*2;
						k2 = (k-n_buffer+1)*2;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}

// ============================================== //
	}                                             //
// ============================================== //




#pragma omp parallel for private(ic0, ic2, i,j,k, j0,j2) schedule(dynamic)

// ---- the same size connection in Y direction ---- //
// ============================================== //
	for (adjY = 1; adjY <= nadjY_eq; adjY++) {    //  
// ============================================== //

		ic0 = adjY_eq[adjY];
		ic2 = adj_number[ic0][1][4];

//// #pragma omp parallel for private(j,k,j0,j2)
		for (i = n_buffer; i <= nx; i++) {
			for (j = 0; j < n_buffer; j++) {
				for (k = n_buffer; k <= nz; k++) {  

					j0 = NcubeY+j;
					j2 = j;

					U1_[ic2][i][j2][k][0] = U1_[ic0][i][j0][k][0];
					U1_[ic2][i][j2][k][1] = U1_[ic0][i][j0][k][1];
					U1_[ic2][i][j2][k][2] = U1_[ic0][i][j0][k][2];
					U1_[ic2][i][j2][k][3] = U1_[ic0][i][j0][k][3];
					U1_[ic2][i][j2][k][4] = U1_[ic0][i][j0][k][4];
					
					j0 = n_buffer+NcubeY+j;
					j2 = n_buffer+j;

					U1_[ic0][i][j0][k][0] = U1_[ic2][i][j2][k][0];
					U1_[ic0][i][j0][k][1] = U1_[ic2][i][j2][k][1];
					U1_[ic0][i][j0][k][2] = U1_[ic2][i][j2][k][2];
					U1_[ic0][i][j0][k][3] = U1_[ic2][i][j2][k][3];
					U1_[ic0][i][j0][k][4] = U1_[ic2][i][j2][k][4];
					
				}
			}
		}
// ============================================== //
	}                                             //
// ============================================== //








// ---- big cube to small cube in [Y+] direction for small cube ---- //



#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0) schedule(dynamic)

// ============================================== //
	for (adjY = 1; adjY <= nadjY_bs_plus; adjY++) {    //  
// ============================================== //

		ic0 = adjY_bs_plus[adjY];

		ic2 = adj_number[ic0][1][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1;
						j0 = ny+j/n_buffer;
						k0 = k/n_buffer+1;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];

					}
				}
			}

		}

		ic2 = adj_number[ic0][2][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1+NcubeX/2;
						j0 = ny+j/n_buffer;
						k0 = k/n_buffer+1;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];


					}
				}
			}

		}

		ic2 = adj_number[ic0][3][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1;
						j0 = ny+j/n_buffer;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];


					}
				}
			}

		}

		ic2 = adj_number[ic0][4][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1+NcubeX/2;
						j0 = ny+j/n_buffer;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];

					}
				}
			}

		}
	
// ============================================== //
	}                                             //
// ============================================== //

// ---- small cube to big cube in [Y+] direction for small cube ---- //


#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,j2,k0) schedule(dynamic)

// ============================================== //
	for (adjY = 1; adjY <= nadjY_bs_minus; adjY++) {    //  
// ============================================== //

		ic2 =adjY_bs_minus[adjY];

		ic0 = adj_number[ic2][1][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1;
						j0 = j+ny+1;
						j2 = n_buffer+j/n_buffer;
						k0 = k/n_buffer+1;

						U1_[ic0][i][j0][k][0] = U1_[ic2][i0][j2][k0][0];
						U1_[ic0][i][j0][k][1] = U1_[ic2][i0][j2][k0][1];
						U1_[ic0][i][j0][k][2] = U1_[ic2][i0][j2][k0][2];
						U1_[ic0][i][j0][k][3] = U1_[ic2][i0][j2][k0][3];
						U1_[ic0][i][j0][k][4] = U1_[ic2][i0][j2][k0][4];


					}
				}
			}

		}


		ic0 = adj_number[ic2][2][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1+NcubeX/2;
						j0 = j+ny+1;
						j2 = n_buffer+j/n_buffer;
						k0 = k/n_buffer+1;

						U1_[ic0][i][j0][k][0] = U1_[ic2][i0][j2][k0][0];
						U1_[ic0][i][j0][k][1] = U1_[ic2][i0][j2][k0][1];
						U1_[ic0][i][j0][k][2] = U1_[ic2][i0][j2][k0][2];
						U1_[ic0][i][j0][k][3] = U1_[ic2][i0][j2][k0][3];
						U1_[ic0][i][j0][k][4] = U1_[ic2][i0][j2][k0][4];


					}
				}
			}

		}


		ic0 = adj_number[ic2][3][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1;
						j0 = j+ny+1;
						j2 = n_buffer+j/n_buffer;
						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic0][i][j0][k][0] = U1_[ic2][i0][j2][k0][0];
						U1_[ic0][i][j0][k][1] = U1_[ic2][i0][j2][k0][1];
						U1_[ic0][i][j0][k][2] = U1_[ic2][i0][j2][k0][2];
						U1_[ic0][i][j0][k][3] = U1_[ic2][i0][j2][k0][3];
						U1_[ic0][i][j0][k][4] = U1_[ic2][i0][j2][k0][4];

					}
				}
			}

		}


		ic0 = adj_number[ic2][4][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						i0 = i/n_buffer+1+NcubeX/2;

						j0 = j+ny+1;
						j2 = n_buffer+j/n_buffer;

						k0 = k/n_buffer+1+NcubeZ/2;

						U1_[ic0][i][j0][k][0] = U1_[ic2][i0][j2][k0][0];
						U1_[ic0][i][j0][k][1] = U1_[ic2][i0][j2][k0][1];
						U1_[ic0][i][j0][k][2] = U1_[ic2][i0][j2][k0][2];
						U1_[ic0][i][j0][k][3] = U1_[ic2][i0][j2][k0][3];
						U1_[ic0][i][j0][k][4] = U1_[ic2][i0][j2][k0][4];


					}
				}
			}

		}

// ============================================== //
	}                                             //
// ============================================== //






// ---- big cube to small cube in [Y+] direction for big cube ---- //

#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,i2,j2,k2) schedule(dynamic)

// ============================================== //
	for (adjY = 1; adjY <= nadjY_bs_plus; adjY++) {    //  
// ============================================== //
		
		ic0 = adjY_bs_plus[adjY];

		ic2 = adj_number[ic0][1][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i;
						j0 = n_buffer+NcubeY+j;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = n_buffer*(j+1);
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic2 = adj_number[ic0][2][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i+NcubeX/2;
						j0 = n_buffer+NcubeY+j;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = n_buffer*(j+1);
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic2 = adj_number[ic0][3][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i;
						j0 = n_buffer+NcubeY+j;
						k0 = k+NcubeZ/2;

						i2 = (i-n_buffer+1)*2;
						j2 = n_buffer*(j+1);
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic2 = adj_number[ic0][4][4];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i+NcubeX/2;
						j0 = n_buffer+NcubeY+j;
						k0 = k+NcubeZ/2;

						i2 = (i-n_buffer+1)*2;
						j2 = n_buffer*(j+1);
						k2 = (k-n_buffer+1)*2;

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);
					}
				}
			}

		}

					
// ============================================== //
	}                                             //
// ============================================== //





// ---- small cube to big cube in [Y+] direction for big cube ---- //

#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,i2,j2,k2) schedule(dynamic)

// ============================================== //
	for (adjY = 1; adjY <= nadjY_bs_minus; adjY++) {    //  
// ============================================== //

		ic2 = adjY_bs_minus[adjY];

		ic0 = adj_number[ic2][1][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i;
						j0 = j;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = ny-2*n_buffer+n_buffer*(j+1)-1;
						k2 = (k-n_buffer+1)*2;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic0 = adj_number[ic2][2][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i+NcubeX/2;
						j0 = j;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = ny-2*n_buffer+n_buffer*(j+1)-1;
						k2 = (k-n_buffer+1)*2;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);


					}
				}
			}

		}



		ic0 = adj_number[ic2][3][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i;
						j0 = j;
						k0 = k+NcubeZ/2;

						i2 = (i-n_buffer+1)*2;
						j2 = ny-2*n_buffer+n_buffer*(j+1)-1;
						k2 = (k-n_buffer+1)*2;

						
						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic0 = adj_number[ic2][4][3];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = 0; j < n_buffer; j++) {
					for (k = n_buffer; k <= (nz-n_buffer+1)/2+1; k++) {  

						i0 = i+NcubeX/2;
						j0 = j;
						k0 = k+NcubeZ/2;

						i2 = (i-n_buffer+1)*2;
						j2 = ny-2*n_buffer+n_buffer*(j+1)-1;
						k2 = (k-n_buffer+1)*2;

						
						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);


					}
				}
			}

		}

// ============================================== //
	}                                             //
// ============================================== //



#pragma omp parallel for private(ic0, ic2, i,j,k, k0,k2) schedule(dynamic)

// ---- the same size connection in Z direction ---- //
// ============================================== //
	for (adjZ = 1; adjZ <= nadjZ_eq; adjZ++) {    //  
// ============================================== //

		ic0 = adjZ_eq[adjZ];
		ic2 = adj_number[ic0][1][6];

// #pragma omp parallel for private(j,k,k0,k2)
		for (i = n_buffer; i <= nx; i++) {
			for (j = n_buffer; j <= ny; j++) {
				for (k = 0; k < n_buffer; k++) {  

					k0 = NcubeZ+k;
					k2 = k;

					U1_[ic2][i][j][k2][0] = U1_[ic0][i][j][k0][0];
					U1_[ic2][i][j][k2][1] = U1_[ic0][i][j][k0][1];
					U1_[ic2][i][j][k2][2] = U1_[ic0][i][j][k0][2];
					U1_[ic2][i][j][k2][3] = U1_[ic0][i][j][k0][3];
					U1_[ic2][i][j][k2][4] = U1_[ic0][i][j][k0][4];
					
					k0 = n_buffer+NcubeZ+k;
					k2 = n_buffer+k;

					U1_[ic0][i][j][k0][0] = U1_[ic2][i][j][k2][0];
					U1_[ic0][i][j][k0][1] = U1_[ic2][i][j][k2][1];
					U1_[ic0][i][j][k0][2] = U1_[ic2][i][j][k2][2];
					U1_[ic0][i][j][k0][3] = U1_[ic2][i][j][k2][3];
					U1_[ic0][i][j][k0][4] = U1_[ic2][i][j][k2][4];
					
				}
			}
		}
// ============================================== //
	}                                             //
// ============================================== //





// ---- big cube to small cube in [Z+] direction for small cube ---- //


#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0) schedule(dynamic)

// ============================================== //
	for (adjZ = 1; adjZ <= nadjZ_bs_plus; adjZ++) {    //  
// ============================================== //
		
		ic0 = adjZ_bs_plus[adjZ];

		ic2 = adj_number[ic0][1][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  

						i0 = i/n_buffer+1;
						j0 = j/n_buffer+1;
						k0 = nz+k/n_buffer;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];

					}
				}
			}

		}


		ic2 = adj_number[ic0][2][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i/n_buffer+1+NcubeX/2;
						j0 = j/n_buffer+1;
						k0 = nz+k/n_buffer;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];



					}
				}
			}

		}


		ic2 = adj_number[ic0][3][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  

						i0 = i/n_buffer+1;
						j0 = j/n_buffer+1+NcubeY/2;
						k0 = nz+k/n_buffer;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];


					}
				}
			}

		}


		ic2 = adj_number[ic0][4][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  

						i0 = i/n_buffer+1+NcubeX/2;
						j0 = j/n_buffer+1+NcubeY/2;
						k0 = nz+k/n_buffer;

						U1_[ic2][i][j][k][0] = U1_[ic0][i0][j0][k0][0];
						U1_[ic2][i][j][k][1] = U1_[ic0][i0][j0][k0][1];
						U1_[ic2][i][j][k][2] = U1_[ic0][i0][j0][k0][2];
						U1_[ic2][i][j][k][3] = U1_[ic0][i0][j0][k0][3];
						U1_[ic2][i][j][k][4] = U1_[ic0][i0][j0][k0][4];

					}
				}
			}

		}

// ============================================== //
	}                                             //
// ============================================== //




// ---- small cube to big cube in [Z+] direction for small cube ---- //


#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,k2) schedule(dynamic)

// ============================================== //
	for (adjZ = 1; adjZ <= nadjZ_bs_minus; adjZ++) {    //  
// ============================================== //

		ic2 = adjZ_bs_minus[adjZ];

		ic0 = adj_number[ic2][1][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  

						i0 = i/n_buffer+1;
						j0 = j/n_buffer+1;

						k0 = nz+k+1;
						k2 = n_buffer+k/n_buffer;

						U1_[ic0][i][j][k0][0] = U1_[ic2][i0][j0][k2][0];
						U1_[ic0][i][j][k0][1] = U1_[ic2][i0][j0][k2][1];
						U1_[ic0][i][j][k0][2] = U1_[ic2][i0][j0][k2][2];
						U1_[ic0][i][j][k0][3] = U1_[ic2][i0][j0][k2][3];
						U1_[ic0][i][j][k0][4] = U1_[ic2][i0][j0][k2][4];


					}
				}
			}

		}


		ic0 = adj_number[ic2][2][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i/n_buffer+1+NcubeX/2;
						j0 = j/n_buffer+1;

						k0 = nz+k+1;
						k2 = n_buffer+k/n_buffer;

						U1_[ic0][i][j][k0][0] = U1_[ic2][i0][j0][k2][0];
						U1_[ic0][i][j][k0][1] = U1_[ic2][i0][j0][k2][1];
						U1_[ic0][i][j][k0][2] = U1_[ic2][i0][j0][k2][2];
						U1_[ic0][i][j][k0][3] = U1_[ic2][i0][j0][k2][3];
						U1_[ic0][i][j][k0][4] = U1_[ic2][i0][j0][k2][4];


					}
				}
			}

		}


		ic0 = adj_number[ic2][3][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i/n_buffer+1;
						j0 = j/n_buffer+1+NcubeY/2;

						k0 = nz+k+1;
						k2 = n_buffer+k/n_buffer;

						U1_[ic0][i][j][k0][0] = U1_[ic2][i0][j0][k2][0];
						U1_[ic0][i][j][k0][1] = U1_[ic2][i0][j0][k2][1];
						U1_[ic0][i][j][k0][2] = U1_[ic2][i0][j0][k2][2];
						U1_[ic0][i][j][k0][3] = U1_[ic2][i0][j0][k2][3];
						U1_[ic0][i][j][k0][4] = U1_[ic2][i0][j0][k2][4];


					}
				}
			}

		}


		ic0 = adj_number[ic2][4][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i/n_buffer+1+NcubeX/2;
						j0 = j/n_buffer+1+NcubeY/2;

						k0 = nz+k+1;
						k2 = n_buffer+k/n_buffer;

						U1_[ic0][i][j][k0][0] = U1_[ic2][i0][j0][k2][0];
						U1_[ic0][i][j][k0][1] = U1_[ic2][i0][j0][k2][1];
						U1_[ic0][i][j][k0][2] = U1_[ic2][i0][j0][k2][2];
						U1_[ic0][i][j][k0][3] = U1_[ic2][i0][j0][k2][3];
						U1_[ic0][i][j][k0][4] = U1_[ic2][i0][j0][k2][4];


					}
				}
			}

		}

// ============================================== //
	}                                             //
// ============================================== //





// ---- big cube to small cube in [Z-] direction ---- //

#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,i2,j2,k2) schedule(dynamic)

// ============================================== //
	for (adjZ = 1; adjZ <= nadjZ_bs_plus; adjZ++) {    //  
// ============================================== //

		ic0 = adjZ_bs_plus[adjZ];

		ic2 = adj_number[ic0][1][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  

						i0 = i;
						j0 = j;
						k0 = n_buffer+NcubeZ+k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = n_buffer*(k+1);

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic2 = adj_number[ic0][2][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i+NcubeX/2;
						j0 = j;
						k0 = n_buffer+NcubeZ+k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = n_buffer*(k+1);

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}


		ic2 = adj_number[ic0][3][6];

		if (ic2 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i;
						j0 = j+NcubeY/2;
						k0 = n_buffer+NcubeZ+k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = n_buffer*(k+1);

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}

		
		ic2 = adj_number[ic0][4][6];

		if (ic2 > 0) {

			//printf("%d\t%d\t%d\n",myid,ic2,ic0);

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i+NcubeX/2;
						j0 = j+NcubeY/2;
						k0 = n_buffer+NcubeZ+k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = n_buffer*(k+1);

						U1_[ic0][i0][j0][k0][0] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2  ][0]+
													  U1_[ic2][i2  ][j2+1][k2  ][0]+
													  U1_[ic2][i2  ][j2  ][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2  ][0]+
													  U1_[ic2][i2+1][j2  ][k2+1][0]+
													  U1_[ic2][i2  ][j2+1][k2+1][0]+
													  U1_[ic2][i2+1][j2+1][k2+1][0]);

						U1_[ic0][i0][j0][k0][1] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2  ][1]+
													  U1_[ic2][i2  ][j2+1][k2  ][1]+
													  U1_[ic2][i2  ][j2  ][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2  ][1]+
													  U1_[ic2][i2+1][j2  ][k2+1][1]+
													  U1_[ic2][i2  ][j2+1][k2+1][1]+
													  U1_[ic2][i2+1][j2+1][k2+1][1]);

						U1_[ic0][i0][j0][k0][2] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2  ][2]+
													  U1_[ic2][i2  ][j2+1][k2  ][2]+
													  U1_[ic2][i2  ][j2  ][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2  ][2]+
													  U1_[ic2][i2+1][j2  ][k2+1][2]+
													  U1_[ic2][i2  ][j2+1][k2+1][2]+
													  U1_[ic2][i2+1][j2+1][k2+1][2]);

						U1_[ic0][i0][j0][k0][3] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2  ][3]+
													  U1_[ic2][i2  ][j2+1][k2  ][3]+
													  U1_[ic2][i2  ][j2  ][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2  ][3]+
													  U1_[ic2][i2+1][j2  ][k2+1][3]+
													  U1_[ic2][i2  ][j2+1][k2+1][3]+
													  U1_[ic2][i2+1][j2+1][k2+1][3]);

						U1_[ic0][i0][j0][k0][4] = 0.125*(U1_[ic2][i2  ][j2  ][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2  ][4]+
													  U1_[ic2][i2  ][j2+1][k2  ][4]+
													  U1_[ic2][i2  ][j2  ][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2  ][4]+
													  U1_[ic2][i2+1][j2  ][k2+1][4]+
													  U1_[ic2][i2  ][j2+1][k2+1][4]+
													  U1_[ic2][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}
		

// ============================================== //
	}                                             //
// ============================================== //




	
// ---- small cube to big cube in [Z+] direction for big cube ---- //

#pragma omp parallel for private(ic0, ic2,i,j,k,i0,j0,k0,i2,j2,k2) schedule(dynamic)

// ============================================== //
	for (adjZ = 1; adjZ <= nadjZ_bs_minus; adjZ++) {    //  
// ============================================== //

		ic2 = adjZ_bs_minus[adjZ];

		ic0 = adj_number[ic2][1][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  

						i0 = i;
						j0 = j;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = nz-2*n_buffer+n_buffer*(k+1)-1;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);


					}
				}
			}

		}


		ic0 = adj_number[ic2][2][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i+NcubeX/2;
						j0 = j;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = nz-2*n_buffer+n_buffer*(k+1)-1;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);
					}
				}
			}

		}


		ic0 = adj_number[ic2][3][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i;
						j0 = j+NcubeY/2;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = nz-2*n_buffer+n_buffer*(k+1)-1;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);
					}
				}
			}

		}


		ic0 = adj_number[ic2][4][5];

		if (ic0 > 0) {

			for (i = n_buffer; i <= (nx-n_buffer+1)/2+1; i++) {
				for (j = n_buffer; j <= (ny-n_buffer+1)/2+1; j++) {
					for (k = 0; k < n_buffer; k++) {  


						i0 = i+NcubeX/2;
						j0 = j+NcubeY/2;
						k0 = k;

						i2 = (i-n_buffer+1)*2;
						j2 = (j-n_buffer+1)*2;
						k2 = nz-2*n_buffer+n_buffer*(k+1)-1;

						U1_[ic2][i0][j0][k0][0] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2  ][0]+
													  U1_[ic0][i2  ][j2+1][k2  ][0]+
													  U1_[ic0][i2  ][j2  ][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2  ][0]+
													  U1_[ic0][i2+1][j2  ][k2+1][0]+
													  U1_[ic0][i2  ][j2+1][k2+1][0]+
													  U1_[ic0][i2+1][j2+1][k2+1][0]);

						U1_[ic2][i0][j0][k0][1] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2  ][1]+
													  U1_[ic0][i2  ][j2+1][k2  ][1]+
													  U1_[ic0][i2  ][j2  ][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2  ][1]+
													  U1_[ic0][i2+1][j2  ][k2+1][1]+
													  U1_[ic0][i2  ][j2+1][k2+1][1]+
													  U1_[ic0][i2+1][j2+1][k2+1][1]);

						U1_[ic2][i0][j0][k0][2] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2  ][2]+
													  U1_[ic0][i2  ][j2+1][k2  ][2]+
													  U1_[ic0][i2  ][j2  ][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2  ][2]+
													  U1_[ic0][i2+1][j2  ][k2+1][2]+
													  U1_[ic0][i2  ][j2+1][k2+1][2]+
													  U1_[ic0][i2+1][j2+1][k2+1][2]);

						U1_[ic2][i0][j0][k0][3] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2  ][3]+
													  U1_[ic0][i2  ][j2+1][k2  ][3]+
													  U1_[ic0][i2  ][j2  ][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2  ][3]+
													  U1_[ic0][i2+1][j2  ][k2+1][3]+
													  U1_[ic0][i2  ][j2+1][k2+1][3]+
													  U1_[ic0][i2+1][j2+1][k2+1][3]);

						U1_[ic2][i0][j0][k0][4] = 0.125*(U1_[ic0][i2  ][j2  ][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2  ][4]+
													  U1_[ic0][i2  ][j2+1][k2  ][4]+
													  U1_[ic0][i2  ][j2  ][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2  ][4]+
													  U1_[ic0][i2+1][j2  ][k2+1][4]+
													  U1_[ic0][i2  ][j2+1][k2+1][4]+
													  U1_[ic0][i2+1][j2+1][k2+1][4]);

					}
				}
			}

		}

		
		
		

				
					
// ============================================== //
	}                                             //
// ============================================== //
		
	
	
	
	
	
	//stop_collection("region4");
	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		MPI_Waitany(ncpu_eq, Rreq_eq, &index, istat);

	}
	
	//start_collection("region5");

	
// ==================================================================================================== //
// ============================== # Unpackage for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	

	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		
		ii = ii+Ncube_Ncpu_eq[icpu_neig_eq];
		count_index = (ii-Ncube_Ncpu_eq[icpu_neig_eq])*zone;
		
	
		L1 = 0*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L2 = 1*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L3 = 2*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L4 = 3*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;
		L5 = 4*Ncube_Ncpu_eq[icpu_neig_eq]*zone+5*count_index;


#pragma omp parallel for private(iicube,icube,i,j,k,iL,k0,j0,i0) 	
		for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_eq[icpu_neig_eq]+icube_send;
			icube = Rcube_Ncpu_eq[iicube];

// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Rdir_eq[iicube] == 6) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k0 = k;
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}

			}    // ---- if (Rdir_eq[iicube] == 6)  ---- //

// ------------------------------------- [Z+] direction ------------------------------------- //


// ------------------------------------- [Z-] direction ------------------------------------- //

			else if (Rdir_eq[iicube] == 5) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k0 = n_buffer+NcubeZ+k;
							
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_eq[iicube] == 5)  ---- //

// ------------------------------------- [Z-] direction ------------------------------------- //


// ------------------------------------- [Y+] direction ------------------------------------- //

			else if (Rdir_eq[iicube] == 4) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							j0 = j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_eq[iicube] == 4) ---- //


// ------------------------------------- [Y+] direction ------------------------------------- //


// ------------------------------------- [Y-] direction ------------------------------------- //

			else if (Rdir_eq[iicube] == 3) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							j0 = n_buffer+NcubeY+j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_eq[iicube] == 3) ---- //


// ------------------------------------- [Y-] direction ------------------------------------- //


// ------------------------------------- [X+] direction ------------------------------------- //

			else if (Rdir_eq[iicube] == 2) {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
							i0 = i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_eq[iicube] == 2) ---- //


// ------------------------------------- [X+] direction ------------------------------------- //


// ------------------------------------- [X-] direction ------------------------------------- //

			else {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
							i0 = n_buffer+NcubeX+i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) ---- //

#pragma omp barrier
		
	}    // ---- for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) ---- //

// ============================== # Unpackage for sending in Current CPU ============================== //
// ==================================================================================================== //


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//stop_collection("region5");

	for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) {

		MPI_Waitany(ncpu_bs, Rreq_bs, &index, istat);

	}
	



	


// ==================================================================================================== //
// ============================== # Unpackage for sending in Current CPU ============================== //

//start_collection("region6");

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	

	for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) {

	
		ii = ii+Ncube_Ncpu_bs[icpu_neig_bs];
		count_index = (ii-Ncube_Ncpu_bs[icpu_neig_bs])*zone/4;
		
		
		L1 = 0*Ncube_Ncpu_bs[icpu_neig_bs]*zone/4+5*count_index;
		L2 = 1*Ncube_Ncpu_bs[icpu_neig_bs]*zone/4+5*count_index;
		L3 = 2*Ncube_Ncpu_bs[icpu_neig_bs]*zone/4+5*count_index;
		L4 = 3*Ncube_Ncpu_bs[icpu_neig_bs]*zone/4+5*count_index;
		L5 = 4*Ncube_Ncpu_bs[icpu_neig_bs]*zone/4+5*count_index;


#pragma omp parallel for private(iicube,icube,i,j,k,iL,k0,j0,i0) 	

		for (icube_send = 1; icube_send <= Ncube_Ncpu_bs[icpu_neig_bs]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_bs[icpu_neig_bs]+icube_send;
			icube = Rcube_Ncpu_bs[iicube];

// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Rdir_bs[iicube] == 6) {

// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 1) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i;
								j0 = j;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//




// -----------------------------------------------------------------------------//				
				if (RadjN_bs[iicube] == 2) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i+NcubeX/2;
								j0 = j;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//




// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 3) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i;
								j0 = j+NcubeY/2;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 4) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i+NcubeX/2;
								j0 = j+NcubeY/2;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								//printf("%f\t",U1_[icube][i0][j0][k0]);
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- if (Rdir_bs[iicube] == 6)  ---- //

// ------------------------------------- [Z+] direction ------------------------------------- //

// ------------------------------------- [Z-] direction ------------------------------------- //

			else if (Rdir_bs[iicube] == 5) {


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 1) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i;
								j0 = j;
								k0 = n_buffer+NcubeZ+k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 2) {


					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i+NcubeX/2;
								j0 = j;
								k0 = n_buffer+NcubeZ+k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 3) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i;
								j0 = j+NcubeY/2;
								k0 = n_buffer+NcubeZ+k;


								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 4) {

					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = n_buffer; j < NcubeY/2+n_buffer; j++) {
							for (k = 0; k < n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*NcubeY/2*n_buffer+(j-n_buffer)*n_buffer+k;

								i0 = i+NcubeX/2;
								j0 = j+NcubeY/2;
								k0 = n_buffer+NcubeZ+k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


			}    // ---- else if (Rdir_bs[iicube] == 5)  ---- //

// ------------------------------------- [Z-] direction ------------------------------------- //


// ------------------------------------- [Y+] direction ------------------------------------- //

			else if (Rdir_bs[iicube] == 4) {


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 1) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = j;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 2) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i+NcubeX/2;
								j0 = j;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 3) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = j;
								k0 = k+NcubeZ/2;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 4) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i+NcubeX/2;
								j0 = j;
								k0 = k+NcubeZ/2;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- else if (Rdir_bs[iicube] == 4) ---- //


// ------------------------------------- [Y+] direction ------------------------------------- //


// ------------------------------------- [Y-] direction ------------------------------------- //

			else if (Rdir_bs[iicube] == 3) {

				
// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 1) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = n_buffer+NcubeY+j;
								k0 = k;

								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 2) {

					

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i+NcubeX/2;
								j0 = n_buffer+NcubeY+j;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 3) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = n_buffer+NcubeY+j;
								k0 = k+NcubeZ/2;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 4) {

					
					for (i = n_buffer; i < NcubeX/2+n_buffer; i++) {
						for (j = 0; j < n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+(i-n_buffer)*n_buffer*NcubeZ/2+j*NcubeZ/2+(k-n_buffer);

								i0 = i+NcubeX/2;
								j0 = n_buffer+NcubeY+j;
								k0 = k+NcubeZ/2;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//

			}    // ---- else if (Rdir_bs[iicube] == 3) ---- //


// ------------------------------------- [Y-] direction ------------------------------------- //


// ------------------------------------- [X+] direction ------------------------------------- //

			else if (Rdir_bs[iicube] == 2) {

				

// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 1) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = j;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 2) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);
								
								i0 = i;
								j0 = j+NcubeY/2;
								k0 = k;

								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 3) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = j;
								k0 = k+NcubeZ/2;

								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 4) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = i;
								j0 = j+NcubeY/2;
								k0 = k+NcubeZ/2;

								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//



			}    // ---- else if (Rdir_bs[iicube] == 2) ---- //


// ------------------------------------- [X+] direction ------------------------------------- //


// ------------------------------------- [X-] direction ------------------------------------- //

			else {

	

// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 1) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = n_buffer+NcubeX+i;
								j0 = j;
								k0 = k;
								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 2) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = n_buffer+NcubeX+i;
								j0 = j+NcubeY/2;
								k0 = k;

								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 3) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = n_buffer+NcubeX+i;
								j0 = j;
								k0 = k+NcubeZ/2;

								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];
								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


// -----------------------------------------------------------------------------//
				if (RadjN_bs[iicube] == 4) {

					
					for (i = 0 ; i < n_buffer; i++) {
						for (j = n_buffer; j < NcubeX/2+n_buffer; j++) {
							for (k = n_buffer; k < NcubeZ/2+n_buffer; k++) {  

								iL = (icube_send-1)*zone/4+i*NcubeY/2*NcubeZ/2+(j-n_buffer)*NcubeZ/2+(k-n_buffer);

								i0 = n_buffer+NcubeX+i;
								j0 = j+NcubeY/2;
								k0 = k+NcubeZ/2;

								
								U1_[icube][i0][j0][k0][0] = recv_data_curr_bs[L1+iL];
								U1_[icube][i0][j0][k0][1] = recv_data_curr_bs[L2+iL];
								U1_[icube][i0][j0][k0][2] = recv_data_curr_bs[L3+iL];
								U1_[icube][i0][j0][k0][3] = recv_data_curr_bs[L4+iL];
								U1_[icube][i0][j0][k0][4] = recv_data_curr_bs[L5+iL];

								
							}
						}
					}

				}
// -----------------------------------------------------------------------------//


			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_bs[icpu_neig_bs]; icube_send++) ---- //

#pragma omp barrier
		
	}    // ---- for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) ---- //


// ============================== # Unpackage for sending in Current CPU ============================== //
// ==================================================================================================== //

	
	

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//stop_collection("region6");





	for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) {

		MPI_Waitany(ncpu_sb, Rreq_sb, &index, istat);

	}
	
	
//start_collection("region7");

	 
	
// ==================================================================================================== //
// ============================== # Unpackage for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	

	for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) {

		
		ii = ii+Ncube_Ncpu_sb[icpu_neig_sb];
		count_index = (ii-Ncube_Ncpu_sb[icpu_neig_sb])*zone;
		
	
		L1 = 0*Ncube_Ncpu_sb[icpu_neig_sb]*zone+5*count_index;
		L2 = 1*Ncube_Ncpu_sb[icpu_neig_sb]*zone+5*count_index;
		L3 = 2*Ncube_Ncpu_sb[icpu_neig_sb]*zone+5*count_index;
		L4 = 3*Ncube_Ncpu_sb[icpu_neig_sb]*zone+5*count_index;
		L5 = 4*Ncube_Ncpu_sb[icpu_neig_sb]*zone+5*count_index;

#pragma omp parallel for private(iicube,icube,i,j,k,iL,k0,j0,i0) 	

		for (icube_send = 1; icube_send <= Ncube_Ncpu_sb[icpu_neig_sb]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_sb[icpu_neig_sb]+icube_send;
			icube = Rcube_Ncpu_sb[iicube];


// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Rdir_sb[iicube] == 6) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							
							k0 = k;
							

							U1_[icube][i][j][k0][0] = recv_data_curr_sb[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_sb[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_sb[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_sb[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_sb[L5+iL];

						}
					}
				}

			}    // ---- if (Rdir_sb[iicube] == 6)  ---- //

// ------------------------------------- [Z+] direction ------------------------------------- //

// ------------------------------------- [Z-] direction ------------------------------------- //

			else if (Rdir_sb[iicube] == 5) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							
							k0 = n_buffer+NcubeZ+k;
							
							U1_[icube][i][j][k0][0] = recv_data_curr_sb[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_sb[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_sb[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_sb[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_sb[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_sb[iicube] == 5)  ---- //

// ------------------------------------- [Z-] direction ------------------------------------- //


// ------------------------------------- [Y+] direction ------------------------------------- //

			else if (Rdir_sb[iicube] == 4) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							
							j0 = j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_sb[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_sb[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_sb[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_sb[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_sb[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_sb[iicube] == 4) ---- //


// ------------------------------------- [Y+] direction ------------------------------------- //


// ------------------------------------- [Y-] direction ------------------------------------- //

			else if (Rdir_sb[iicube] == 3) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+(i-n_buffer)*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							
							j0 = n_buffer+NcubeY+j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_sb[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_sb[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_sb[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_sb[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_sb[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_sb[iicube] == 3) ---- //


// ------------------------------------- [Y-] direction ------------------------------------- //


// ------------------------------------- [X+] direction ------------------------------------- //

			else if (Rdir_sb[iicube] == 2) {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
							
							i0 = i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_sb[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_sb[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_sb[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_sb[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_sb[L5+iL];
							
						}
					}
				}

			}    // ---- else if (Rdir_sb[iicube] == 2) ---- //


// ------------------------------------- [X+] direction ------------------------------------- //


// ------------------------------------- [X-] direction ------------------------------------- //

			else {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*NcubeZ+(j-n_buffer)*NcubeZ+(k-n_buffer);
							
							i0 = n_buffer+NcubeX+i;

							U1_[icube][i0][j][k][0] = recv_data_curr_sb[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_sb[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_sb[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_sb[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_sb[L5+iL];

						}
					}
				}

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_sb[icpu_neig_sb]; icube_send++) ---- //
		
#pragma omp barrier

	}    // ---- for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) ---- //

// ============================== # Unpackage for sending in Current CPU ============================== //
// ==================================================================================================== //


//stop_collection("region7");

	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		MPI_Waitany(ncpu_eq, Sreq_eq, &index, istat);

	}

	
	for (icpu_neig_sb = 0;  icpu_neig_sb < ncpu_sb; icpu_neig_sb++) {     
				
		MPI_Waitany(ncpu_sb, Sreq_sb, &index, istat);

	}
	
	
	for (icpu_neig_bs = 0;  icpu_neig_bs < ncpu_bs; icpu_neig_bs++) {

		MPI_Waitany(ncpu_bs, Sreq_bs, &index, istat);

	}
	
	
	

	
#pragma omp parallel for private(i,j,k)
	for (icube = 1; icube < ncube; icube++) { 

		for (i = n_buffer-1; i <= nxx; i++) {

			U1_[icube][i][n_buffer-1][n_buffer-1][0] =  0.5*(U1_[icube][i][n_buffer-1][n_buffer][0]+U1_[icube][i][n_buffer][n_buffer-1][0]);
			U1_[icube][i][n_buffer-1][n_buffer-1][1] =  0.5*(U1_[icube][i][n_buffer-1][n_buffer][1]+U1_[icube][i][n_buffer][n_buffer-1][1]);
			U1_[icube][i][n_buffer-1][n_buffer-1][2] =  0.5*(U1_[icube][i][n_buffer-1][n_buffer][2]+U1_[icube][i][n_buffer][n_buffer-1][2]);
			U1_[icube][i][n_buffer-1][n_buffer-1][3] =  0.5*(U1_[icube][i][n_buffer-1][n_buffer][3]+U1_[icube][i][n_buffer][n_buffer-1][3]);
			U1_[icube][i][n_buffer-1][n_buffer-1][4] =  0.5*(U1_[icube][i][n_buffer-1][n_buffer][4]+U1_[icube][i][n_buffer][n_buffer-1][4]);

			U1_[icube][i][nyy][n_buffer-1][0] =  0.5*(U1_[icube][i][nyy][n_buffer][0]+U1_[icube][i][ny][n_buffer-1][0]);
			U1_[icube][i][nyy][n_buffer-1][1] =  0.5*(U1_[icube][i][nyy][n_buffer][1]+U1_[icube][i][ny][n_buffer-1][1]);
			U1_[icube][i][nyy][n_buffer-1][2] =  0.5*(U1_[icube][i][nyy][n_buffer][2]+U1_[icube][i][ny][n_buffer-1][2]);
			U1_[icube][i][nyy][n_buffer-1][3] =  0.5*(U1_[icube][i][nyy][n_buffer][3]+U1_[icube][i][ny][n_buffer-1][3]);
			U1_[icube][i][nyy][n_buffer-1][4] =  0.5*(U1_[icube][i][nyy][n_buffer][4]+U1_[icube][i][ny][n_buffer-1][4]);

			U1_[icube][i][n_buffer-1][nzz][0] =  0.5*(U1_[icube][i][n_buffer-1][nz][0]+U1_[icube][i][n_buffer][nzz][0]);
			U1_[icube][i][n_buffer-1][nzz][1] =  0.5*(U1_[icube][i][n_buffer-1][nz][1]+U1_[icube][i][n_buffer][nzz][1]);
			U1_[icube][i][n_buffer-1][nzz][2] =  0.5*(U1_[icube][i][n_buffer-1][nz][2]+U1_[icube][i][n_buffer][nzz][2]);
			U1_[icube][i][n_buffer-1][nzz][3] =  0.5*(U1_[icube][i][n_buffer-1][nz][3]+U1_[icube][i][n_buffer][nzz][3]);
			U1_[icube][i][n_buffer-1][nzz][4] =  0.5*(U1_[icube][i][n_buffer-1][nz][4]+U1_[icube][i][n_buffer][nzz][4]);

			U1_[icube][i][nyy][nzz][0] =  0.5*(U1_[icube][i][nyy][nz][0]+U1_[icube][i][ny][nzz][0]);
			U1_[icube][i][nyy][nzz][1] =  0.5*(U1_[icube][i][nyy][nz][1]+U1_[icube][i][ny][nzz][1]);
			U1_[icube][i][nyy][nzz][2] =  0.5*(U1_[icube][i][nyy][nz][2]+U1_[icube][i][ny][nzz][2]);
			U1_[icube][i][nyy][nzz][3] =  0.5*(U1_[icube][i][nyy][nz][3]+U1_[icube][i][ny][nzz][3]);
			U1_[icube][i][nyy][nzz][4] =  0.5*(U1_[icube][i][nyy][nz][4]+U1_[icube][i][ny][nzz][4]);

		}




		for (j = n_buffer-1; j <= nyy; j++) {

			U1_[icube][n_buffer-1][j][n_buffer-1][0] =  0.5*(U1_[icube][n_buffer-1][j][n_buffer][0]+U1_[icube][n_buffer][j][n_buffer-1][0]);
			U1_[icube][n_buffer-1][j][n_buffer-1][1] =  0.5*(U1_[icube][n_buffer-1][j][n_buffer][1]+U1_[icube][n_buffer][j][n_buffer-1][1]);
			U1_[icube][n_buffer-1][j][n_buffer-1][2] =  0.5*(U1_[icube][n_buffer-1][j][n_buffer][2]+U1_[icube][n_buffer][j][n_buffer-1][2]);
			U1_[icube][n_buffer-1][j][n_buffer-1][3] =  0.5*(U1_[icube][n_buffer-1][j][n_buffer][3]+U1_[icube][n_buffer][j][n_buffer-1][3]);
			U1_[icube][n_buffer-1][j][n_buffer-1][4] =  0.5*(U1_[icube][n_buffer-1][j][n_buffer][4]+U1_[icube][n_buffer][j][n_buffer-1][4]);

			U1_[icube][n_buffer-1][j][nzz][0] =  0.5*(U1_[icube][n_buffer][j][nzz][0]+U1_[icube][n_buffer-1][j][nz][0]);
			U1_[icube][n_buffer-1][j][nzz][1] =  0.5*(U1_[icube][n_buffer][j][nzz][1]+U1_[icube][n_buffer-1][j][nz][1]);
			U1_[icube][n_buffer-1][j][nzz][2] =  0.5*(U1_[icube][n_buffer][j][nzz][2]+U1_[icube][n_buffer-1][j][nz][2]);
			U1_[icube][n_buffer-1][j][nzz][3] =  0.5*(U1_[icube][n_buffer][j][nzz][3]+U1_[icube][n_buffer-1][j][nz][3]);
			U1_[icube][n_buffer-1][j][nzz][4] =  0.5*(U1_[icube][n_buffer][j][nzz][4]+U1_[icube][n_buffer-1][j][nz][4]);

			U1_[icube][nxx][j][n_buffer-1][0] =  0.5*(U1_[icube][nx][j][n_buffer-1][0]+U1_[icube][nxx][j][n_buffer][0]);
			U1_[icube][nxx][j][n_buffer-1][1] =  0.5*(U1_[icube][nx][j][n_buffer-1][1]+U1_[icube][nxx][j][n_buffer][1]);
			U1_[icube][nxx][j][n_buffer-1][2] =  0.5*(U1_[icube][nx][j][n_buffer-1][2]+U1_[icube][nxx][j][n_buffer][2]);
			U1_[icube][nxx][j][n_buffer-1][3] =  0.5*(U1_[icube][nx][j][n_buffer-1][3]+U1_[icube][nxx][j][n_buffer][3]);
			U1_[icube][nxx][j][n_buffer-1][4] =  0.5*(U1_[icube][nx][j][n_buffer-1][4]+U1_[icube][nxx][j][n_buffer][4]);

			U1_[icube][nxx][j][nzz][0] =  0.5*(U1_[icube][nx][j][nzz][0]+U1_[icube][nxx][j][nz][0]);
			U1_[icube][nxx][j][nzz][1] =  0.5*(U1_[icube][nx][j][nzz][1]+U1_[icube][nxx][j][nz][1]);
			U1_[icube][nxx][j][nzz][2] =  0.5*(U1_[icube][nx][j][nzz][2]+U1_[icube][nxx][j][nz][2]);
			U1_[icube][nxx][j][nzz][3] =  0.5*(U1_[icube][nx][j][nzz][3]+U1_[icube][nxx][j][nz][3]);
			U1_[icube][nxx][j][nzz][4] =  0.5*(U1_[icube][nx][j][nzz][4]+U1_[icube][nxx][j][nz][4]);

		}



		for (k = n_buffer-1; k <= nzz; k++) {

			U1_[icube][n_buffer-1][n_buffer-1][k][0] =  0.5*(U1_[icube][n_buffer-1][n_buffer][k][0]+U1_[icube][n_buffer][n_buffer-1][k][0]);
			U1_[icube][n_buffer-1][n_buffer-1][k][1] =  0.5*(U1_[icube][n_buffer-1][n_buffer][k][1]+U1_[icube][n_buffer][n_buffer-1][k][1]);
			U1_[icube][n_buffer-1][n_buffer-1][k][2] =  0.5*(U1_[icube][n_buffer-1][n_buffer][k][2]+U1_[icube][n_buffer][n_buffer-1][k][2]);
			U1_[icube][n_buffer-1][n_buffer-1][k][3] =  0.5*(U1_[icube][n_buffer-1][n_buffer][k][3]+U1_[icube][n_buffer][n_buffer-1][k][3]);
			U1_[icube][n_buffer-1][n_buffer-1][k][4] =  0.5*(U1_[icube][n_buffer-1][n_buffer][k][4]+U1_[icube][n_buffer][n_buffer-1][k][4]);

			U1_[icube][n_buffer-1][nyy][k][0] =  0.5*(U1_[icube][n_buffer][nyy][k][0]+U1_[icube][n_buffer-1][ny][k][0]);
			U1_[icube][n_buffer-1][nyy][k][1] =  0.5*(U1_[icube][n_buffer][nyy][k][1]+U1_[icube][n_buffer-1][ny][k][1]);
			U1_[icube][n_buffer-1][nyy][k][2] =  0.5*(U1_[icube][n_buffer][nyy][k][2]+U1_[icube][n_buffer-1][ny][k][2]);
			U1_[icube][n_buffer-1][nyy][k][3] =  0.5*(U1_[icube][n_buffer][nyy][k][3]+U1_[icube][n_buffer-1][ny][k][3]);
			U1_[icube][n_buffer-1][nyy][k][4] =  0.5*(U1_[icube][n_buffer][nyy][k][4]+U1_[icube][n_buffer-1][ny][k][4]);

			U1_[icube][nxx][n_buffer-1][k][0] =  0.5*(U1_[icube][nx][n_buffer-1][k][0]+U1_[icube][nxx][n_buffer][k][0]);
			U1_[icube][nxx][n_buffer-1][k][1] =  0.5*(U1_[icube][nx][n_buffer-1][k][1]+U1_[icube][nxx][n_buffer][k][1]);
			U1_[icube][nxx][n_buffer-1][k][2] =  0.5*(U1_[icube][nx][n_buffer-1][k][2]+U1_[icube][nxx][n_buffer][k][2]);
			U1_[icube][nxx][n_buffer-1][k][3] =  0.5*(U1_[icube][nx][n_buffer-1][k][3]+U1_[icube][nxx][n_buffer][k][3]);
			U1_[icube][nxx][n_buffer-1][k][4] =  0.5*(U1_[icube][nx][n_buffer-1][k][4]+U1_[icube][nxx][n_buffer][k][4]);

			U1_[icube][nxx][nyy][k][0] =  0.5*(U1_[icube][nx][nyy][k][0]+U1_[icube][nxx][ny][k][0]);
			U1_[icube][nxx][nyy][k][1] =  0.5*(U1_[icube][nx][nyy][k][1]+U1_[icube][nxx][ny][k][1]);
			U1_[icube][nxx][nyy][k][2] =  0.5*(U1_[icube][nx][nyy][k][2]+U1_[icube][nxx][ny][k][2]);
			U1_[icube][nxx][nyy][k][3] =  0.5*(U1_[icube][nx][nyy][k][3]+U1_[icube][nxx][ny][k][3]);
			U1_[icube][nxx][nyy][k][4] =  0.5*(U1_[icube][nx][nyy][k][4]+U1_[icube][nxx][ny][k][4]);
			
		}
		
	}
	

}                                          








