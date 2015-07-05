



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

extern int Ncpu_eq;
	
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
	
	
	MPI_Request (*Sreq_eq) = new MPI_Request[ncpu_eq];
	MPI_Request (*Rreq_eq) = new MPI_Request[ncpu_eq];

	
	double (*tempFWS)[X_size][Y_size][Z_size] = new double[ncube][X_size][Y_size][Z_size]; 

	
	for (icube = 1; icube < ncube; icube++) {  

		if (csl[icube] == 0)
#pragma omp parallel for private(j,k)
			
			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
					for (k = 0; k <= nzzz; k++) { 
					
							tempFWS[icube][i][j][k] = FWS[icube][i][j][k];  

					}
				}
			}

	}


	



// ================================================================================================== //
// ============================== # Package for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		L1 = count_index;
		
		iL = -1;
		
		for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) {  

			iicube = iicube + 1;
			icube = Scube_Ncpu_eq[iicube];

// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Sdir_eq[iicube] == 6) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							count_index = count_index+1;
							iL = iL + 1;
							k2 = NcubeZ+k;
							
							send_data_curr_eq[L1+iL] = tempFWS[icube][i][j][k2];
					
							
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

							count_index = count_index+1;
							iL = iL + 1;
							k2 = n_buffer+k;
							
							send_data_curr_eq[L1+iL] = tempFWS[icube][i][j][k2];
							
							
							//if (myid == 1 & icube == 1) printf("%f\t",U1_[icube][i][j][k2]);


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

							count_index = count_index+1;
							iL = iL + 1;
							j2 = NcubeY+j;
							
							send_data_curr_eq[L1+iL] = tempFWS[icube][i][j2][k];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							j2 = n_buffer+j;
							
							send_data_curr_eq[L1+iL] = tempFWS[icube][i][j2][k];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							i2 = NcubeX+i;
							
							send_data_curr_eq[L1+iL] = tempFWS[icube][i2][j][k];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							i2 = n_buffer+i;
							
							send_data_curr_eq[L1+iL] = tempFWS[icube][i2][j][k];
							
							
						}
					}
				}

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) ---- //

	}    // ---- for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) ---- //

	

// ============================== # Package for sending in Current CPU ============================== //
// ================================================================================================== //



	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Status istat[8];
	MPI_Request requ1, reqps1;

	
	

	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		istart = ist_eq[icpu_neig_eq]*NcubeX*NcubeY*n_buffer;

		idest = neighbor_cpu_eq[icpu_neig_eq];

		icount = Ncube_Ncpu_eq[icpu_neig_eq]*NcubeX*NcubeY*n_buffer;

		itag = 0;

		MPI_Isend((void *)&send_data_curr_eq[istart], icount, MPI_DOUBLE, idest, itag, comm, &Sreq_eq[icpu_neig_eq]);

	}


	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		istart = ist_eq[icpu_neig_eq]*NcubeX*NcubeY*n_buffer;

		isrc = neighbor_cpu_eq[icpu_neig_eq];

		icount = Ncube_Ncpu_eq[icpu_neig_eq]*NcubeX*NcubeY*n_buffer;

		itag = 0;

		MPI_Irecv((void *)&recv_data_curr_eq[istart], icount, MPI_DOUBLE, isrc, itag, comm, &Rreq_eq[icpu_neig_eq]);


	}

	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		int index;

		MPI_Waitany(ncpu_eq, Rreq_eq, &index, istat);

	}
	
	

// ==================================================================================================== //
// ============================== # Unpackage for sending in Current CPU ============================== //

	count_index = 0;
	ii = jj = kk = iicube = iL =0;
	L1 = L2 = L3 = L4 = L5 = 0;
	

	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		L1 = count_index;
		
		iL = -1;

		for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) {  

			iicube = iicube + 1;
			icube = Rcube_Ncpu_eq[iicube];


// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Rdir_eq[iicube] == 6) {

				for (i = n_buffer; i <= nx; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							count_index = count_index+1;
							iL = iL + 1;
							k0 = k;
							
							tempFWS[icube][i][j][k0] = recv_data_curr_eq[L1+iL];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							k0 = n_buffer+NcubeZ+k;
							
							
							tempFWS[icube][i][j][k0] = recv_data_curr_eq[L1+iL];
							
							
							
							//if (myid == 0 & icube == 2) printf("%f\t",U1_[icube][i][j][k0]);

							
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

							count_index = count_index+1;
							iL = iL + 1;
							j0 = j;
							
							tempFWS[icube][i][j0][k] = recv_data_curr_eq[L1+iL];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							j0 = n_buffer+NcubeY+j;
							
							tempFWS[icube][i][j0][k] = recv_data_curr_eq[L1+iL];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							i0 = i;
							
							tempFWS[icube][i0][j][k] = recv_data_curr_eq[L1+iL];
							
							
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

							count_index = count_index+1;
							iL = iL + 1;
							i0 = n_buffer+NcubeX+i;
							
							tempFWS[icube][i0][j][k] = recv_data_curr_eq[L1+iL];
							
							
						}
					}
				}

			}    // ---- else ---- //

// ------------------------------------- [X-] direction ------------------------------------- //

		}    // ---- for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) ---- //

	}    // ---- for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) ---- //
	
// ============================== # Unpackage for sending in Current CPU ============================== //
// ==================================================================================================== //

	
	//MPI_Wait(Sreq_eq,istat);

	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		int index;

		MPI_Waitany(ncpu_eq, Sreq_eq, &index, istat);

	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// ---- the same size connection in X direction ---- //
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
					
					tempFWS[ic2][i2][j][k] = tempFWS[ic0][i0][j][k];
					
					
					i0 = n_buffer+NcubeX+i;
					i2 = n_buffer+i;

					tempFWS[ic0][i0][j][k] = tempFWS[ic2][i2][j][k];
					
				}
			}
		}
// ============================================== //
	}                                             //
// ============================================== //





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

					tempFWS[ic2][i][j2][k] = tempFWS[ic0][i][j0][k];
					
					
					j0 = n_buffer+NcubeY+j;
					j2 = n_buffer+j;

					tempFWS[ic0][i][j0][k] = tempFWS[ic2][i][j2][k];
					
					
				}
			}
		}
// ============================================== //
	}                                             //
// ============================================== //










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

					tempFWS[ic2][i][j][k2] = tempFWS[ic0][i][j][k0];
					
					
					k0 = n_buffer+NcubeZ+k;
					k2 = n_buffer+k;

					tempFWS[ic0][i][j][k0] = tempFWS[ic2][i][j][k2];
					
					
				}
			}
		}
// ============================================== //
	}                                             //
// ============================================== //


	for (icube = 1; icube < ncube; icube++) { 

		if (csl[icube] == 0) {

			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
					for (k = 0; k <= nzzz; k++) { 

						if (tempFWS[icube][i][j][k] > 0.5)

							FWS[icube][i][j][k] = 1;

						else if (tempFWS[icube][i][j][k] < -0.5)

							FWS[icube][i][j][k] = -1;

						else

							FWS[icube][i][j][k] = 0;

					}
				}
			}

		}  // ---- if (csl[icube] == 0) ---- //

	}  







// -------------------------------------------- Hole treatment -------------------------------------------- //
	
	int end = 0;
	count_index = 0;

	for (icube = 1; icube < ncube; icube++) { 

		if (csl[icube] == 0) {

			for (int iFWS = 1; iFWS < end; iFWS ++) {

				
				for (i = n_buffer-1; i <= nxx; i++) {
					for (j = n_buffer-1; j <= nyy; j++) {
						for (k = n_buffer-1; k <= nzz; k++) { 

							if ( FWS[icube][i][j][k] == IFLUID  && (
								(FWS[icube][i+1][j][k] + FWS[icube][i-1][j][k] <  IFLUID ) |
								(FWS[icube][i][j+1][k] + FWS[icube][i][j-1][k] <  IFLUID ) |
								(FWS[icube][i][j][k+1] + FWS[icube][i][j][k-1] <  IFLUID ) 
								)) {
								FWS[icube][i][j][k] = ISOLID;
								count_index = 1;
								
							}

						}
					}
				}

				if(count_index > 0) {
					
					end = end;
					count_index = 0;

				}

			}    // ---- for (int iFWS = 1; iFWS < 100000; iFWS ++) ---- //


	//		end = 2;
	//		for (int iFWS = 1; iFWS < end; iFWS ++) {

	//			
	//			for (i = n_buffer-1; i <= nx; i++) {
	//				for (j = n_buffer-1; j <= ny; j++) {
	//					for (k = n_buffer-1; k <= nz; k++) { 

	//						if ( FWS[icube][i][j][k] == 1 && (
	//							(FWS[icube][i+2][j][k] + FWS[icube][i-1][j][k] ==  20 ) |
	//							(FWS[icube][i][j+2][k] + FWS[icube][i][j-1][k] ==  20 ) |
	//							(FWS[icube][i][j][k+2] + FWS[icube][i][j][k-1] ==  20 ) 
	//							)) {
	//							FWS[icube][i][j][k] = 10;
	//							tempFWS[icube][i][j][k] = 100;
	//							end = end + 1;
	//						}

	//					}
	//				}
	//			}
	//			

	//		}    // ---- for (int iFWS = 1; iFWS < 100000; iFWS ++) ---- //




			// for (i = n_buffer-1; i <= nxx; i++) {
				// for (j = n_buffer-1; j <= nyy; j++) {
					// for (k = n_buffer-1; k <= nzz; k++) { 
					
						// if ( tempFWS[icube][i][j][k] == 100)

						// FWS[icube][i][j][k] = 100;

					// }
				// }
			// }
			

		}  // ---- if (csl[icube] == 0) ---- //

	}  
	



// -------------------------------------------- Hole treatment -------------------------------------------- //



// -------------------------------------------- Corner treatment -------------------------------------------- //

	delete[] tempFWS;

}                                          








