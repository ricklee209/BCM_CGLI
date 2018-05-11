



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
	
void BCM_Interface_EDGE
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

double (*send_data_curr_eq) = new double[5*NcubeX*n_buffer*n_buffer*4*Max_nei_eq+1],
double (*recv_data_curr_eq) = new double[5*NcubeX*n_buffer*n_buffer*4*Max_nei_eq+1],

double (*send_data_neig_eq) = new double[5*NcubeX*n_buffer*n_buffer*4*Max_nei_eq+1],
double (*recv_data_neig_eq) = new double[5*NcubeX*n_buffer*n_buffer*4*Max_nei_eq+1],

int (*Sdir_eq) = new int[Max_nei_eq+1],
int (*Rdir_eq) = new int[Max_nei_eq+1],

int (*ist_eq) = new int[Ncpu_eq],

int (*csl) = new int[Ncube],

int (*adj_number)[5][7] = new int[Ncube][5][7],

int (*adjX_eq) = new int[Ncube],
int (*adjY_eq) = new int[Ncube],
int (*adjZ_eq) = new int[Ncube],


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

	int icpu_neig_eq;

	int icube_send, icube_recv;

	int iL, iadj;

	int ii,jj,kk,_ii,_jj,_kk,ii_,jj_,kk_;

	int L1, L2, L3, L4, L5;

	int zone = NcubeX*n_buffer*n_buffer*4;
    
    int subzone = NcubeX*n_buffer*n_buffer;

	
	int index;

	MPI_Request (*Sreq_eq) = new MPI_Request[ncpu_eq];
	MPI_Request (*Rreq_eq) = new MPI_Request[ncpu_eq];

	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Status istat[8];
	MPI_Request requ1, reqps1;

	


// ================================================================================================== //
// ============================== # Package for sending in Current CPU ============================== //

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

		
#pragma omp parallel for private(iicube,icube,i,j,k,iL,k2,j2,i2)

		for (icube_send = 1; icube_send <= Ncube_Ncpu_eq[icpu_neig_eq]; icube_send++) {  

			iicube = ii-Ncube_Ncpu_eq[icpu_neig_eq]+icube_send;
			icube = Scube_Ncpu_eq[iicube];
			
			

// ------------------------------------- [Z+] direction ------------------------------------- //

			if (Sdir_eq[iicube] == 6) {

				for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

                            k2 = NcubeZ+k;
                
							iL = (icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
                            
						}
					}
				}
                
                for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

                            k2 = NcubeZ+k;
                
							iL = subzone+(icube_send-1)*zone+(i-nxx)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
                            
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  

                            k2 = NcubeZ+k;
                
							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
                            
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							
						}
					}
				}
                
                for (i = n_buffer;i <= nx; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = 0; k < n_buffer; k++) {  

                            k2 = NcubeZ+k;
                
							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+(j-nyy)*n_buffer+k;
                            
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


                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  
                
                            k2 = n_buffer+k;    
				
							iL = (icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							

						}
					}
				}
                
                for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  
               
                            k2 = n_buffer+k;    
				
							iL = subzone+(icube_send-1)*zone+(i-nxx)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							

						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  
               
                            k2 = n_buffer+k;    
				
							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j][k2][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j][k2][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j][k2][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j][k2][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j][k2][4];
							

						}
					}
				}
                
                for (i = n_buffer;i <= nx; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = 0; k < n_buffer; k++) {  
               
                            k2 = n_buffer+k;    
				
							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+(j-nyy)*n_buffer+k;
							
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

				for (i = 0; i < n_buffer; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+k;
							j2 = NcubeY+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}
                
                
               for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+(i-nxx)*n_buffer*NcubeZ+j*NcubeZ+k;
							j2 = NcubeY+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
							j2 = NcubeY+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
                    for (j = 0; j < n_buffer; j++) {
                        for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+(k-nzz);
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
                
                for (i = 0; i < n_buffer; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+k;
							j2 = n_buffer+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}
                
               for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+(i-nxx)*n_buffer*NcubeZ+j*NcubeZ+k;
							j2 = n_buffer+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
							j2 = n_buffer+j;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i][j2][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i][j2][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i][j2][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i][j2][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i][j2][k][4];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
                    for (j = 0; j < n_buffer; j++) {
                        for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+(k-nzz);
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
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							i2 = NcubeX+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+i*n_buffer*NcubeZ+(j-nyy)*NcubeZ+(k-n_buffer);
							i2 = NcubeX+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							i2 = NcubeX+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+(k-nzz);
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
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							i2 = n_buffer+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+i*n_buffer*NcubeZ+(j-nyy)*NcubeZ+(k-n_buffer);
							i2 = n_buffer+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							i2 = n_buffer+i;
							
							send_data_curr_eq[L1+iL] = U1_[icube][i2][j][k][0];
							send_data_curr_eq[L2+iL] = U1_[icube][i2][j][k][1];
							send_data_curr_eq[L3+iL] = U1_[icube][i2][j][k][2];
							send_data_curr_eq[L4+iL] = U1_[icube][i2][j][k][3];
							send_data_curr_eq[L5+iL] = U1_[icube][i2][j][k][4];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+(k-nzz);
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

	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		istart = ist_eq[icpu_neig_eq]*5*NcubeX*n_buffer*n_buffer*4;

		idest = neighbor_cpu_eq[icpu_neig_eq];

		icount = Ncube_Ncpu_eq[icpu_neig_eq]*5*NcubeX*n_buffer*n_buffer*4;

		itag = 0;

		MPI_Isend((void *)&send_data_curr_eq[istart], icount, MPI_DOUBLE, idest, itag, comm, &Sreq_eq[icpu_neig_eq]);

	}
	
	
	
	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		istart = ist_eq[icpu_neig_eq]*5*NcubeX*n_buffer*n_buffer*4;

		isrc = neighbor_cpu_eq[icpu_neig_eq];

		icount = Ncube_Ncpu_eq[icpu_neig_eq]*5*NcubeX*n_buffer*n_buffer*4;

		itag = 0;

		MPI_Irecv((void *)&recv_data_curr_eq[istart], icount, MPI_DOUBLE, isrc, itag, comm, &Rreq_eq[icpu_neig_eq]);

	}
	
    
// ---- the same size connection in X direction ---- //

#pragma omp parallel for private(ic0, ic2, i,j,k, i0,i2) schedule(dynamic)
// ============================================== //
	for (adjX = 1; adjX <= nadjX_eq; adjX++) {    //  
// ============================================== //

		ic0 = adjX_eq[adjX];
		ic2 = adj_number[ic0][1][2];
		
		for (i = 0; i < n_buffer; i++) {
			for (j = 0; j < n_buffer; j++) {
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
        
        for (i = 0; i < n_buffer; i++) {
			for (j = nyy; j < nyy+n_buffer+1; j++) {
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
        
        for (i = 0; i < n_buffer; i++) {
			for (j = n_buffer; j <= ny; j++) {
				for (k = 0; k < n_buffer; k++) {  

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
        
        for (i = 0; i < n_buffer; i++) {
			for (j = n_buffer; j <= ny; j++) {
				for (k = nzz; k < nzz+n_buffer+1; k++) {  

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





#pragma omp parallel for private(ic0, ic2, i,j,k, j0,j2) schedule(dynamic)

// ---- the same size connection in Y direction ---- //
// ============================================== //
	for (adjY = 1; adjY <= nadjY_eq; adjY++) {    //  
// ============================================== //

		ic0 = adjY_eq[adjY];
		ic2 = adj_number[ic0][1][4];
        
		for (i = 0; i < n_buffer; i++) {
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
        
        for (i = nxx; i < nxx+n_buffer+1; i++) {
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
        
        for (i = n_buffer; i <= nx; i++) {
			for (j = 0; j < n_buffer; j++) {
				for (k = 0; k < n_buffer; k++) {  

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
        
        for (i = n_buffer; i <= nx; i++) {
            for (j = 0; j < n_buffer; j++) {
                for (k = nzz; k < nzz+n_buffer+1; k++) {  

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






#pragma omp parallel for private(ic0, ic2, i,j,k, k0,k2) schedule(dynamic)

// ---- the same size connection in Z direction ---- //
// ============================================== //
	for (adjZ = 1; adjZ <= nadjZ_eq; adjZ++) {    //  
// ============================================== //

		ic0 = adjZ_eq[adjZ];
		ic2 = adj_number[ic0][1][6];

		for (i = 0; i < n_buffer; i++) {
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
        
        
        
		for (i = nxx; i < nxx+n_buffer+1; i++) {
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
        
        
		for (i = n_buffer; i <= nx; i++) {
			for (j = 0; j < n_buffer; j++) {
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
        
        
		for (i = n_buffer; i <= nx; i++) {
			for (j = nyy; j < nyy+n_buffer+1; j++) {
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





	
	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		MPI_Waitany(ncpu_eq, Rreq_eq, &index, istat);

	}
	

	
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

                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = (icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k0 = k;
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							
							
						}
					}
				}
                
                for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = subzone+(icube_send-1)*zone+(i-nxx)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
                            k0 = k;
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
                            k0 = k;
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = n_buffer;i <= nx; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+(j-nyy)*n_buffer+k;
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

                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  
                
							iL = (icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k0 = n_buffer+NcubeZ+k;
							
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							

						}
					}
				}
                
                for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  
               
							iL = subzone+(icube_send-1)*zone+(i-nxx)*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							k0 = n_buffer+NcubeZ+k;
							
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							

						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  
               
							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
							k0 = n_buffer+NcubeZ+k;
							
							
							U1_[icube][i][j][k0][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j][k0][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j][k0][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j][k0][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j][k0][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = n_buffer;i <= nx; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = 0; k < n_buffer; k++) {  
               
							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+(j-nyy)*n_buffer+k;
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

                
                
				for (i = 0; i < n_buffer; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+k;
							j0 = j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
                            
						}
					}
				}
                
               for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+(i-nxx)*n_buffer*NcubeZ+j*NcubeZ+k;
							j0 = j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
							j0 = j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
                    for (j = 0; j < n_buffer; j++) {
                        for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+(k-nzz);
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

                
                for (i = 0; i < n_buffer; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+k;
							j0 = n_buffer+NcubeY+j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
               for (i = nxx; i < nxx+n_buffer+1; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+(i-nxx)*n_buffer*NcubeZ+j*NcubeZ+k;
							j0 = n_buffer+NcubeY+j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = n_buffer; i <= nx; i++) {
					for (j = 0; j < n_buffer; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+k;
							j0 = n_buffer+NcubeY+j;
							
							U1_[icube][i][j0][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i][j0][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i][j0][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i][j0][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i][j0][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
               for (i = n_buffer; i <= nx; i++) {
                    for (j = 0; j < n_buffer; j++) {
                        for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+(i-n_buffer)*n_buffer*n_buffer+j*n_buffer+(k-nzz);
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
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							i0 = i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+i*n_buffer*NcubeZ+(j-nyy)*NcubeZ+(k-n_buffer);
							i0 = i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							i0 = i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+(k-nzz);
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
					for (j = 0; j < n_buffer; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = (icube_send-1)*zone+i*n_buffer*NcubeZ+j*NcubeZ+(k-n_buffer);
							i0 = n_buffer+NcubeX+i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = nyy; j < nyy+n_buffer+1; j++) {
						for (k = n_buffer; k <= nz; k++) {  

							iL = subzone+(icube_send-1)*zone+i*n_buffer*NcubeZ+(j-nyy)*NcubeZ+(k-n_buffer);
							i0 = n_buffer+NcubeX+i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = 0; k < n_buffer; k++) {  

							iL = 2*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+k;
							i0 = n_buffer+NcubeX+i;
							
							U1_[icube][i0][j][k][0] = recv_data_curr_eq[L1+iL];
							U1_[icube][i0][j][k][1] = recv_data_curr_eq[L2+iL];
							U1_[icube][i0][j][k][2] = recv_data_curr_eq[L3+iL];
							U1_[icube][i0][j][k][3] = recv_data_curr_eq[L4+iL];
							U1_[icube][i0][j][k][4] = recv_data_curr_eq[L5+iL];
							
						}
					}
				}
                
                for (i = 0; i < n_buffer; i++) {
					for (j = n_buffer; j <= ny; j++) {
						for (k = nzz; k < nzz+n_buffer+1; k++) {  

							iL = 3*subzone+(icube_send-1)*zone+i*NcubeY*n_buffer+(j-n_buffer)*n_buffer+(k-nzz);
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


	for (icpu_neig_eq = 0;  icpu_neig_eq < ncpu_eq; icpu_neig_eq++) {

		MPI_Waitany(ncpu_eq, Sreq_eq, &index, istat);

	}

	
}                                          








