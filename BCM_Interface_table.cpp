



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <string.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;    
extern int N_wallcube;    
extern int MPI_Nadj;  

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
)

{

#include "BCM.h"
#include "prm.h"

	int count_index = 0;

	int acube;

	int adjX,adjY,adjZ;
	int ic0,ic1,ic2,ic3,ic4,ic5,ic6;
	int i0,i1,i2,i3,i4,i5,i6;
	int j0,j1,j2,j3,j4,j5,j6;
	int k0,k1,k2,k3,k4,k5,k6;



	*nadjX_eq = *nadjX_bs_plus = *nadjX_sb_plus = *nadjX_bs_minus = *nadjX_sb_minus = 0;
	*nadjY_eq = *nadjY_bs_plus = *nadjY_sb_plus = *nadjY_bs_minus = *nadjY_sb_minus = 0;
	*nadjZ_eq = *nadjZ_bs_plus = *nadjZ_sb_plus = *nadjZ_bs_minus = *nadjZ_sb_minus = 0;


// --------------------------------------------- for each CPU --------------------------------------------- //

// ============================================ //
	for (icube = 1; icube < ncube; icube++) {   //
// ============================================ //

		if (adj_number[icube][1][2] > 0 |
			adj_number[icube][2][2] > 0 |
			adj_number[icube][3][2] > 0 |
			adj_number[icube][4][2] > 0) {
				
				acube = max(max(max(adj_number[icube][1][2],adj_number[icube][2][2]),adj_number[icube][3][2]),adj_number[icube][4][2]);
				
				// ------------------------------------------------------------- //
				if (csl[icube] == csl[acube]) {

					*nadjX_eq = *nadjX_eq+1;
					adjX_eq[*nadjX_eq] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				else if (csl[icube] < csl[acube]) {

					*nadjX_sb_plus = *nadjX_sb_plus+1;
					adjX_sb_plus[*nadjX_sb_plus] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				else if (csl[icube] > csl[acube]) {

					*nadjX_bs_plus = *nadjX_bs_plus+1;
					adjX_bs_plus[*nadjX_bs_plus] = icube;

				}
				// ------------------------------------------------------------- //
				
		}
	
		if (adj_number[icube][1][4] > 0 |
			adj_number[icube][2][4] > 0 |
			adj_number[icube][3][4] > 0 |
			adj_number[icube][4][4] > 0) {

				acube = max(max(max(adj_number[icube][1][4],adj_number[icube][2][4]),adj_number[icube][3][4]),adj_number[icube][4][4]);

				// ------------------------------------------------------------- //
				if (csl[icube] == csl[acube]) {

					*nadjY_eq = *nadjY_eq+1;
					adjY_eq[*nadjY_eq] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				else if (csl[icube] < csl[acube]){

					*nadjY_sb_plus = *nadjY_sb_plus+1;
					adjY_sb_plus[*nadjY_sb_plus] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				else if (csl[icube] > csl[acube]){

					*nadjY_bs_plus = *nadjY_bs_plus+1;
					adjY_bs_plus[*nadjY_bs_plus] = icube;

				}
				// ------------------------------------------------------------- //

		}




		if (adj_number[icube][1][6] > 0 |
			adj_number[icube][2][6] > 0 |
			adj_number[icube][3][6] > 0 |
			adj_number[icube][4][6] > 0) {

				acube = max(max(max(adj_number[icube][1][6],adj_number[icube][2][6]),adj_number[icube][3][6]),adj_number[icube][4][6]);

				// ------------------------------------------------------------- //
				if (csl[icube] == csl[acube]) {

					*nadjZ_eq = *nadjZ_eq+1;
					adjZ_eq[*nadjZ_eq] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				else if (csl[icube] < csl[acube]){

					*nadjZ_sb_plus = *nadjZ_sb_plus+1;
					adjZ_sb_plus[*nadjZ_sb_plus] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				else if (csl[icube] > csl[acube]){

					*nadjZ_bs_plus = *nadjZ_bs_plus+1;
					adjZ_bs_plus[*nadjZ_bs_plus] = icube;

				}
				// ------------------------------------------------------------- //

		}






		if (adj_number[icube][1][1] > 0 |
			adj_number[icube][2][1] > 0 |
			adj_number[icube][3][1] > 0 |
			adj_number[icube][4][1] > 0) {
				
				acube = max(max(max(adj_number[icube][1][1],adj_number[icube][2][1]),adj_number[icube][3][1]),adj_number[icube][4][1]);
				
				
				// ------------------------------------------------------------- //
				if (csl[icube] < csl[acube]) {

					*nadjX_sb_minus = *nadjX_sb_minus+1;
					adjX_sb_minus[*nadjX_sb_minus] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				if (csl[icube] > csl[acube]) {

					*nadjX_bs_minus = *nadjX_bs_minus+1;
					adjX_bs_minus[*nadjX_bs_minus] = icube;

				}
				// ------------------------------------------------------------- //
				
		}
	
		if (adj_number[icube][1][3] > 0 |
			adj_number[icube][2][3] > 0 |
			adj_number[icube][3][3] > 0 |
			adj_number[icube][4][3] > 0) {

				acube = max(max(max(adj_number[icube][1][3],adj_number[icube][2][3]),adj_number[icube][3][3]),adj_number[icube][4][3]);

				
				// ------------------------------------------------------------- //
				if (csl[icube] < csl[acube]){

					*nadjY_sb_minus = *nadjY_sb_minus+1;
					adjY_sb_minus[*nadjY_sb_minus] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				if (csl[icube] > csl[acube]) {

					*nadjY_bs_minus = *nadjY_bs_minus+1;
					adjY_bs_minus[*nadjY_bs_minus] = icube;

				}
				// ------------------------------------------------------------- //

		}




		if (adj_number[icube][1][5] > 0 |
			adj_number[icube][2][5] > 0 |
			adj_number[icube][3][5] > 0 |
			adj_number[icube][4][5] > 0) {

				acube = max(max(max(adj_number[icube][1][5],adj_number[icube][2][5]),adj_number[icube][3][5]),adj_number[icube][4][5]);

				// ------------------------------------------------------------- //
				if (csl[icube] < csl[acube]){

					*nadjZ_sb_minus = *nadjZ_sb_minus+1;
					adjZ_sb_minus[*nadjZ_sb_minus] = icube;

				}
				// ------------------------------------------------------------- //

				// ------------------------------------------------------------- //
				if (csl[icube] > csl[acube]){

					*nadjZ_bs_minus = *nadjZ_bs_minus+1;
					adjZ_bs_minus[*nadjZ_bs_minus] = icube;

				}
				// ------------------------------------------------------------- //

		}




// ============================================ //
	}                                           //
// ============================================ //
	
}                                          








