



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int Ncube;    
extern int N_wallcube;    

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
)

{

#include "BCM.h"
#include "prm.h"

	int Ncell = NcubeX*NcubeY*NcubeZ;

	int *jkl_switch;
	jkl_switch = new int [Ncell+1];
	
	int ii, i1, i2, pp, qq, n_change, flag, ch_pos, count_index;


	char str[1024];
	int BCM_ncube, BCM_n_wallcube, BCM_n;

	FILE *fptr;
	char file_name[100];
	sprintf(file_name,"BCMgrid""%0.5d"".dat",myid);   
	fptr = fopen(file_name,"r"); 

	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);

// ======================================================= //
	fscanf(fptr,"%d\n",&BCM_ncube);    //  how many cubes  //
// ======================================================= //

	fscanf(fptr,"%[^\n]\n",str);

// ======================================================================================================= //
	for (icube = 1; icube < ncube; icube++) {                                                              //
		                                                                                                   //
		fscanf(fptr,"%lf\t%lf\t%lf\t%lf\n",&cube_size[icube],&Xcube[icube],&Ycube[icube],&Zcube[icube]);   //
		                                                                                                   //
	}                                                                                                      //
// ======================================================================================================= //

	fscanf(fptr,"%[^\n]\n",str);

// =================================================================================== //
	for (icube = 1; icube < ncube; icube++) {                                          //
		for (int direction_index = 1; direction_index <= 6;  direction_index++) {      //
			                                                                           //
			fscanf(fptr,"%d\t%d\t%d\t%d\n",&adj_number[icube][1][direction_index],     //
				                           &adj_number[icube][2][direction_index],     //
										   &adj_number[icube][3][direction_index],     //
										   &adj_number[icube][4][direction_index]);    //
		}                                                                              //
	}                                                                                  //
// =================================================================================== //

	fscanf(fptr,"%[^\n]\n",str);

// ============================================================= //
	fscanf(fptr,"%d\n",&BCM_n);    //  how many cells in a cube  //
// ============================================================= //

	fscanf(fptr,"%[^\n]\n",str);

// ================================================================ //
	fscanf(fptr,"%d\n",&BCM_n_wallcube);    //  how many wallcubes  //
// ================================================================ //

	fscanf(fptr,"%[^\n]\n",str);

// =============================================================== //
	for (iwallcube = 1; iwallcube < n_wallcube; iwallcube++) {     //
		                                                           //
		fscanf(fptr,"%d\n",&wallcube[iwallcube]);                  //
		                                                           //
	}                                                              //
// =============================================================== //

	fscanf(fptr,"%[^\n]\n",str);

	 
// ========================================================= //
	for (icube = 1; icube < ncube; icube++) { 

		dx = cube_size[icube]/NcubeX;
		dy = cube_size[icube]/NcubeY;
		dz = cube_size[icube]/NcubeZ;

		#pragma omp parallel for 
		for (i = 0; i <= nxxx; i++) {

			Xcnt[icube][i] = Xcube[icube]+(i-n_buffer+0.5)*dx;

		}

		#pragma omp parallel for
		for (j = 0; j <= nyyy; j++) {

			Ycnt[icube][j] = Ycube[icube]+(j-n_buffer+0.5)*dy;

		}

		#pragma omp parallel for
		for (k = 0; k <= nzzz; k++) {

			Zcnt[icube][k] = Zcube[icube]+(k-n_buffer+0.5)*dz;

		}
		
	}
// ========================================================= //



double size_min = MAX;

// ===================================================================================== //
	for (icube = 1; icube < ncube; icube++) {                                            //
		                                                                                 //
		size_min = min(size_min,cube_size[icube]);                                       //
		                                                                                 //
	}																					 //
	                                                                                     //
	for (icube = 1; icube < ncube; icube++) {                                            //
																					     //
	  csl[icube] = static_cast<int>(0.5+log( cube_size[icube]/size_min ) / log(2.0));    //
																						 //
	}																					 //
		                                                                                 //
// ===================================================================================== //

	for (icube = 1; icube < ncube; icube++) {    
#pragma omp parallel for private(j,k)
		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {
				for (k = 0; k <= nzzz; k++) {  

					FWS[icube][i][j][k] = IFLUID;   // ---- all fluid ---- //

				}
			}
		}
	}
	

	


// ---- run length ---- //
// =============================================================== //
	for (iwallcube = 1; iwallcube < n_wallcube; iwallcube++) {     //
// =============================================================== //

		icube = wallcube[iwallcube];

		fscanf(fptr,"%d\n",&n_change);    // ---- how many cells needed to be changed the flag ---- //
		
		fscanf(fptr,"%d\n",&flag);    // ---- the first cell's flag ---- //

		if (n_change == 0) {    // ---- in this cube, all flags is the same ---- //

#pragma omp parallel for private(j,k)
			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) {  

						FWS[icube][i][j][k] = flag; 

					}    // ---- nx ---- //
				}    // ---- ny ---- //
			}    // ---- nz ---- //
#pragma omp barrier
		} 

		else {

			jkl_switch[1] = flag;    // ---- the initial value of the flag ---- //

			fscanf(fptr,"%d",&ch_pos);    // ---- the position of the cell which changes flag ---- //
			
			count_index = 0;

			for (int iflag = 2; iflag <= Ncell; iflag++) {

				if (iflag == ch_pos) {

					count_index = count_index+1;

					jkl_switch[iflag] = max(abs(jkl_switch[iflag-1]-1),0);

					if (count_index < n_change) fscanf(fptr,"%d\n",&ch_pos);    

				}

				else jkl_switch[iflag] = jkl_switch[iflag-1];

			}    // ---- iflag ---- //

#pragma omp parallel for private(j,i,pp,qq,i1,i2,ii)
			for (k = 1; k <= NcubeZ; k++) {  
				for (j = 1; j <= NcubeY; j++) {
					for (i = 1; i <= NcubeX; i++) {

						pp = j%2;
						qq = k%2;
						i1 = NcubeX*(j-1)+pp*i+(1-pp)*(NcubeX+1-i);
						i2 = NcubeX*(NcubeY-j)+pp*(NcubeX+1-i)+(1-pp)*i;
						ii = qq*i1+(1-qq)*i2+NcubeX*NcubeY*(k-1);

						FWS[icube][i+n_buffer-1][j+n_buffer-1][k+n_buffer-1] = jkl_switch[ii];

					}
				}
			}
#pragma omp barrier

		}    // ---- else ---- //


// =============================================================== //
	}    // ---- nwallcube ---- //                                 //
// =============================================================== //

	
	for (icube = 1; icube < ncube; icube++) {    
#pragma omp parallel for private(j,k)
		
		for (i = n_buffer; i <= nx; i++) {
			for (j = n_buffer; j <= ny; j++) {
				for (k = n_buffer; k <= nz; k++) { 
		
					if (FWS[icube][i][j][k] == 0)
						FWS[icube][i][j][k] = -1;  

				}
			}
		}
	}

	
//	int (*tempFWS)[X_size][Y_size][Z_size] = new int[ncube][X_size][Y_size][Z_size]; 
//
//
//	for (icube = 1; icube < ncube; icube++) {  
//
//		if (csl[icube] == 0) {
//#pragma omp parallel for private(j,k)
//
//			for (i = n_buffer; i <= nx; i++) {
//				for (j = n_buffer; j <= ny; j++) {
//					for (k = n_buffer; k <= nz; k++) { 
//						
//						if (FWS[icube][i][j][k] == -1 & 
//							(FWS[icube][i+1][j  ][k  ] == IFLUID |
//							FWS[icube][i-1][j  ][k  ] == IFLUID |
//							FWS[icube][i  ][j-1][k+1] == IFLUID |
//							FWS[icube][i+1][j+1][k  ] == IFLUID |
//							FWS[icube][i  ][j  ][k-1] == IFLUID |
//							FWS[icube][i  ][j  ][k+1] == IFLUID )
//							)
//							
//
//							tempFWS[icube][i  ][j  ][k  ] = -2;  
//						
//					}
//				}
//			}
//
//		}
//	}
//
//	
//	for (icube = 1; icube < ncube; icube++) {  
//
//		
//		if (csl[icube] == 0)
//#pragma omp parallel for private(j,k)
//			
//			for (i = n_buffer; i <= nx; i++) {
//				for (j = n_buffer; j <= ny; j++) {
//					for (k = n_buffer; k <= nz; k++) { 
//					
//			
//						if (tempFWS[icube][i][j][k] == -2) {
//							FWS[icube][i  ][j  ][k  ] = IGHOST; 
//						
//
//						}
//					}
//				}
//			}
//
//	}    // ---- for (icube = 1; icube < ncube; icube++) ---- //
//	
	



	delete[] jkl_switch;

	fclose(fptr);

}