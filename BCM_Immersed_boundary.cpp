



#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <mpi.h>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>

using namespace std;


unsigned int iswap(unsigned int i) {    
  union
  {
    unsigned int i;
    unsigned char b[4];
  } dat1, dat2;
  
  dat1.i = i;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.i;
} 


float fswap(float f) {    
  union
  {
    float f;
    unsigned char b[8];
  } dat1, dat2;
  
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  
  return dat2.f;
  
}


double dswap(double f) {    
  union
  {
    double f;
    unsigned char b[8];
  } dat1, dat2;
  
  dat1.f = f;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.f;
}



struct Prop
{

	int Node_index;

	double Vx; 
	double Vy;
	double Vz;

};

struct Array
{
	Prop *content;
	unsigned int size;
};


int compareX( const void *a, const void *b )
{

	return ( (Prop*)a)->Vx > ( (Prop*)b)->Vx ?1:-1;

};

int compareY( const void *a, const void *b )
{

	return ( (Prop*)a)->Vy > ( (Prop*)b)->Vy ?1:-1;

};


int compareZ( const void *a, const void *b )
{

	return ( (Prop*)a)->Vz > ( (Prop*)b)->Vz ?1:-1;

};




#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"
#include "BCM_Immersed_boundary.h"

extern int Ncube;    
extern int N_wallcube;    

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
)

{

#include "BCM.h"
#include "prm.h"

	char file_name[100];

	double V[64]; // Vandermonde matrix //

	double sur_x[9];
	double sur_y[9];
	double sur_z[9];

	double weight[8];


	int itri, Ntri, count_index, Scount_index, Ecount_index;

	int N_share_node, i_share_node;

	double temp;

	double Nx, Ny, Nz, V1x, V1y,V1z, V2x, V2y, V2z, V3x, V3y, V3z;

	int pitri;

	double pV1x, pV1y, pV1z, pV2x, pV2y, pV2z, pV3x, pV3y, pV3z;

	char str[1024];

	double size_min;

	char str_endloop[1024] = "endloop";
	char str_normal[1024] = "normal";


	FILE *fptr;
	FILE *fptr_plus,*fptr_minus;
	fptr = fopen("BCM_STL.stl","r");
	//fptr = fopen("BCM_STL_bin.stl","rb");


	count_index = 0;

	// ---- calculate how many vertices ---- //
	// ========================================= //

	
	Ntri = 0;

	while (!feof(fptr)) {
		
		fscanf(fptr,"%s",str); 
		if(strstr(str_endloop,str) != NULL) Ntri = Ntri+1;
		
		fgets(str, 1024, fptr);

	}

	rewind(fptr);    // ---- back to the beginning of the file ---- //

	// ========================================= //

	int Ntemp, NNtemp, N_line = 15;

	double *tri_ = new double[Ntri*N_line+1];  
	
	itri = 0;


	// ---- ( Nx, Ny, Nz, V1x, V1y,V1z, V2x, V2y, V2z, V3x, V3y, V3z) * nver ---- //
	
	while (!feof(fptr)) {

		fscanf(fptr,"%s",str); 

		if(strstr(str_normal,str) != NULL) {

			itri = itri + 1;

			fscanf(fptr,"%lf%lf%lf%s%s%s%lf%lf%lf%s%lf%lf%lf%s%lf%lf%lf",
				&Nx,&Ny,&Nz,str,str,str,&V1x,&V1y,&V1z,str,&V2x,&V2y,&V2z,str,&V3x,&V3y,&V3z);

			// ==== ( Nx, Ny, Nz, V1x, V1y,V1z, V2x, V2y, V2z, V3x, V3y, V3z) ==== //

			Ntemp = (itri-1)*N_line;

			tri_[Ntemp+1] = Nx;
			tri_[Ntemp+2] = Ny;
			tri_[Ntemp+3] = Nz;

			tri_[Ntemp+4] = 0.;  // ---- switch turn on ---- //

			tri_[Ntemp+5] = V1x;
			tri_[Ntemp+6] = V1y;
			tri_[Ntemp+7] = V1z;

			tri_[Ntemp+8] = 0.;  // ---- switch turn on ---- //

			tri_[Ntemp+ 9] = V2x;
			tri_[Ntemp+10] = V2y;
			tri_[Ntemp+11] = V2z;

			tri_[Ntemp+12] = 0.;  // ---- switch turn on ---- //

			tri_[Ntemp+13] = V3x;
			tri_[Ntemp+14] = V3y;
			tri_[Ntemp+15] = V3z;

		}

	}

	fclose(fptr);

	
	
	 //Ntri = 0;

	 //char title[80];

	 //fread(title, 80, 1, fptr);
	 //fread(&Ntri, 4, 1, fptr);

	 //printf("triangle number = %d\n\n",Ntri);

	 //// ========================================= //

	 //int Ntemp, NNtemp, N_line = 15;

	 //float v[12];

	 //double *tri_ = new double[Ntri*N_line+1];  

	 //itri = 0;
	
	 //float enidan_temp;


	 //// ---- ( Nx, Ny, Nz, V1x, V1y,V1z, V2x, V2y, V2z, V3x, V3y, V3z) * nver ---- //

	 //unsigned short uint16; 

	 //// Every Face is 50 Bytes: Normal(3*float), Vertices(9*float), 2 Bytes Spacer

	 //for (itri = 0; itri < Ntri; ++itri) {

		// for (i = 0; i < 12; ++i) {

		//	 fread((void*)&enidan_temp, sizeof(float), 1, fptr);

		//	 v[i] = fswap(enidan_temp); 
		//	
		// }

		// Ntemp = (itri)*N_line;

		// tri_[Ntemp+1] = v[0];
		// tri_[Ntemp+2] = v[1];
		// tri_[Ntemp+3] = v[2];

		// tri_[Ntemp+5] = v[3];
		// tri_[Ntemp+6] = v[4];
		// tri_[Ntemp+7] = v[5];

		// tri_[Ntemp+9] = v[6];
		// tri_[Ntemp+10] = v[7];
		// tri_[Ntemp+11] = v[8];

		// tri_[Ntemp+13] = v[9];
		// tri_[Ntemp+14] = v[10];
		// tri_[Ntemp+15] = v[11];

		// fread((void*)&uint16, sizeof(unsigned short), 1, fptr); // spacer between successive faces

	 //}

	 //fclose(fptr);



	// ----------------------------------------------------------- //
	// -------------- tri_num = 0 inside the domain -------------- //

	double boxcenter[3];
	double boxhalfsize[3];
	double triverts[3][3];

	int *tri_num = new int[Ntri+1];

	#pragma omp parallel for 
	for (itri = 1; itri <= Ntri; itri++) {

		tri_num[itri] = -1;

	}


	//#pragma omp parallel for private

	count_index = 0;

	for (itri = 1; itri <= Ntri; itri++) {

		Ntemp = (itri-1)*N_line;

		for (icube = 1; icube < ncube; icube++) {  

			dx = dy = dz = cube_size[icube]/NcubeX;
			
			boxhalfsize[0] = boxhalfsize[1] = boxhalfsize[2] = 0.5*NcubeX*dx+2*dx;

			boxcenter[0] = Xcube[icube]+(NcubeX/2)*dx;
			boxcenter[1] = Ycube[icube]+(NcubeY/2)*dy;
			boxcenter[2] = Zcube[icube]+(NcubeZ/2)*dz;


			triverts[0][0] = tri_[Ntemp+5];
			triverts[0][1] = tri_[Ntemp+6];
			triverts[0][2] = tri_[Ntemp+7];

			triverts[1][0] = tri_[Ntemp+9];
			triverts[1][1] = tri_[Ntemp+10];
			triverts[1][2] = tri_[Ntemp+11];

			triverts[2][0] = tri_[Ntemp+13];
			triverts[2][1] = tri_[Ntemp+14];
			triverts[2][2] = tri_[Ntemp+15];

			if (triBoxOverlap(boxcenter,boxhalfsize,triverts) == 1) {
				
				tri_num[itri] = 0;    /* box and triangle overlaps */

				count_index = count_index+1;

			}
		}

	}


	int *cube_trinum = new int[ncube];    // ---- how many triangles in each cube ---- //
	int *tri_table = new int[count_index+1];



	//#pragma omp barrier

	
	// -------------- tri_num = 0 inside the domain -------------- //
	// ----------------------------------------------------------- //




	// ------------------------------------------------------------------------------- //
	// ------------------ calculate how many triangles in this rank ------------------ //

	count_index = 0;

#pragma omp parallel for reduction(+:count_index)

	for (itri = 1; itri <= Ntri; itri++) {

		if (tri_num[itri] == 0) count_index = count_index+1;

	}

	// ------------------ calculate how many triangles in this rank ------------------ //
	// ------------------------------------------------------------------------------- //


	
	// ----------------------------------------------------------------------------- //
	// ------------------ remove the triangles outside the domain ------------------ //

	double *tri = new double[count_index*N_line+1];  

	count_index = 0;

	for (itri = 1; itri <= Ntri; itri++) {

		NNtemp = (itri-1)*N_line;

		if (tri_num[itri] == 0) {

			count_index = count_index+1;
			Ntemp = (count_index-1)*N_line;

			tri[Ntemp+1] = tri_[NNtemp+1];
			tri[Ntemp+2] = tri_[NNtemp+2];
			tri[Ntemp+3] = tri_[NNtemp+3];

			tri[Ntemp+4] = 0.;  // ---- switch turn on ---- //

			tri[Ntemp+5] = tri_[NNtemp+5];
			tri[Ntemp+6] = tri_[NNtemp+6];
			tri[Ntemp+7] = tri_[NNtemp+7];

			tri[Ntemp+8] = 0.;  // ---- switch turn on ---- //

			tri[Ntemp+9] = tri_[NNtemp+9];
			tri[Ntemp+10] = tri_[NNtemp+10];
			tri[Ntemp+11] = tri_[NNtemp+11];

			tri[Ntemp+12] = 0.;  // ---- switch turn on ---- //

			tri[Ntemp+13] = tri_[NNtemp+13];
			tri[Ntemp+14] = tri_[NNtemp+14];
			tri[Ntemp+15] = tri_[NNtemp+15];


		}

	}
	

	// ------------------ remove the triangles outside the domain ------------------ //
	// ----------------------------------------------------------------------------- //
	



	int Np = 3*count_index;    // ---- Each triangle has 3 vertexies ---- //

	Ntri = count_index;    // ---- triangle numbers inside the domain ---- //

	int itri_index = 0;
	int *pNode = new int[Np+1];    // ---- nodes serial number (share the same vertex) ---- //
	int *pNode_inf = new int[Ntri+1];    // ---- nodes information (end) ---- //

	int index = 0;
	int *ptri_number = new int[Ntri+1];    // ---- how many triangles share the same vertex ---- //

	double *pvertex = new double[Ntri*3+1];    // ---- vertex serial number ---- //

	double tri0,tri1,tri2,tri3;

	int N0,N1,N2,N3;


	// --------------------------------------------------------------------------------- //
	// ----------------------- detect the triangle in which cube ----------------------- //

	
	count_index = 0;

	for (icube = 1; icube < ncube; icube++) {  
		
		dx = dy = dz = cube_size[icube]/NcubeX;

		boxhalfsize[0] = boxhalfsize[1] = boxhalfsize[2] = 0.5*NcubeX*dx+2*dx;
		
		boxcenter[0] = Xcube[icube]+(NcubeX/2)*dx;
		boxcenter[1] = Ycube[icube]+(NcubeY/2)*dy;
		boxcenter[2] = Zcube[icube]+(NcubeZ/2)*dz;

		for (itri = 1; itri <= Ntri; itri++) {

			Ntemp = (itri-1)*N_line;

			triverts[0][0] = tri[Ntemp+5];
			triverts[0][1] = tri[Ntemp+6];
			triverts[0][2] = tri[Ntemp+7];

			triverts[1][0] = tri[Ntemp+9];
			triverts[1][1] = tri[Ntemp+10];
			triverts[1][2] = tri[Ntemp+11];
			
			triverts[2][0] = tri[Ntemp+13];
			triverts[2][1] = tri[Ntemp+14];
			triverts[2][2] = tri[Ntemp+15];

			
			if (triBoxOverlap(boxcenter,boxhalfsize,triverts) == 1) {

				count_index = count_index+1;

				tri_table[count_index] = itri;

			}  // ---- if (triBoxOverlap(boxcenter,boxhalfsize,triverts) == 1) ---- //

		}

		cube_trinum[icube] = count_index; 

	}


	delete[] tri_num;
	delete[] tri_;

	
	// ----------------------- detect the triangle in which cube ----------------------- //
	// --------------------------------------------------------------------------------- //


	int ii,jj,kk;

	int Ngc = 0;


	for (icube = 1; icube < ncube; icube++) {  

		if (csl[icube] == 0) {

			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
					for (k = 0; k <= nzzz; k++) { 

						if (FWS[icube][i][j][k] == IGHOST)
							Ngc = Ngc+1;

					}
				}
			}

		}
	}


	// ======================================== //




	int igc, jgc, kgc;

	double dmin = MAX;
	double dis;


	int *GC = new int[Ngc*4+1];
	double *GCcnt = new double[Ngc*3+1];

	// ---- (icube, i, j, k, xcnt, ycnt ,zctn) ---- //


	//// ================  calculate the ghost cells ================ //


	int icount;

	Ngc = 0;
	for (icube = 1; icube < ncube; icube++) {  

		if (csl[icube] == 0) {

			//for (i = 0; i <= nxxx; i++) {
				//for (j = 0; j <= nyy; j++) {
					//for (k = 0; k <= nzz; k++) { 

			for (i = n_buffer; i <= nx; i++) {
				for (j = n_buffer; j <= ny; j++) {
					for (k = n_buffer; k <= nz; k++) { 


						if (FWS[icube][i][j][k] == IGHOST) {

							Ngc = Ngc+1;

							Ntemp = (Ngc-1)*4;

							GC[Ntemp+1] = icube;

							GC[Ntemp+2] = i;
							GC[Ntemp+3] = j;
							GC[Ntemp+4] = k;

							Ntemp = (Ngc-1)*3;

							GCcnt[Ntemp+1] = Xcnt[icube][i];
							GCcnt[Ntemp+2] = Ycnt[icube][j];
							GCcnt[Ntemp+3] = Zcnt[icube][k];

							//FWS[icube][i][j][k] = IFLUID;

						}    // ---- if (FWS[icube][i][j][k] == IGHOST) ---- //

					}    // ---- for (i = 0; i < NcubeX; i++) ---- //
				}    // ---- for (j = 0; j <= nyyy; j++) ---- //
			}    // ---- for (k = 0; k <= nzzz; k++) ---- //

		}    // ---- if (csl[icube] == 0) ---- //

	}



	// ================  calculate the ghost cells ================ //



	// ============ immersed boundary file ============= //

	//	char file_name[100];
	sprintf(file_name,"IBM_plus""%0.5d"".dat",myid);    
	fptr_plus = fopen(file_name,"w"); 

	sprintf(file_name,"IBM_minus""%0.5d"".dat",myid);    
	fptr_minus = fopen(file_name,"w"); 


	// ============ immersed boundary file ============ //




	//fprintf(fptr,"icube i j k and weight of surrounding points\n");


	double orig[3];
	double dir[3];
	double vert0[3];
	double vert1[3];
	double vert2[3];
	double nuvec[3];


	double tmin, ttmin;

	double Ndis;
	double Ndis_min = MAX;

	double dotp;
	double BIx, BIy, BIz, IPx, IPy, IPz;
	double BIp0,BIp1,BIp2;
	double dir0, dir1, dir2;

	double xcnt, ycnt, zcnt;

	int pNgc = 0;
	int	itemp = 0;
	int	pBI;


	int iflag, iflag_total;
	int ilarge = 0;
	int nlarge[300];

	int iNgc_plus = 0;
	int iNgc_minus = 0;

	double s1, s2, s3;    // ---- surrounding points Sx Sy Sz ---- //

	double w1,w2,w3,w4,w5,w6,w7,w8;

	double vert0_iflag_total0_x,vert0_iflag_total0_y,vert0_iflag_total0_z;

    
	for (igc = 1; igc <= Ngc; igc++) {

		if(Ntri == 0) continue;    // ---- cube contains ghost cells but no triangles ---- //

		Ntemp = (igc-1)*3;

		xcnt = GCcnt[Ntemp+1];
		ycnt = GCcnt[Ntemp+2];
		zcnt = GCcnt[Ntemp+3];

		NNtemp = (igc-1)*4;

		icube = GC[NNtemp+1];
		i = GC[NNtemp+2];
		j = GC[NNtemp+3];
		k = GC[NNtemp+4];

		dx = cube_size[icube]/NcubeX;
		dy = cube_size[icube]/NcubeY;
		dz = cube_size[icube]/NcubeZ;

		iflag_total = 0;
		tmin = MAX;

		Ndis_min = MAX;

		itemp = -1;
		ilarge = 0;


		for (itri = 1; itri <= cube_trinum[icube]; itri++) {

			Ntemp = (tri_table[itri]-1)*N_line;

			BI_detect(myid, Ntemp, N_line, &Ndis, &Ndis_min, &iflag_total, &tmin, &dotp, &dir0, &dir1, &dir2, pNode_inf, ptri_number, pNode, 
				tri, vert0, vert1, vert2, nuvec, dir, orig, xcnt, ycnt, zcnt, &itemp, &ilarge, nlarge);

		}

		// ---- if (dotp > 0) the point is ouside the object ---- //


		if (iflag_total > 0) {

			if (dotp == 0) {

				BIx = IPx = orig[0];
				BIy = IPy = orig[1];
				BIz = IPz = orig[2];

			}  // ---- if (dotp == 0) ---- //

			else if (dotp > 0){

				BIx = orig[0] - tmin*dir0;
				BIy = orig[1] - tmin*dir1;
				BIz = orig[2] - tmin*dir2;

				IPx = orig[0] + tmin*dir0;
				IPy = orig[1] + tmin*dir1;
				IPz = orig[2] + tmin*dir2;

			}  // ---- else if (dotp > 0) ---- //

			else {

				BIx = orig[0] + tmin*dir0;
				BIy = orig[1] + tmin*dir1;
				BIz = orig[2] + tmin*dir2;

				IPx = orig[0] - tmin*dir0;
				IPy = orig[1] - tmin*dir1;
				IPz = orig[2] - tmin*dir2;

			}

			dis = (BIx-orig[0])*(BIx-orig[0])+(BIy-orig[1])*(BIy-orig[1])+(BIz-orig[2])*(BIz-orig[2]);
			
		}  // ---- if (iflag_total > 0) ---- //




		if (dis > dx*dx | iflag_total == 0) {

			Ndis_min = MAX;

			for (itri =1; itri <= cube_trinum[icube]; itri++) {

				Ntemp = (tri_table[itri]-1)*N_line;

				BI_detect_iflag0(myid, Ntemp, &BIp0, &BIp1, &BIp2, N_line, &Ndis, &Ndis_min,  pNode_inf, ptri_number, pNode, 
					tri, vert0, vert1, vert2,xcnt, ycnt, zcnt);

			}

			BIx = BIp0;
			BIy = BIp1;
			BIz = BIp2;

			IPx = 2*xcnt-BIx;
			IPy = 2*ycnt-BIy;
			IPz = 2*zcnt-BIz;

			dis = (BIx-xcnt)*(BIx-xcnt)+(BIy-ycnt)*(BIy-ycnt)+(BIz-zcnt)*(BIz-zcnt);

			if (dis > dx*dx) {
				
				FWS[icube][i][j][k] = IFLUID; 
				continue;
			
			}

		}
		


		
		//if(myid==0 && icube ==5 && i == 4 && j == 13 && k == 15) printf("2===%f\n",dis);
		

		/*
		for (itri = 0; itri < Ntri; itri++) {

			Ntemp = (itri)*N_line;

			BI_detect(myid, Ntemp, N_line, &Ndis, &Ndis_min, &iflag_total, &tmin, &dotp, &dir0, &dir1, &dir2, pNode_inf, ptri_number, pNode, 
				tri, vert0, vert1, vert2, nuvec, dir, orig, xcnt, ycnt, zcnt, &itemp, &ilarge, nlarge);


		}
		*/
		

		//// ---- if (dotp > 0) the point is ouside the object ---- //


		//if (iflag_total == 0) {

		//	Ndis_min = MAX;

		//	for (itri = 0; itri < Ntri; itri++) {

		//		Ntemp = (itri)*N_line;

		//		BI_detect_iflag0(myid, Ntemp, &BIp0, &BIp1, &BIp2, N_line, &Ndis, &Ndis_min,  pNode_inf, ptri_number, pNode, 
		//			tri, vert0, vert1, vert2,xcnt, ycnt, zcnt);

		//	}
		//}



		icount = 0;

		iicube = icube;
		ii = i;
		jj = j;
		kk = k;

		for (i = 0; i <= nxx; i++) {
			if (IPx >= Xcnt[icube][i] & IPx <= Xcnt[icube][i]+dx) 
			{

				icount = icount+1;
				break;

			}
		}   // ---- for (i = n_buffer; i <= nx; i++) ---- //

		for (j = 0; j <= nyy; j++) {
			if (IPy >= Ycnt[icube][j] & IPy <= Ycnt[icube][j]+dy) 
			{

				icount = icount+1;
				break;

			}
		}   // ---- for (j = 0; j < NcubeY; j++) { ---- //

		for (k = 0; k <= nzz; k++) {
			if ( IPz >= Zcnt[icube][k] & IPz <= Zcnt[icube][k]+dz) 
			{
				icount = icount+1;
				break;

			}
		}   // ---- for (k = 0; k <= NcubeZ; k++) ---- //


		if (icount < 3) {

			continue;

		}

		

		// ---- 000 => 100 => 010 => 001 => 110 => 101 => 011 => 111 ---- //

		s1 = IPx-Xcnt[icube][i];
		s2 = IPy-Ycnt[icube][j];
		s3 = IPz-Zcnt[icube][k];


		sur_x[1] = 0;
		sur_x[2] = Xcnt[icube][i+1]-Xcnt[icube][i];
		sur_x[3] = 0;
		sur_x[4] = 0;
		sur_x[5] = Xcnt[icube][i+1]-Xcnt[icube][i];
		sur_x[6] = Xcnt[icube][i+1]-Xcnt[icube][i];
		sur_x[7] = 0;
		sur_x[8] = Xcnt[icube][i+1]-Xcnt[icube][i];

		sur_y[1] = 0;
		sur_y[2] = 0;
		sur_y[3] = Ycnt[icube][j+1]-Ycnt[icube][j];
		sur_y[4] = 0;
		sur_y[5] = Ycnt[icube][j+1]-Ycnt[icube][j];
		sur_y[6] = 0;
		sur_y[7] = Ycnt[icube][j+1]-Ycnt[icube][j];
		sur_y[8] = Ycnt[icube][j+1]-Ycnt[icube][j];

		sur_z[1] = 0;
		sur_z[2] = 0;
		sur_z[3] = 0;
		sur_z[4] = Zcnt[icube][k+1]-Zcnt[icube][k];
		sur_z[5] = 0;
		sur_z[6] = Zcnt[icube][k+1]-Zcnt[icube][k];
		sur_z[7] = Zcnt[icube][k+1]-Zcnt[icube][k];
		sur_z[8] = Zcnt[icube][k+1]-Zcnt[icube][k];


		for (int vi = 0; vi < 8; vi++) {

			V[vi*8+7] = 1;

			V[vi*8+6] = sur_z[vi+1];

			V[vi*8+5] = sur_y[vi+1];

			V[vi*8+4] = sur_x[vi+1];

			V[vi*8+3] = sur_y[vi+1]*sur_z[vi+1];

			V[vi*8+2] = sur_x[vi+1]*sur_z[vi+1];

			V[vi*8+1] = sur_x[vi+1]*sur_y[vi+1];

			V[vi*8+0] = sur_x[vi+1]*sur_y[vi+1]*sur_z[vi+1];

		}


		GetInverseMatrix( V, V, 8 );

		for (int vi = 0; vi < 8; vi++) {

			weight[vi] = V[vi]*s1*s2*s3+V[vi+8]*s1*s2+V[vi+16]*s1*s3+V[vi+24]*s2*s3+V[vi+32]*s1+V[vi+40]*s2+V[vi+48]*s3+V[vi+56];

		}



		temp = weight[0]+weight[1]+weight[2]+weight[3]+weight[4]+weight[5]+weight[6]+weight[7];


		if ( 
			(weight[0] >=(-minimum) && weight[0] <= (1+minimum)) &&
			(weight[1] >=(-minimum) && weight[1] <= (1+minimum)) &&
			(weight[2] >=(-minimum) && weight[2] <= (1+minimum)) &&
			(weight[3] >=(-minimum) && weight[3] <= (1+minimum)) &&
			(weight[4] >=(-minimum) && weight[4] <= (1+minimum)) &&
			(weight[5] >=(-minimum) && weight[5] <= (1+minimum)) &&
			(weight[6] >=(-minimum) && weight[6] <= (1+minimum)) &&
			(weight[7] >=(-minimum) && weight[7] <= (1+minimum)) 
			)

		{ 
			s1 = orig[0]-BIx;
			s2 = orig[1]-BIy;
			s3 = orig[2]-BIz;

			tmin = sqrt(s1*s1+s2*s2+s3*s3)+0.0000001;

			if (dotp >= 0) {

				Ntemp = (igc-1)*4;
				fprintf(fptr_minus,"%d\t%d\t%d\t%d\n",GC[Ntemp+1],GC[Ntemp+2],GC[Ntemp+3],GC[Ntemp+4]);


				iNgc_minus = iNgc_minus+1;
				fprintf(fptr_minus,"%d\t%d\t%d\t%d\n",icube,i,j,k);
				fprintf(fptr_minus,"%f\t%f\t%f\n",s1/tmin,s2/tmin,s3/tmin);
				fprintf(fptr_minus,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",sqrt(dis),weight[0],weight[1],weight[2],weight[3],weight[4],weight[5],weight[6],weight[7]);

			}
			else {

				Ntemp = (igc-1)*4;
				fprintf(fptr_plus,"%d\t%d\t%d\t%d\n",GC[Ntemp+1],GC[Ntemp+2],GC[Ntemp+3],GC[Ntemp+4]);

				iNgc_plus = iNgc_plus+1;
				fprintf(fptr_plus,"%d\t%d\t%d\t%d\n",icube,i,j,k);
				fprintf(fptr_plus,"%f\t%f\t%f\n",s1/tmin,s2/tmin,s3/tmin);
				fprintf(fptr_plus,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",sqrt(dis),weight[0],weight[1],weight[2],weight[3],weight[4],weight[5],weight[6],weight[7]);

			}


		}


		/*
		else {

		FWS[iicube][ii][jj][kk] = -50;
		printf("%d\t%d\t%d\t%d\n",GC[Ntemp+1],GC[Ntemp+2],GC[Ntemp+3],GC[Ntemp+4]);
		printf("iflag_total = %d\n dotp = %.16f\n ",iflag_total,dotp);
		printf("BI = %f\t%f\t%f\n",BIx,BIy,BIz);
		printf("IP = %f\t%f\t%f\n",IPx,IPy,IPz);
		printf("orig = %f\t%f\t%f\n",orig[0],orig[1],orig[2]);
		printf("vert = %f\t%f\t%f\n",vert0_iflag_total0_x,vert0_iflag_total0_y,vert0_iflag_total0_z);

		}
		*/
	}    // ---- for (igc = 1; igc <= Ngc; igc++) ---- //

	fprintf(fptr_minus,"Number of Boundary Cells\n");
	fprintf(fptr_minus,"%d",iNgc_minus);

	fprintf(fptr_plus,"Number of Boundary Cells\n");
	fprintf(fptr_plus,"%d",iNgc_plus);


	// ==== immersed boundary file ===== //
	fclose(fptr_plus);           //
	fclose(fptr_minus);			 //
	// ==== immersed boundary file ===== //

	delete[] tri_table;
	delete[] cube_trinum;
	delete[] tri;
	delete[] pvertex;
	delete[] pNode;
	delete[] ptri_number;
	delete[] GC;
	delete[] GCcnt;

}



int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3])

{

  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) this gives 3x3=9 more tests */

   double v0[3],v1[3],v2[3];

//   float axis[3];

   double min,max,p0,p1,p2,rad,fex,fey,fez;		// -NJMP- "d" local variable removed

   double normal[3],e0[3],e1[3],e2[3];



   /* This is the fastest branch on Sun */
   /* move everything so that the boxcenter is in (0,0,0) */

   SUB(v0,triverts[0],boxcenter);
   SUB(v1,triverts[1],boxcenter);
   SUB(v2,triverts[2],boxcenter);

   /* compute triangle edges */

   SUB(e0,v1,v0);      /* tri edge 0 */
   SUB(e1,v2,v1);      /* tri edge 1 */
   SUB(e2,v0,v2);      /* tri edge 2 */


   /* Bullet 3:  */

   /*  test the 9 tests first (this was faster) */

   fex = fabsf(e0[X]);
   fey = fabsf(e0[Y]);
   fez = fabsf(e0[Z]);

   AXISTEST_X01(e0[Z], e0[Y], fez, fey);
   AXISTEST_Y02(e0[Z], e0[X], fez, fex);
   AXISTEST_Z12(e0[Y], e0[X], fey, fex);



   fex = fabsf(e1[X]);
   fey = fabsf(e1[Y]);
   fez = fabsf(e1[Z]);

   AXISTEST_X01(e1[Z], e1[Y], fez, fey);
   AXISTEST_Y02(e1[Z], e1[X], fez, fex);
   AXISTEST_Z0(e1[Y], e1[X], fey, fex);



   fex = fabsf(e2[X]);
   fey = fabsf(e2[Y]);
   fez = fabsf(e2[Z]);

   AXISTEST_X2(e2[Z], e2[Y], fez, fey);
   AXISTEST_Y1(e2[Z], e2[X], fez, fex);
   AXISTEST_Z12(e2[Y], e2[X], fey, fex);



   /* Bullet 1: */
   /*  first test overlap in the {x,y,z}-directions */
   /*  find min, max of the triangle each direction, and test for overlap in */
   /*  that direction -- this is equivalent to testing a minimal AABB around */
   /*  the triangle against the AABB */



   /* test in X-direction */

   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;


   /* test in Y-direction */

   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;


   /* test in Z-direction */

   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;



   /* Bullet 2: */
   /*  test if the box intersects the plane of the triangle */
   /*  compute plane equation of triangle: normal*x+d=0 */

   CROSS(normal,e0,e1);

   // -NJMP- (line removed here)

   if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;	// -NJMP-

   return 1;   /* box and triangle overlaps */

}



int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3])	// -NJMP-

{

  int q;

  float vmin[3],vmax[3],v;

  for(q=X;q<=Z;q++)

  {

    v=vert[q];					// -NJMP-

    if(normal[q]>0.0f)

    {

      vmin[q]=-maxbox[q] - v;	// -NJMP-

      vmax[q]= maxbox[q] - v;	// -NJMP-

    }

    else

    {

      vmin[q]= maxbox[q] - v;	// -NJMP-

      vmax[q]=-maxbox[q] - v;	// -NJMP-

    }

  }

  if(DOT(normal,vmin)>0.0f) return 0;	// -NJMP-

  if(DOT(normal,vmax)>=0.0f) return 1;	// -NJMP-


  return 0;

}








void BI_detect
	(
	// ======================= //
	int myid,

	int Ntemp,

	int N_line,

	double *Ndis,

	double *Ndis_min,

	int *iflag_total,

	double *tmin,

	double *dotp,

	double *dir0,
	double *dir1,
	double *dir2,

	int pNode_inf[],
	int ptri_number[],
	int pNode[],

	double tri[],

	double vert0[3],
	double vert1[3],
	double vert2[3],

	double nuvec[3],

	double dir[3],

	double orig[3],

	double xcnt,
	double ycnt,
	double zcnt,

	int *itemp,
	int *ilarge,
	int nlarge[]

// ======================= //
)
{


	int ishare;
	int iflag;


	double v12[3], v13[3];
	double dis;


	double t, u, v;



	vert0[0] = tri[Ntemp+5];
	vert0[1] = tri[Ntemp+6];
	vert0[2] = tri[Ntemp+7];

	vert1[0] = tri[Ntemp+9];
	vert1[1] = tri[Ntemp+10];
	vert1[2] = tri[Ntemp+11];

	vert2[0] = tri[Ntemp+13];
	vert2[1] = tri[Ntemp+14];
	vert2[2] = tri[Ntemp+15];

	v12[0] = vert1[0] - vert0[0];
	v12[1] = vert1[1] - vert0[1];
	v12[2] = vert1[2] - vert0[2];

	v13[0] = vert2[0] - vert0[0];
	v13[1] = vert2[1] - vert0[1];
	v13[2] = vert2[2] - vert0[2];


	nuvec[0] = v12[1] * v13[2] - v12[2] * v13[1];
	nuvec[1] = v12[2] * v13[0] - v12[0] * v13[2];
	nuvec[2] = v12[0] * v13[1] - v12[1] * v13[0];

	dis = sqrt(nuvec[0]*nuvec[0]+nuvec[1]*nuvec[1]+nuvec[2]*nuvec[2]);

	dir[0] = -nuvec[0]/dis;
	dir[1] = -nuvec[1]/dis;
	dir[2] = -nuvec[2]/dis;

	orig[0] = xcnt;
	orig[1] = ycnt;
	orig[2] = zcnt;

	t = 0;
	iflag = intersect_triangle(orig, dir, vert0, vert1, vert2, &t, &u, &v);
	*Ndis = fabs(t)*(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);

	

	if (iflag == 1) {

		*iflag_total = *iflag_total + 1;

		if (*Ndis < *Ndis_min)
		{
			*Ndis_min = *Ndis;
			*tmin = fabs(t);
			*dir0 = dir[0];
			*dir1 = dir[1];
			*dir2 = dir[2];
			*dotp = nuvec[0]*(vert0[0]-xcnt)+nuvec[1]*(vert0[1]-ycnt)+nuvec[2]*(vert0[2]-zcnt);

			/*if(  myid == 0 & *itemp == 100) {

				printf("%d\t%f\t%f\t%f\t%f\n",Ntemp,tri[Ntemp+5],tri[Ntemp+6],tri[Ntemp+7],*Ndis);
				printf("%d\t%f\t%f\t%f\t%f\n",Ntemp,tri[Ntemp+9],tri[Ntemp+10],tri[Ntemp+11],*Ndis);
				printf("%d\t%f\t%f\t%f\t%f\n\n",Ntemp,tri[Ntemp+13],tri[Ntemp+14],tri[Ntemp+15],*Ndis);

			}*/


		}
	}    // ---- if (iflag == 1) ---- //


}


void BI_detect_iflag0
	(

	// ======================= //
	int myid,

	int Ntemp,

	double *BIp0,
	double *BIp1,
	double *BIp2,

	int N_line,

	double *Ndis,

	double *Ndis_min,

	int pNode_inf[],
	int ptri_number[],
	int pNode[],

	double tri[],

	double vert0[3],
	double vert1[3],
	double vert2[3],

	double xcnt,
	double ycnt,
	double zcnt

	// ======================= //
	)

{

	int ishare;
	int iflag;


	double v12[3], v13[3];
	double dis;
	double t, u, v;
	double p0,p1,p2;

	vert0[0] = tri[Ntemp+5];
	vert0[1] = tri[Ntemp+6];
	vert0[2] = tri[Ntemp+7];

	vert1[0] = tri[Ntemp+9];
	vert1[1] = tri[Ntemp+10];
	vert1[2] = tri[Ntemp+11];

	vert2[0] = tri[Ntemp+13];
	vert2[1] = tri[Ntemp+14];
	vert2[2] = tri[Ntemp+15];

	ClosestPTPoinTrangle(myid, &p0,  &p1,  &p2, vert0, vert1, vert2, xcnt, ycnt, zcnt);

	*Ndis = (p0-xcnt)*(p0-xcnt)+(p1-ycnt)*(p1-ycnt)+(p2-zcnt)*(p2-zcnt);

	if (*Ndis <= *Ndis_min) {

		*BIp0 = p0;
		*BIp1 = p1;
		*BIp2 = p2;

		*Ndis_min = *Ndis;

	}


}


void ClosestPTPoinTrangle
	(

	// ======================= //
	int myid,

	double *p0,
	double *p1,
	double *p2,

	double vert0[3],
	double vert1[3],
	double vert2[3],

	double xcnt,
	double ycnt,
	double zcnt

	// ======================= //
	)
{

	int ishare;
	int Ntemp;
	int iflag;


	double ab[3], ac[3], bc[3], ap[3], bp[3], cp[3];
	double d1,d2,d3,d4,d5,d6,vc,v,vb,w,va,denom;

	ab[0] = vert1[0] - vert0[0];
	ab[1] = vert1[1] - vert0[1];
	ab[2] = vert1[2] - vert0[2];

	ac[0] = vert2[0] - vert0[0];
	ac[1] = vert2[1] - vert0[1];
	ac[2] = vert2[2] - vert0[2];

	bc[0] = vert2[0] - vert1[0];
	bc[1] = vert2[1] - vert1[1];
	bc[2] = vert2[2] - vert1[2];

	ap[0] = xcnt - vert0[0];
	ap[1] = ycnt - vert0[1];
	ap[2] = zcnt - vert0[2];


	d1 = ab[0]*ap[0]+ab[1]*ap[1]+ab[2]*ap[2];
	d2 = ac[0]*ap[0]+ac[1]*ap[1]+ac[2]*ap[2];

	if( d1 <= 0. && d2 <= 0. ) {

		*p0 = vert0[0];
		*p1 = vert0[1];
		*p2 = vert0[2];

		return;

	}

	bp[0] = xcnt - vert1[0];
	bp[1] = ycnt - vert1[1];
	bp[2] = zcnt - vert1[2];

	d3 = ab[0]*bp[0]+ab[1]*bp[1]+ab[2]*bp[2];
	d4 = ac[0]*bp[0]+ac[1]*bp[1]+ac[2]*bp[2];

	if( d3 >= 0. && d4 <= d3 ) {

		*p0 = vert1[0];
		*p1 = vert1[1];
		*p2 = vert1[2];

		return;

	}

	vc = d1*d4-d3*d2;

	if(vc <= 0. && d1 >= 0.&& d3 <= 0.) {

		v = d1/(d1-d3);

		*p0 = vert0[0]+v*ab[0];
		*p1 = vert0[1]+v*ab[1];
		*p2 = vert0[2]+v*ab[2];

		return;

	}


	cp[0] = xcnt - vert2[0];
	cp[1] = ycnt - vert2[1];
	cp[2] = zcnt - vert2[2];

	d5 = ab[0]*cp[0]+ab[1]*cp[1]+ab[2]*cp[2];
	d6 = ac[0]*cp[0]+ac[1]*cp[1]+ac[2]*cp[2];

	if(d6 >= 0. && d5 <= d6){

		*p0 = vert2[0];
		*p1 = vert2[1];
		*p2 = vert2[2];

		return;

	}


	vb = d5*d2-d1*d6;
	if(vb <= 0. && d2 >= 0. && d6 <= 0.) {

		w = d2/(d2-d6);


		*p0 = vert0[0]+w*ac[0];
		*p1 = vert0[1]+w*ac[1];
		*p2 = vert0[2]+w*ac[2];

		return;

	}

	va = d3*d6-d5*d4;
	if(va <= 0. && (d4-d3) >= 0. && (d5-d6) >= 0.) {

		w = (d4-d3)/((d4-d3)+(d5-d6));

		*p0 = vert1[0]+w*bc[0];
		*p1 = vert1[1]+w*bc[1];
		*p2 = vert1[2]+w*bc[2];

		return;

	}

	denom = 1./(va+vb+vc);
	v = vb*denom;
	w = vc+denom;

	*p0 = vert0[0]+ab[0]*v+ac[0]*w;
	*p1 = vert0[1]+ab[1]*v+ac[1]*w;
	*p2 = vert0[2]+ab[2]*v+ac[2]*w;


}




int intersect_triangle
	(
	// ======================= //
	double orig[3],
	double dir[3],
	double vert0[3],
	double vert1[3],
	double vert2[3],
	double *t,
	double *u,
	double *v
	// ======================= //
	)
{

	double edge1[3], edge2[3], tvec[3], pvec[3], qvec[3];
	double det, inv_det;



	SUB(edge1, vert1, vert0);
	SUB(edge2, vert2, vert0);

	CROSS(pvec, dir, edge2);

	det = DOT(edge1, pvec);



	if (det < 0.00000000000001)
		return 0;

	SUB(tvec, orig, vert0);

	*u = DOT(tvec, pvec);

	if (*u < 0.0 || *u > det)
		return 0;

	CROSS(qvec, tvec, edge1);

	*v = DOT(dir, qvec);

	if (*v < 0.0 || *u+*v > det)
		return 0;

	*t = DOT(edge2, qvec);
	inv_det = 1.0 / det;

	*t = DOT(edge2, qvec);

	inv_det = 1.0 / det;
	*t *= inv_det;
	*u *= inv_det;
	*v *= inv_det;

	return 1;

}


void GetInverseMatrix( const double* inMat, double* outMat, const int n )
{
	double tmp;
	int i, j, k;
	double **a;
	double **inv_a;


	a = (double**)malloc(sizeof(double*)*n);
	a[0] = (double*)malloc(sizeof(double)*(n*n));
	for( i=1; i<n; i++ ){
		a[i] = a[i-1] + n;
	}
	inv_a = (double**)malloc(sizeof(double*)*n);
	inv_a[0] = (double*)malloc(sizeof(double)*(n*n));
	for( i=1;i<n;i++){
		inv_a[i] = inv_a[i-1] + n;
	}

	for( j=0; j<n; j++ ){
		for( i=0; i<n; i++ ){
			a[j][i] = inMat[i + n*j];
		}
	}

	for( j=0; j<n; j++ ){
		for( i=0; i<n; i++ ){
			inv_a[j][i] = ( i==j ) ? 1.0 : 0.0;
		}
	}

	for( k=0; k<n; k++){

		int max = k;
		for( j=k+1; j<n; j++){
			if( fabs(a[j][k]) > fabs(a[max][k]) ){
				max = j;
			}
		}

		if( max != k ){
			for( i=0; i<n; i++ ){

				tmp = a[max][i];
				a[max][i] = a[k][i];
				a[k][i] = tmp;

				tmp = inv_a[max][i];
				inv_a[max][i] = inv_a[k][i];
				inv_a[k][i] = tmp;
			}
		}

		tmp = a[k][k];
		for(i=0;i<n;i++){
			a[k][i] /= tmp;
			inv_a[k][i] /= tmp;
		}

		for( j=0;j<n;j++ ){
			if( j != k ){
				tmp =   a[j][k] / a[k][k];
				for(i=0;i<n;i++){
					a[j][i] = a[j][i] - a[k][i] * tmp;
					inv_a[j][i] = inv_a[j][i] - inv_a[k][i] * tmp;
				}
			}
		}

	}
	for( j=0; j<n; j++ ){
		for( i=0; i<n; i++ ){
			outMat[ i+n*j ] = inv_a[j][i];

		}
	}

	free(a[0]);
	free(a);
	free(inv_a[0]);
	free(inv_a);
}

