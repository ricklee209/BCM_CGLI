



#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <string>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "main.h"
#include "prm.h"

int main(int argc, char **argv)
{   

	//==== MPI start ====//

	int myid;
	int nproc = np;

	// ============================================ //
	MPI_Status istat[8];					    //
	//
	MPI_Comm comm;								//
	//
	MPI_Init (&argc, &argv);					//
	//
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);		//
	//
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);		//
	//
	comm=MPI_COMM_WORLD;						//
	// ============================================ //


#include "Resolution.h"

	int statistic_step = 1089;    // ---- periodic step ---- //

	int start_step = 6000;    // ---- how many steps for reaching the quasi steady ---- //

	int dp_step = 500;    // ---- how many steps for periodically outputing the dp ---- //

	int iteration_end_step = 1;
	int output_step = 1;
	int count = 1;	
	int step;



	double deltaT = 0.0001;
	double deltaTau = deltaT/200.;
	double e = 0.1;


	int NBC,NBC_plus, NBC_minus, Ntemp, gicube, gi, gj, gk, mp_switch;

	double (*er) = new double[8];

	double E1 = 0;
	double E2 = 0;
	double E3 = 0;
	double E4 = 0;
	double E5 = 0;

#include "BCM.h"
#include "Resolution.h"

	FILE *fptr;
	char file_name[100];
	char str[1024];
	char str_ncube[1024] = "#_of_cubes in the computational region >>";
	char str_n_wallcube[1024] = "#_of_wall_cubes >>";

	int Ncube;    
	int N_wallcube;    

	sprintf(file_name,"BCMgrid""%0.5d"".dat",myid);    
	fptr = fopen(file_name,"r"); 

	// ---- Cube numbers ---- //
	while( strcmp(str,str_ncube) != 0 )  { fscanf(fptr,"%[^\n]\n",str); }
	fscanf(fptr,"%d\n",&Ncube); 
	Ncube = Ncube + 1;    

	// ---- Wall cube numbers ---- // 
	while( strcmp(str,str_n_wallcube) != 0 )  { fscanf(fptr,"%[^\n]\n",str); }
	fscanf(fptr,"%d\n",&N_wallcube);
	N_wallcube = N_wallcube + 1;    

	fclose (fptr);

#include "Array.h"    // ---- allocate the memory ----//


	BCM_Grid(myid, Ncube, N_wallcube, cube_size, csl, Xcube, Ycube, Zcube, Xcnt, Ycnt, Zcnt, FWS, adj_number, wallcube);


	// -------------------------------------------------------------------------------------------------------------------- //
	// ------------------------------------------------- Detect boundary  ------------------------------------------------- //

	double Xmax = MIN;    /**** initialize ****/
	double Xmin = MAX;    /**** initialize ****/

	double Ymax = MIN;    /**** initialize ****/
	double Ymin = MAX;    /**** initialize ****/

	double Zmax = MIN;    /**** initialize ****/
	double Zmin = MAX;    /**** initialize ****/


	for (icube = 1; icube < Ncube; icube++) {  

		Xmin = min(Xmin,Xcube[icube]);
		Xmax = max(Xmax,Xcube[icube]+cube_size[icube]);

		Ymin = min(Ymin,Ycube[icube]);
		Ymax = max(Ymax,Ycube[icube]+cube_size[icube]);

		Zmin = min(Zmin,Zcube[icube]);
		Zmax = max(Zmax,Zcube[icube]+cube_size[icube]);

	}

	double gXmax = MIN;    /**** initialize ****/
	double gXmin = MAX;    /**** initialize ****/

	double gYmax = MIN;    /**** initialize ****/
	double gYmin = MAX;    /**** initialize ****/

	double gZmax = MIN;    /**** initialize ****/
	double gZmin = MAX;    /**** initialize ****/

	MPI_Allreduce ((void*)&Xmax, (void*)&gXmax, 1, MPI_DOUBLE, MPI_MAX, comm );
	MPI_Allreduce ((void*)&Xmin, (void*)&gXmin, 1, MPI_DOUBLE, MPI_MIN, comm );
	MPI_Allreduce ((void*)&Ymax, (void*)&gYmax, 1, MPI_DOUBLE, MPI_MAX, comm );
	MPI_Allreduce ((void*)&Ymin, (void*)&gYmin, 1, MPI_DOUBLE, MPI_MIN, comm );
	MPI_Allreduce ((void*)&Zmax, (void*)&gZmax, 1, MPI_DOUBLE, MPI_MAX, comm );
	MPI_Allreduce ((void*)&Zmin, (void*)&gZmin, 1, MPI_DOUBLE, MPI_MIN, comm );

	int NXbc_l = 0;
	int NXbc_u = 0;

	int NYbc_l = 0;
	int NYbc_u = 0;

	int NZbc_l = 0;
	int NZbc_u = 0;

	for (icube = 1; icube < Ncube; icube++) {  

		if (fabs(Xcube[icube]+cube_size[icube]-gXmax) < minimum) NXbc_u = NXbc_u+1;
		if (fabs(Ycube[icube]+cube_size[icube]-gYmax) < minimum) NYbc_u = NYbc_u+1;
		if (fabs(Zcube[icube]+cube_size[icube]-gZmax) < minimum) NZbc_u = NZbc_u+1;

		if (fabs(Xcube[icube]-gXmin) < minimum) NXbc_l = NXbc_l+1;
		if (fabs(Ycube[icube]-gYmin) < minimum) NYbc_l = NYbc_l+1;
		if (fabs(Zcube[icube]-gZmin) < minimum) NZbc_l = NZbc_l+1;

	};

	int (*Xbc_u) = new int[NXbc_u+1]; 
	int (*Xbc_l) = new int[NXbc_l+1]; 

	int (*Ybc_u) = new int[NYbc_u+1]; 
	int (*Ybc_l) = new int[NYbc_l+1]; 

	int (*Zbc_u) = new int[NZbc_u+1]; 
	int (*Zbc_l) = new int[NZbc_l+1]; 

	NXbc_l = NXbc_u = NYbc_l = NYbc_u = NZbc_l = NZbc_u = 0;

	for (icube = 1; icube < Ncube; icube++) {  

		if (fabs(Xcube[icube]+cube_size[icube]-gXmax) < minimum) { NXbc_u = NXbc_u+1; Xbc_u[NXbc_u] = icube; }
		if (fabs(Ycube[icube]+cube_size[icube]-gYmax) < minimum) { NYbc_u = NYbc_u+1; Ybc_u[NYbc_u] = icube; }
		if (fabs(Zcube[icube]+cube_size[icube]-gZmax) < minimum) { NZbc_u = NZbc_u+1; Zbc_u[NZbc_u] = icube; }

		if (fabs(Xcube[icube]-gXmin) < minimum) { NXbc_l = NXbc_l+1; Xbc_l[NXbc_l] = icube; }
		if (fabs(Ycube[icube]-gYmin) < minimum) { NYbc_l = NYbc_l+1; Ybc_l[NYbc_l] = icube; }
		if (fabs(Zcube[icube]-gZmin) < minimum) { NZbc_l = NZbc_l+1; Zbc_l[NZbc_l] = icube; }

	};

	// ------------------------------------------------- Detect boundary  ------------------------------------------------- //	
	// -------------------------------------------------------------------------------------------------------------------- //





	// ---------------------------- for each CPU ---------------------------- //

	int nadjX_eq, nadjY_eq, nadjZ_eq;

	int nadjX_bs_plus, nadjX_sb_plus, nadjX_bs_minus, nadjX_sb_minus;    // ---- bs means big to small ---- //
	int nadjY_bs_plus, nadjY_sb_plus, nadjY_bs_minus, nadjY_sb_minus;    // ---- sb means small to big ---- //
	int nadjZ_bs_plus, nadjZ_sb_plus, nadjZ_bs_minus, nadjZ_sb_minus;    // ---- eq means equal ---- //

	// ---------------------------- for each CPU ---------------------------- //


	// ------------------ reading communication_table for MPI ------------------ //

	int send_size, recv_size;

	int MPI_Nadj;

	fptr = fopen("MPI_communication_table_in_turn.dat","r"); 

	char str_MPIcon[1024] = "Number of Cubes for connection  >>";

	while( strcmp(str,str_MPIcon) != 0 )  { fscanf(fptr,"%[^\n]\n",str); }
	fscanf(fptr,"%d\n",&MPI_Nadj);
	MPI_Nadj = MPI_Nadj + 1;

	fclose(fptr);


	int (*MPI_cube) = new int[MPI_Nadj];
	int (*MPI_cpu) = new int[MPI_Nadj];

	int (*MPI_cube_adj) = new int[MPI_Nadj];
	int (*MPI_cpu_adj) = new int[MPI_Nadj];

	int (*MPI_direction) = new int[MPI_Nadj];
	int (*MPI_interface) = new int[MPI_Nadj];


	// ------------------ reading communication_table for MPI ------------------ //




	// ---------------------------- for MPI ---------------------------- //

	int MadjX_bs, MadjX_eq, MadjX_sb;    // ---- bs means big to small ---- //
	int MadjY_bs, MadjY_eq, MadjY_sb;    // ---- sb means small to big ---- //
	int MadjZ_bs, MadjZ_eq, MadjZ_sb;    // ---- eq means equal ---- //

	int Ncpu_bs, Ncpu_eq, Ncpu_sb;    // ---- number of neighbor CPU ---- //

	int Scube_bs, Scube_eq, Scube_sb;    // ---- how many cube in each CPU for Sending ---- //

	int Rcube_bs, Rcube_eq, Rcube_sb;    // ---- how many cube in each CPU for Receiving ---- //

	int Max_nei_bs, Max_nei_eq, Max_nei_sb;

	// ---------------------------- for MPI ---------------------------- //




	// ------------------ reading communication_table each CPU ------------------ //

	sprintf(file_name,"Communication_irank.""%0.5d"".dat",myid);   
	fptr = fopen(file_name,"r"); 

	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);



	// ======================================================================================================= //
	// ====================================== # The same size adjancent ====================================== //


	fscanf(fptr,"%d\n",&Ncpu_eq);


	int (*neighbor_cpu_eq) = new int[Ncpu_eq];
	int (*Ncube_Ncpu_eq) = new int[Ncpu_eq];    // ---- how many cubes in each CPU --- //

	fscanf(fptr,"%[^\n]\n",str);

	for (int neicpu = 0; neicpu < Ncpu_eq; neicpu++) {

		fscanf(fptr,"%d\n",&neighbor_cpu_eq[neicpu]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	for (int neicpu = 0; neicpu < Ncpu_eq; neicpu++) {

		fscanf(fptr,"%d\n",&Ncube_Ncpu_eq[neicpu]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	Max_nei_eq = 0;
	for (int neicpu = 0; neicpu < Ncpu_eq; neicpu++) {

		Max_nei_eq = Max_nei_eq+Ncube_Ncpu_eq[neicpu];

	}

	int (*Scube_Ncpu_eq) = new int[Max_nei_eq+1];    // ---- Cube number in each CPU for Sending ---- //
	int (*Rcube_Ncpu_eq) = new int[Max_nei_eq+1];    // ---- Cube number in each CPU for Receiving ---- //


	double (*send_data_curr_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1];
	double (*recv_data_curr_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1];

	double (*send_data_neig_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1];
	double (*recv_data_neig_eq) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_eq+1];


	for (int icube = 1; icube < Max_nei_eq+1; icube++) {

		fscanf(fptr,"%d\n",&Scube_Ncpu_eq[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	for (int icube = 1; icube < Max_nei_eq+1; icube++) {

		fscanf(fptr,"%d\n",&Rcube_Ncpu_eq[icube]);

	}


	int (*Sdir_eq) = new int[Max_nei_eq+1];    // ---- direction index in each CPU for Sending ---- //
	int (*Rdir_eq) = new int[Max_nei_eq+1];    // ---- direction index in each CPU for Sending ---- //

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_eq+1; icube++) {

		fscanf(fptr,"%d\n",&Sdir_eq[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_eq+1; icube++) {

		fscanf(fptr,"%d\n",&Rdir_eq[icube]);

	}

	// ====================================== # The same size adjancent ====================================== //
	// ======================================================================================================= //


	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);


	// ============================================================================================ //
	// ====================================== # Small to big ====================================== //


	fscanf(fptr,"%d\n",&Ncpu_sb);


	int (*neighbor_cpu_sb) = new int[Ncpu_sb];
	int (*Ncube_Ncpu_sb) = new int[Ncpu_sb];    // ---- how many cubes in each CPU --- //

	fscanf(fptr,"%[^\n]\n",str);

	for (int neicpu = 0; neicpu < Ncpu_sb; neicpu++) {

		fscanf(fptr,"%d\n",&neighbor_cpu_sb[neicpu]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	for (int neicpu = 0; neicpu < Ncpu_sb; neicpu++) {

		fscanf(fptr,"%d\n",&Ncube_Ncpu_sb[neicpu]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	Max_nei_sb = 0;
	for (int neicpu = 0; neicpu < Ncpu_sb; neicpu++) {

		Max_nei_sb = Max_nei_sb+Ncube_Ncpu_sb[neicpu];

	}

	int (*Scube_Ncpu_sb) = new int[Max_nei_sb+1];    // ---- Cube number in each CPU for Sending ---- //
	int (*Rcube_Ncpu_sb) = new int[Max_nei_sb+1];    // ---- Cube number in each CPU for Receiving ---- //


	double (*send_data_curr_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1];
	double (*recv_data_curr_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1];

	double (*send_data_neig_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1];
	double (*recv_data_neig_sb) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_sb+1];



	for (int icube = 1; icube < Max_nei_sb+1; icube++) {

		fscanf(fptr,"%d\n",&Scube_Ncpu_sb[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	for (int icube = 1; icube < Max_nei_sb+1; icube++) {

		fscanf(fptr,"%d\n",&Rcube_Ncpu_sb[icube]);

	}

	int (*Sdir_sb) = new int[Max_nei_sb+1];    // ---- direction index in each CPU for Sending ---- //
	int (*Rdir_sb) = new int[Max_nei_sb+1];    // ---- direction index in each CPU for Sending ---- //

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_sb+1; icube++) {

		fscanf(fptr,"%d\n",&Sdir_sb[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_sb+1; icube++) {

		fscanf(fptr,"%d\n",&Rdir_sb[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);
	int (*adjN_sb) = new int[Max_nei_sb+1];    // ---- Adjacent order of current cube(small) for adjacent cube(big) ---- //


	for (int icube = 1; icube < Max_nei_sb+1; icube++) {

		fscanf(fptr,"%d\n",&adjN_sb[icube]);

	}


	// ====================================== # Small to big ====================================== //
	// ============================================================================================ //


	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);


	// ====================================== # Big to small ====================================== //
	// ============================================================================================ //

	fscanf(fptr,"%d\n",&Ncpu_bs);


	int (*neighbor_cpu_bs) = new int[Ncpu_bs];
	int (*Ncube_Ncpu_bs) = new int[Ncpu_bs];    // ---- how many cubes in each CPU --- //

	fscanf(fptr,"%[^\n]\n",str);

	for (int neicpu = 0; neicpu < Ncpu_bs; neicpu++) {

		fscanf(fptr,"%d\n",&neighbor_cpu_bs[neicpu]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	for (int neicpu = 0; neicpu < Ncpu_bs; neicpu++) {

		fscanf(fptr,"%d\n",&Ncube_Ncpu_bs[neicpu]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	Max_nei_bs = 0;
	for (int neicpu = 0; neicpu < Ncpu_bs; neicpu++) {

		Max_nei_bs = Max_nei_bs+Ncube_Ncpu_bs[neicpu];

	}

	int (*Scube_Ncpu_bs) = new int[Max_nei_bs+1];    // ---- Cube number in each CPU for Sending ---- //
	int (*Rcube_Ncpu_bs) = new int[Max_nei_bs+1];    // ---- Cube number in each CPU for Receiving ---- //

	double (*send_data_curr_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1];
	double (*recv_data_curr_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1];

	double (*send_data_neig_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1];
	double (*recv_data_neig_bs) = new double[5*NcubeX*NcubeY*n_buffer*Max_nei_bs+1];


	for (int icube = 1; icube < Max_nei_bs+1; icube++) {

		fscanf(fptr,"%d\n",&Scube_Ncpu_bs[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);

	for (int icube = 1; icube < Max_nei_bs+1; icube++) {

		fscanf(fptr,"%d\n",&Rcube_Ncpu_bs[icube]);

	}

	int (*Sdir_bs) = new int[Max_nei_bs+1];    // ---- direction index in each CPU for Sending ---- //
	int (*Rdir_bs) = new int[Max_nei_bs+1];    // ---- direction index in each CPU for Recving ---- //

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_bs+1; icube++) {

		fscanf(fptr,"%d\n",&Sdir_bs[icube]);

	}

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_bs+1; icube++) {

		fscanf(fptr,"%d\n",&Rdir_bs[icube]);

	}


	int (*RadjN_bs) = new int[Max_nei_bs+1];    // ---- Adjacent order of adjacent cube(small) for current cube(big) ---- //

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_bs+1; icube++) {

		fscanf(fptr,"%d\n",&RadjN_bs[icube]);

	}


	int (*SadjN_bs) = new int[Max_nei_bs+1];    // ---- Adjacent order of adjacent cube(small) for current cube(big) ---- //

	fscanf(fptr,"%[^\n]\n",str);
	for (int icube = 1; icube < Max_nei_bs+1; icube++) {

		fscanf(fptr,"%d\n",&SadjN_bs[icube]);

	}

	// ====================================== # Big to small ====================================== //
	// ============================================================================================ //

	fclose(fptr);

	// ------------------ reading communication_table each CPU ------------------ //




	// ---------------------- for each CPU ---------------------- //	
	BCM_Interface_table(myid, Ncube, N_wallcube,
		&nadjX_eq, &nadjY_eq, &nadjZ_eq,
		&nadjX_bs_plus, &nadjX_sb_plus, &nadjX_bs_minus, &nadjX_sb_minus,
		&nadjY_bs_plus, &nadjY_sb_plus, &nadjY_bs_minus, &nadjY_sb_minus,
		&nadjZ_bs_plus, &nadjZ_sb_plus, &nadjZ_bs_minus, &nadjZ_sb_minus,
		csl, 
		adj_number,
		adjX_eq, adjY_eq, adjZ_eq,
		adjX_bs_plus, adjX_sb_plus, adjX_bs_minus, adjX_sb_minus,
		adjY_bs_plus, adjY_sb_plus, adjY_bs_minus, adjY_sb_minus,
		adjZ_bs_plus, adjZ_sb_plus, adjZ_bs_minus, adjZ_sb_minus);
	// ---------------------- for each CPU ---------------------- //							




	int (*rank_map)[MPI_Ncube] = new int[2][MPI_Ncube];

	// ---- rank_map [Cubes(Global)] [Ranks] [Cubes(local)] ---- //

	// =========================================================== //

	sprintf(file_name,"rank_map.dat");             

	fptr = fopen(file_name,"r");    

	fscanf(fptr,"%[^\n]\n",str);

	for (icube = 0; icube < MPI_Ncube; icube++) {

		fscanf(fptr,"%d\t%d\t%d\n",&i,&rank_map[0][icube],&rank_map[1][icube]);

	}

	fclose(fptr);    // ---- close metis file ---- //

	// =========================================================== //



	// ---- calculate istart for MPI ---- //
	// =========================================================== //
	int (*ist_eq) = new int[Ncpu_eq];
	int (*ist_sb) = new int[Ncpu_sb];
	int (*ist_bs) = new int[Ncpu_bs];

	ist_eq[0] = 0;
	ist_sb[0] = 0;
	ist_bs[0] = 0;

	for (int neicpu = 1; neicpu < Ncpu_eq; neicpu++) {

		ist_eq[neicpu] = ist_eq[neicpu-1]+Ncube_Ncpu_eq[neicpu-1];

	}

	for (int neicpu = 1; neicpu < Ncpu_sb; neicpu++) {

		ist_sb[neicpu] = ist_sb[neicpu-1]+Ncube_Ncpu_sb[neicpu-1];

	}

	for (int neicpu = 1; neicpu < Ncpu_bs; neicpu++) {

		ist_bs[neicpu] = ist_bs[neicpu-1]+Ncube_Ncpu_bs[neicpu-1];

	}
	// =========================================================== //




	BCM_FWS_Interface(myid,Ncube, 
		MPI_Nadj,
		Ncpu_eq, 
		Max_nei_eq,
		nadjX_eq, nadjY_eq, nadjZ_eq,
		rank_map,
		MPI_cpu, MPI_cube, MPI_cpu_adj, MPI_cube_adj, MPI_direction, MPI_interface,
		neighbor_cpu_eq, Ncube_Ncpu_eq, 
		Scube_Ncpu_eq, Rcube_Ncpu_eq, send_data_curr_eq, recv_data_curr_eq, send_data_neig_eq, recv_data_neig_eq, Sdir_eq, Rdir_eq, 
		ist_eq,
		csl, 
		adj_number, 
		adjX_eq, adjY_eq, adjZ_eq,
		FWS);





	BCM_Initial_condition(myid, Ncube, N_wallcube,
		rank_map,
		FWS,
		U1_, 
		U1,  
		U1q);


	for (icube = 1; icube < Ncube; icube++) {  

#pragma omp parallel for private(j,k)

		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {
				for (k = 0; k <= nzzz; k++) {  

					U1_[icube][i][j][k][0] = U1[icube][i][j][k][0]; 
					U1_[icube][i][j][k][1] = U1[icube][i][j][k][1]; 
					U1_[icube][i][j][k][2] = U1[icube][i][j][k][2]; 
					U1_[icube][i][j][k][3] = U1[icube][i][j][k][3]; 
					U1_[icube][i][j][k][4] = U1[icube][i][j][k][4]; 

					U1q[icube][i][j][k][0] = U1[icube][i][j][k][0]; 
					U1q[icube][i][j][k][1] = U1[icube][i][j][k][1]; 
					U1q[icube][i][j][k][2] = U1[icube][i][j][k][2]; 
					U1q[icube][i][j][k][3] = U1[icube][i][j][k][3]; 
					U1q[icube][i][j][k][4] = U1[icube][i][j][k][4]; 
					
				}
			}
		}

	}

	
	// --------------------------------------------------------------------------------------- //
	// -------------------------------- detect the Ghost cell -------------------------------- //

	int ii,jj,kk;


	int (*tempFWS)[X_size][Y_size][Z_size] = new int[Ncube][X_size][Y_size][Z_size]; 


	for (icube = 1; icube < Ncube; icube++) {    
#pragma omp parallel for private(j,k)
		for (i = 0; i <= nxxx; i++) {
			for (j = 0; j <= nyyy; j++) {
				for (k = 0; k <= nzzz; k++) {  

					tempFWS[icube][i][j][k] = FWS[icube][i][j][k];;

				}
			}
		}
	}



	for (icube = 1; icube < Ncube; icube++) {  

		if (csl[icube] == 0) {
#pragma omp parallel for private(j,k)

			for (i = 1; i <= nxx; i++) {
				for (j = 1; j <= nyy; j++) {
					for (k = 1; k <= nzz; k++) { 

						if (FWS[icube][i][j][k] == ISOLID & 
							(FWS[icube][i+1][j  ][k  ] == IFLUID |
							FWS[icube][i-1][j  ][k  ] == IFLUID |
							FWS[icube][i  ][j-1][k  ] == IFLUID |
							FWS[icube][i  ][j+1][k  ] == IFLUID |
							FWS[icube][i  ][j  ][k-1] == IFLUID |
							FWS[icube][i  ][j  ][k+1] == IFLUID )
							)

							tempFWS[icube][i  ][j  ][k  ] = -2;  

					}
				}
			}

		}    // ---- if (csl[icube] == 0) ---- // 
	}




	// for (icube = 1; icube < Ncube; icube++) {  

	// if (csl[icube] == 0)
	// #pragma omp parallel for private(j,k)

	// for (i = 0; i <= nxxx; i++) {
	// for (j = 0; j <= nyyy; j++) {
	// for (k = 0; k <= nzzz; k++) { 

	// if (tempFWS[icube][i][j][k] == -2) 

	// FWS[icube][i][j][k] = IGHOST;



	// }
	// }
	// }

	// }    // ---- for (icube = 1; icube < Ncube; icube++) ---- //



	for (icube = 1; icube < Ncube; icube++) {  

		if (csl[icube] == 0)
#pragma omp parallel for private(j,k)

			for (i = 0; i <= nxxx; i++) {
				for (j = 0; j <= nyyy; j++) {
					for (k = 0; k <= nzzz; k++) { 


						if (tempFWS[icube][i][j][k] == -2) {

							//FWS[icube][i][j][k] = IGHOST;


							for (ii = i-0; ii <= i+0; ii++) {
								for (jj = j-0; jj <= j+0; jj++) {
									for (kk = k-0; kk <= k+0; kk++) { 

										if (ii < 0 | ii > nxxx | jj < 0 | jj > nyyy | kk < 0 | kk > nzzz)
											continue;

										FWS[icube][ii][jj][kk] = IGHOST;

									}
								}
							}

						}    // ---- if (tempFWS[icube][i][j][k] == -2) ---- //

					}
				}
			}

	}    // ---- for (icube = 1; icube < Ncube; icube++) ---- //



	delete[] tempFWS;

	// -------------------------------- detect the Ghost cell -------------------------------- //
	// --------------------------------------------------------------------------------------- //	


	
	
	BCM_Immersed_boundary(myid, Ncube, N_wallcube, &NBC, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, 
		cube_size, csl, Xcube, Ycube, Zcube, Xcnt, Ycnt, Zcnt, FWS, wallcube);				  
	
	


	BCM_Interface(myid,Ncube, 

		MPI_Nadj,

		Ncpu_bs, Ncpu_eq, Ncpu_sb,
		Max_nei_bs,Max_nei_eq,Max_nei_bs,

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
		U1_);


	// ---------------------------------------------------------- //
	// ---------------------------------------------------------- //


	char Number_of_BC[1024] = "Number of Boundary Cells";



	sprintf(file_name,"IBM_minus""%0.5d"".dat",myid);   
	fptr = fopen(file_name,"r"); 

	while( strcmp(str,Number_of_BC) != 0)  { fscanf(fptr,"%[^\n]\n",str); }

	fscanf(fptr,"%d\n",&NBC_minus);

	fclose(fptr);	


	printf("NBC_minus = %d\n",NBC_minus);



	NBC = NBC_minus;

	double (*weight_minus) = new double[NBC*8+1];
	int (*GCindex_minus) = new int[NBC*4+1];
	int (*IPsur_minus) = new int[NBC*4+1];

	mp_switch = -1;
	BCM_Reading_IBM(myid, &NBC, mp_switch, weight_minus, GCindex_minus, IPsur_minus);

	double (*bU1_minus) = new double[NBC_minus+1];
	double (*bU2_minus) = new double[NBC_minus+1];
	double (*bU3_minus) = new double[NBC_minus+1];
	double (*bU4_minus) = new double[NBC_minus+1];
	double (*bU5_minus) = new double[NBC_minus+1];

	double (*bU1_minus_store) = new double[NBC_minus+1];
	double (*bU2_minus_store) = new double[NBC_minus+1];
	double (*bU3_minus_store) = new double[NBC_minus+1];
	double (*bU4_minus_store) = new double[NBC_minus+1];
	double (*bU5_minus_store) = new double[NBC_minus+1];

	// =========================================================== //


	char strp[1024];

	sprintf(file_name,"IBM_plus""%0.5d"".dat",myid);   
	fptr = fopen(file_name,"r"); 

	while( strcmp(strp,Number_of_BC) != 0)  { fscanf(fptr,"%[^\n]\n",strp); }

	fscanf(fptr,"%d\n",&NBC_plus);


	fclose(fptr);	


	printf("NBC_plus = %d\n",NBC_plus);


	NBC = NBC_plus;


	double (*weight_plus) = new double[NBC*8+1];
	int (*GCindex_plus) = new int[NBC*4+1];
	int (*IPsur_plus) = new int[NBC*4+1];

	mp_switch = 1;
	BCM_Reading_IBM(myid, &NBC, mp_switch, weight_plus, GCindex_plus, IPsur_plus);

	double (*bU1_plus) = new double[NBC_plus+1];
	double (*bU2_plus) = new double[NBC_plus+1];
	double (*bU3_plus) = new double[NBC_plus+1];
	double (*bU4_plus) = new double[NBC_plus+1];
	double (*bU5_plus) = new double[NBC_plus+1];

	double (*bU1_plus_store) = new double[NBC_plus+1];
	double (*bU2_plus_store) = new double[NBC_plus+1];
	double (*bU3_plus_store) = new double[NBC_plus+1];
	double (*bU4_plus_store) = new double[NBC_plus+1];
	double (*bU5_plus_store) = new double[NBC_plus+1];

	NBC = max(NBC_plus,NBC_minus);




#pragma omp parallel for private(Ntemp,gicube,gi,gj,gk)
	for (int iNBC = 1; iNBC <= NBC_minus; iNBC++) {

		Ntemp = (iNBC-1)*4;

		gicube = GCindex_minus[Ntemp+1];
		gi = GCindex_minus[Ntemp+2];
		gj = GCindex_minus[Ntemp+3];
		gk = GCindex_minus[Ntemp+4];

		bU1_minus[iNBC] = U1_[gicube][gi][gj][gk][0];
		bU2_minus[iNBC] = U1_[gicube][gi][gj][gk][1];
		bU3_minus[iNBC] = U1_[gicube][gi][gj][gk][2];
		bU4_minus[iNBC] = U1_[gicube][gi][gj][gk][3];
		bU5_minus[iNBC] = U1_[gicube][gi][gj][gk][4];

	}


#pragma omp parallel for private(Ntemp,gicube,gi,gj,gk)
	for (int iNBC = 1; iNBC <= NBC_plus; iNBC++) {

		Ntemp = (iNBC-1)*4;

		gicube = GCindex_plus[Ntemp+1];
		gi = GCindex_plus[Ntemp+2];
		gj = GCindex_plus[Ntemp+3];
		gk = GCindex_plus[Ntemp+4];

		bU1_plus[iNBC] = U1_[gicube][gi][gj][gk][0];
		bU2_plus[iNBC] = U1_[gicube][gi][gj][gk][1];
		bU3_plus[iNBC] = U1_[gicube][gi][gj][gk][2];
		bU4_plus[iNBC] = U1_[gicube][gi][gj][gk][3];
		bU5_plus[iNBC] = U1_[gicube][gi][gj][gk][4];

	}



	// ---------------------------------------------------------- //
	// ---------------------------------------------------------- //





	// for (icube = 1; icube < Ncube; icube++) {  

	// #pragma omp parallel for private(j,k)
	// for (i = 0; i <= nxxx; i++) {
	// for (j = 0; j <= nyyy; j++) {
	// for (k = 0; k <= nzzz; k++) {  

	// if ( FWS[icube][i][j][k] != ISOLID)

	// FWS[icube][i][j][k] = IFLUID;

	// }
	// }
	// }

	// }



#pragma omp parallel for private(Ntemp,gicube,gi,gj,gk)
	for (int iNBC = 1; iNBC <= NBC_minus; iNBC++) {

		Ntemp = (iNBC-1)*4;

		gicube = GCindex_minus[Ntemp+1];
		gi = GCindex_minus[Ntemp+2];
		gj = GCindex_minus[Ntemp+3];
		gk = GCindex_minus[Ntemp+4];

		FWS[gicube][gi][gj][gk] = IGHOST;

	}


#pragma omp parallel for private(Ntemp,gicube,gi,gj,gk)
	for (int iNBC = 1; iNBC <= NBC_plus; iNBC++) {

		Ntemp = (iNBC-1)*4;

		gicube = GCindex_plus[Ntemp+1];
		gi = GCindex_plus[Ntemp+2];
		gj = GCindex_plus[Ntemp+3];
		gk = GCindex_plus[Ntemp+4];

		FWS[gicube][gi][gj][gk] = IGHOST;

	}











	BCM_Interface(myid,Ncube, 

		MPI_Nadj,

		Ncpu_bs, Ncpu_eq, Ncpu_sb,
		Max_nei_bs,Max_nei_eq,Max_nei_bs,

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
		U1_);


	BCM_FWS_Interface(myid,Ncube, 
		MPI_Nadj,
		Ncpu_eq, 
		Max_nei_eq,
		nadjX_eq, nadjY_eq, nadjZ_eq,
		rank_map,
		MPI_cpu, MPI_cube, MPI_cpu_adj, MPI_cube_adj, MPI_direction, MPI_interface,
		neighbor_cpu_eq, Ncube_Ncpu_eq, 
		Scube_Ncpu_eq, Rcube_Ncpu_eq, send_data_curr_eq, recv_data_curr_eq, send_data_neig_eq, recv_data_neig_eq, Sdir_eq, Rdir_eq, 
		ist_eq,
		csl, 
		adj_number, 
		adjX_eq, adjY_eq, adjZ_eq,
		FWS);



#pragma omp parallel for private(Ntemp,gicube,gi,gj,gk)
	for (int iNBC = 1; iNBC <= NBC_minus; iNBC++) {

		Ntemp = (iNBC-1)*4;

		gicube = GCindex_minus[Ntemp+1];
		gi = GCindex_minus[Ntemp+2];
		gj = GCindex_minus[Ntemp+3];
		gk = GCindex_minus[Ntemp+4];

		//FWS[gicube][gi][gj][gk] = -10;

		U1[gicube][gi][gj][gk][0] = rho0;
		U1[gicube][gi][gj][gk][1] = 0;
		U1[gicube][gi][gj][gk][2] = 0;
		U1[gicube][gi][gj][gk][3] = 0;
		U1[gicube][gi][gj][gk][4] = P0/(K-1);

		U1_[gicube][gi][gj][gk][0] = rho0;
		U1_[gicube][gi][gj][gk][1] = 0;
		U1_[gicube][gi][gj][gk][2] = 0;
		U1_[gicube][gi][gj][gk][3] = 0;
		U1_[gicube][gi][gj][gk][4] = P0/(K-1);

	}


#pragma omp parallel for private(Ntemp,gicube,gi,gj,gk)
	for (int iNBC = 1; iNBC <= NBC_plus; iNBC++) {

		Ntemp = (iNBC-1)*4;

		gicube = GCindex_plus[Ntemp+1];
		gi = GCindex_plus[Ntemp+2];
		gj = GCindex_plus[Ntemp+3];
		gk = GCindex_plus[Ntemp+4];

		//FWS[gicube][gi][gj][gk] = 10;


		U1[gicube][gi][gj][gk][0] = rho0;
		U1[gicube][gi][gj][gk][1] = 0;
		U1[gicube][gi][gj][gk][2] = 0;
		U1[gicube][gi][gj][gk][3] = 0;
		U1[gicube][gi][gj][gk][4] = P0/(K-1);

		U1_[gicube][gi][gj][gk][0] = rho0;
		U1_[gicube][gi][gj][gk][1] = 0;
		U1_[gicube][gi][gj][gk][2] = 0;
		U1_[gicube][gi][gj][gk][3] = 0;
		U1_[gicube][gi][gj][gk][4] = P0/(K-1);

	}





// =============================================== //
	for (step = 1 ; step <= count; step++) {       //
// =============================================== //



// ============================================================================ //
		for (int iteration = 1; iteration < 50000; iteration++) {               //
// ============================================================================ //		


	// ------------------------------------------ //
			for (int RK = 1; RK < 4; ++RK) {      //
	// ------------------------------------------ //


				BCM_X_boundary_condition(myid, NXbc_l, NXbc_u, Xbc_l, Xbc_u, U1_);

				BCM_Y_boundary_condition(NYbc_l, NYbc_u, Ybc_l, Ybc_u, U1_);

				BCM_Z_boundary_condition(NZbc_l, NZbc_u, Zbc_l, Zbc_u, U1_);

				
				for (int ig = 1; ig <= 1; ig++) {

					BCM_Ghostcell_minus(myid, &NBC_minus, weight_minus, GCindex_minus, IPsur_minus, FWS, U1_);

					BCM_Ghostcell_plus(myid, &NBC_plus, weight_plus, GCindex_plus, IPsur_plus, FWS, U1_);

				}
				

				BCM_Interface(myid,Ncube, 

					MPI_Nadj,

					Ncpu_bs, Ncpu_eq, Ncpu_sb,
					Max_nei_bs,Max_nei_eq,Max_nei_bs,

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
					U1_);
				
				
				BCM_Flux_XYZ_Viscous_Runge_kutta(myid, Ncube, RK, deltaT, deltaTau, e, FWS, csl, cube_size,
					U1_,U1 ,U1q,U1p1,U1p2,
					Rku1,Residual1,
					er);
					


				BCM_Interface(myid,Ncube, 

					MPI_Nadj,

					Ncpu_bs, Ncpu_eq, Ncpu_sb,
					Max_nei_bs,Max_nei_eq,Max_nei_bs,

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
					U1_);


		// ------------------------------------------ //
			}    //---- Runge-Kutaa end ---- //       //
		// ------------------------------------------ //


			E1 = max(E1, er[1]);
			E2 = max(E2, er[2]);
			E3 = max(E3, er[3]);
			E4 = max(E4, er[4]);
			E5 = max(E5, er[5]);

			er[1] = er[1]/E1;
			er[2] = er[2]/E2;
			er[3] = er[3]/E3;
			er[4] = er[4]/E4;
			er[5] = er[5]/E5;


			// if (myid == 0) {

			// printf("%d\t%f\t%f\t%f\t%f\t%f\n",iteration,er[1],er[2],er[3],er[4],er[5]);

			// }



// ==================================================================================================================== //
			if ((er[1]<0.02 & er[2]<0.02 & er[3]<0.02 & er[4]<0.02 & er[5]<0.02) | iteration == iteration_end_step) {   //
// ==================================================================================================================== //


				// if (step >= start_step) 
				// BCM_Statistic(myid, Ncube, step, start_step, statistic_step, dp_step, rank_map, U1_,U2_,U3_,U4_,U5_);



				//// ============================================================================== ////
				if (myid == 0) {

					//printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",iteration,step,er[1],er[2],er[3],er[4],er[5]);
					printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",iteration,-er[6],-er[7],er[1],er[2],er[3],er[4],er[5]);


				}
				//// ============================================================================= ////



				for (icube = 1; icube < Ncube; icube++) {  

#pragma omp parallel for private(j,k)
					for (i = n_buffer; i <= nx; i++) {
						for (j = n_buffer; j <= ny; j++) {
							for (k = n_buffer; k <= nz; k++) {  

								U1q[icube][i][j][k][0] = U1[icube][i][j][k][0]; 
								U1q[icube][i][j][k][1] = U1[icube][i][j][k][1]; 
								U1q[icube][i][j][k][2] = U1[icube][i][j][k][2]; 
								U1q[icube][i][j][k][3] = U1[icube][i][j][k][3]; 
								U1q[icube][i][j][k][4] = U1[icube][i][j][k][4];

								U1[icube][i][j][k][0] = U1_[icube][i][j][k][0]; 
								U1[icube][i][j][k][1] = U1_[icube][i][j][k][1]; 
								U1[icube][i][j][k][2] = U1_[icube][i][j][k][2]; 
								U1[icube][i][j][k][3] = U1_[icube][i][j][k][3]; 
								U1[icube][i][j][k][4] = U1_[icube][i][j][k][4];    

							}
						}
					}

				}

				break;  /**** jump out of the iteration ****/

// ==================================================================================================================== //
			}	/**** if-end ****/																						//
// ==================================================================================================================== //



// ============================================================================ //
		}    /**** iteration-end ****/ 											//
// ============================================================================ //

		/*
		for (icube = 1; icube < Ncube; icube++) {  

		#pragma omp parallel for private(j,k)
		for (i = n_buffer-1; i <= nxx; i++) {
		for (j = n_buffer-1; j <= nyy; j++) {
		for (k = n_buffer-1; k <= nzz; k++) {  

		U1_[icube][i][j][k][0] = FWS[icube][i][j][k]; 

		}
		}
		}

		}
		*/

		if ( step%output_step == 0 ) {

			BCM_Output(myid, Ncube, step, rank_map, U1_,U1q,cube_size, Xcnt, Ycnt, Zcnt);


		}




// =============================================== //
	}    /**** step-end ****/                      // 
// =============================================== //


	printf("\n--------------------------------------------\n");


	delete []MPI_cube;
	delete []MPI_cpu;
	delete []MPI_cube_adj;
	delete []MPI_cpu_adj;
	delete []MPI_direction;
	delete []MPI_interface;


	delete []neighbor_cpu_eq;
	delete []Ncube_Ncpu_eq;

	delete []Scube_Ncpu_eq;
	delete []Rcube_Ncpu_eq;


	delete []neighbor_cpu_sb;
	delete []Ncube_Ncpu_sb;

	delete []Scube_Ncpu_sb;
	delete []Rcube_Ncpu_sb;


	delete []neighbor_cpu_bs;
	delete []Ncube_Ncpu_bs;

	delete []Scube_Ncpu_bs;
	delete []Rcube_Ncpu_bs;

	delete []Xbc_u;
	delete []Xbc_l;

	delete []Ybc_u;
	delete []Ybc_l;

	delete []Zbc_u;
	delete []Zbc_l;


// ============================================ //
	MPI_Finalize();    /**** MPI-end ****/      //
// ============================================ //

}
