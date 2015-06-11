



#define CROSS(dest, v1, v2) \
	dest[0] = v1[1]*v2[2]-v1[2]*v2[1]; \
	dest[1] = v1[2]*v2[0]-v1[0]*v2[2]; \
	dest[2] = v1[0]*v2[1]-v1[1]*v2[0]; 

#define DOT(v1, v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest, v1, v2) \
	dest[0] = v1[0]-v2[0]; \
	dest[1] = v1[1]-v2[1]; \
	dest[2] = v1[2]-v2[2];


// -----------------------------------------------------------------------------------//
// --------------------------- detect riganle Cube overlap ---------------------------//



#define X 0
#define Y 1
#define Z 2


#define FINDMINMAX(x0,x1,x2,min,max) \
	min = max = x0;   \
	if(x1<min) min=x1;\
	if(x1>max) max=x1;\
	if(x2<min) min=x2;\
	if(x2>max) max=x2;


/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
	if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;


#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;


/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
	if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;


#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return 0;


/*======================== Z-tests ========================*/



#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
	if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;


#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
	if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;


int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);

int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3]);


// --------------------------- detect riganle Cube overlap ---------------------------//
// -----------------------------------------------------------------------------------//





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
);


void BI_detect_iflag0
	(

	// ======================= //
	int myid,

	int pBI,

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
	);


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
	);


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
	);

void GetInverseMatrix( const double* inMat, double* outMat, const int n );

void swap(double *a, double *b);


void Dirichlet_Vandermonde_matrix(
// ================================================================================ //	
double dx,
double dy,
double dz,
double van[64],
double BIvan[8]
// ================================================================================ //	
);


void Neumann_Vandermonde_matrix(
// ================================================================================ //	
double dx,
double dy,
double dz,
double van[64],
double BIvan[8]
// ================================================================================ //	
);



