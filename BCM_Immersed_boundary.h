



#define CROSS(dest, v1, v2) \
          dest[0] = v1[1]*v2[2]-v1[2]*v2[1]; \
		  dest[1] = v1[2]*v2[0]-v1[0]*v2[2]; \
		  dest[2] = v1[0]*v2[1]-v1[1]*v2[0]; 

#define DOT(v1, v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest, v1, v2) \
           dest[0] = v1[0]-v2[0]; \
		   dest[1] = v1[1]-v2[1]; \
		   dest[2] = v1[2]-v2[2];


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

