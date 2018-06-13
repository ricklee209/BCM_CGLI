



#define K 1.4
#define R 287 

#define Cv 717.5                      /**** specific constant volume heat capacity J/(kg.k) ****/
#define mu_L 0.0000185
#define Pr_L 0.72
#define R 287

double Cp = Cv*K;                     /**** Cp=Cv*K ****/
double lambda_L = mu_L*Cv*K/Pr_L;


/**** computational parameters ****/

#define Ep 1.

#define MAX 100000000
#define MIN -100000000

#define minimum 1.0e-12

double psi = 0.0001;

/**** computational parameters-end ****/

double rho0 = 1.1842;

double U0 = 69.0;

double Uin = 69.0;

double V0 = 0.;

double W0 = 0;

double P0 = 101300;

double T0 = P0/rho0/R;

