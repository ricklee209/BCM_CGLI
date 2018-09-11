// ---- PROE=0 : CGLI Roe scheme ---- //

// ---- PROE=1 : Weiss and Smith Roe scheme ---- //

// ---- APROE=2 : Roe modified by X. Li ---- //

// ---- LMROE=3 : Roe modified by F. Rieper ---- //

// ---- TDROE=4 : Roe modified by B. Thornber ---- //

// ---- ROEAMTP=5 : Roe modified by X. Li ---- //


#define ROE 0


// ---- limiter : limiter function ---- //
// ---- Nolimter : No limiter function ---- //

#define Nolimiter 





// ---- Kcomputer : K ---- //
// ---- PC : my desktop ---- //

#define PC


// ---- DT : dual time stepping ---- //
// ---- NODT : No dual time stepping ---- //

#define DT 


// ---- 0 : no output coarse result ---- //
// ---- 1 : output coarse result ---- //
#define CR 1


// ---- RK3NODT : 3rd Runge-Kutta ---- //
// ---- RKO5 : Bogey and Bailly 5th Runge-Kutta ---- //
// ---- DPLUSGS ---- //
// ---- LUSGS ---- //

#define LUSGS 


// ---- DTau : DT includes deltaTau ---- //
// ---- DTau_CFL : CFL_number for deltaTau ---- //
// ---- DTau_fix : DT includes deltaTau and deltaTau is fixed---- //
// ---- NODTau : DT doesn't include deltaTau (Aeroacoustic) ---- //
// ---- DTauCAA : Can run Aeroacoustic ---- //

#define NODTau
#define DTau_CFL 1000


// ---- ILES : implicit turbulence model ---- //
// ---- DNS : implicit turbulence model ---- //
// ---- NOUPD : No UpWinding term ---- //
  
#define NOUPD 


// ---- Wall_model : Wall model ---- //
// ---- No_Wall_model : Wall model ---- //

#define No_Wall_model


// ---- K_endian_transfer ---- //
// ---- No_endian_transfer ---- //

#define No_endian_transfer




// ---- Limiter function ---- //
#define rho_limit 0.8
#define rholimit  1.8

#define U_limit   -40.0
#define Ulimit    120.0

#define V_limit   -40.0
#define Vlimit    120.0

#define W_limit   -40.0
#define Wlimit    120.0

#define P_limit   101300*0.8
#define Plimit    101300*1.2




#define Char_D 0.1
