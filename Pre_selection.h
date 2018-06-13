
// ---- PROE=1 : Weiss and Smith Roe scheme ---- //

// ---- APROE=2 : Roe modified by X. Li ---- //

// ---- LMROE=3 : Roe modified by F. Rieper ---- //

// ---- TDROE=4 : Roe modified by B. Thornber ---- //

// ---- ROEAMTP=5 : Roe modified by X. Li ---- //


#define ROE 3 


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
#define CR 0


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

#define DTauCAA 
#define DTau_CFL 1000


// ---- ILES : implicit turbulence model ---- //
// ---- DNS : implicit turbulence model ---- //

#define DNS 


// ---- Wall_model : Wall model ---- //
// ---- No_Wall_model : Wall model ---- //

#define No_Wall_model


// ---- K_endian_transfer ---- //
// ---- No_endian_transfer ---- //

#define No_endian_transfer




// ---- Limiter function ---- //
#define rho_limit 0.5
#define rholimit  1.5

#define U_limit   -100.0
#define Ulimit    200.0

#define V_limit   -100.0
#define Vlimit    200.0

#define W_limit   -100.0
#define Wlimit    200.0

#define P_limit   101300*0.5
#define Plimit    101300*2.0




#define Char_D 4.5e-5 
