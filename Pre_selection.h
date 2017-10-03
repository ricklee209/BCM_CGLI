
// ---- PROE=1 : Weiss and Smith Roe scheme ---- //

// ---- APROE=2 : Roe modified by X. Li ---- //

// ---- LMROE=3 : Roe modified by F. Rieper ---- //

// ---- TDROE=4 : Roe modified by B. Thornber ---- //

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

// ---- RK3NODT : 3rd Runge-Kutta ---- //
// ---- RKO5 : Bogey and Bailly 5th Runge-Kutta ---- //
// ---- DPLUSGS ---- //
// ---- LUSGS ---- //

#define LUSGS 


// ---- ILES : implicit turbulence model ---- //
// ---- DNS : implicit turbulence model ---- //

#define DNS 


// ---- Wall_model : Wall model ---- //
// ---- No_Wall_model : Wall model ---- //

#define No_Wall_model


// ---- K_endian_transfer ---- //
// ---- No_endian_transfer ---- //

#define No_endian_transfer


#define Char_D 0.1
