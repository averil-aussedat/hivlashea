#ifndef SELA_VP_1D1V_CART_ADVECT
#define SELA_VP_1D1V_CART_ADVECT

#include <math.h>             // function  floor
// #include "adv1d_non_periodic_lag.h" // type      adv1d_non_periodic_lag_t
//                                     // functions adv1d_non_periodic_lag_compute, 
//                                     //           adv1d_non_periodic_lag_init
//                                     // (used not directly in this file but
//                                     //  used by simulations that include this file)
// #include "adv1d_periodic_lag.h" 	// type      adv1d_periodic_lag_t
//                                     // functions adv1d_periodic_lag_compute, 
//                                     //			 adv1d_periodic_init
//                                     // (used not directly in this file but
//                                     //  used by simulations that include this file)

#include "remap.h"                  // type      parallel_stuff
                                    // functions exchange_parallelizations, local_to_global_2d
#include "parameter_reader.h"       // type      PC_tree_t
                                    // functions PC_get, PC_double, PC_int
#include "string_helpers.h"         // macro     ERROR_MESSAGE

#define adv1d_x_t adv1d_lag_t
#define adv1d_x_init adv1d_init

#define adv1d_v_t adv1d_lag_t
#define adv1d_v_init adv1d_init

// --- Table of content ---
// Useful little functions
// Periodic world
// Non-periodic world
// Common world

/////////////////////////////////////////////////////////////////
// Useful little functions
/////////////////////////////////////////////////////////////////

/*
 * Computes Lagrange coefficients.
 *
 * @param[in]  x displacement
 * @param[in]  d degree of the interpolator is 2d+1
 * @param[out] lag coefficients of the interpolator (expected size : 2d+2)
 */
void adv1d_compute_lag(double x, int d, double* lag) {
    int i;
    double a;
    
    if (d>0) { // if more than linear interpolation
        a = 1.;
        // Because of parity : (x-i)*(x+i) = x^2-i^2
        for (i = 2; i <= d; i++)
            a *= (x*x-((double)i)*((double)i))/(((double)d)*((double)d));
        a *= (x+1.)/((double)d);
        a *= (x-((double)d)-1.)/((double)d);
        lag[d]   = a*(x-1.)/((double)d);
        lag[d+1] = a*x/((double)d);
        a *= x*(x-1.)/(((double)d)*((double)d));
        for (i = -d; i <= -1; i++)
            lag[i+d] = a/((x-(double)i)/((double)d));
        for (i = 2; i <= d+1; i++)
            lag[i+d] = a/((x-(double)i)/((double)d));
        a = 1.;
        for (i=-d; i <= d+1; i++) {
            lag[i+d] *= a;
            a *= (double)d/((double)(d+i+1));
        }
        a = 1.;
        for (i = d+1; i >= -d; i--) {
            lag[i+d] *= a;
            a *= (double)d/((double)(i-1-d-1));
        }
    } else { // if linear interpolation
        lag[0] = 1. - x;
        lag[1] = x;
    }
}

/*
 * Computes the location where the advection makes us go, and converts
 * it into integer part and fractional part.
 *
 * @param[in] coeff  : advection parameter.
 * @param[in] dt     : time step.
 * @param[in] xmin   : minimum value of the mesh.
 * @param[in] xmax   : maximum value of the mesh.
 * @param[in] N      : number of cells in the mesh.
 * @param[out] i0    : integer part.
 * @param[out] alpha : fractional part.
 */
void adv1d_compute_i0_and_alpha(double coeff, double dt, double xmin,
        double xmax, int N, int* i0, double* alpha) {
    double dx = (xmax - xmin) / ((double)N); // N is the number of intervals
    double displacement = -coeff*dt/dx; // coeff is the advection speed
    *i0 = (int)floor(displacement);
    *alpha = displacement - (double)(*i0);
}

/////////////////////////////////////////////////////////////////
// Periodic world
/////////////////////////////////////////////////////////////////

typedef struct adv1d_periodic_lag_t {
    int d;
    double min;
    double max;
    int N;
    double v;
    double* lag; // Lagrange coefficients
    double* buf; // Stores f points
} adv1d_periodic_lag_t;

void adv1d_periodic_lag_init(adv1d_periodic_lag_t* *adv, PC_tree_t conf, double* x, int sizex, int mpi_rank) {
    long tmp;
    double val;
    if (PC_get(conf, ".adv1d_periodic_lag")) {
        if (PC_get(PC_get(conf,".adv1d_periodic_lag"), ".d")) {
            PC_int(PC_get(PC_get(conf,".adv1d_periodic_lag"), ".d"), &tmp);
        } else {
            ERROR_MESSAGE("#Error in advection %s: missing d.\n", conf->key);
        }
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing adv1d_periodic_lag.\n", conf->key);
    }
    if (PC_get(conf, ".v")) {
        PC_double(PC_get(conf,".v"), &val);
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing v.\n", conf->key);
    }
    *adv = malloc(sizeof(adv1d_periodic_lag_t));
    int d = (int)tmp;
    (*adv)->d = d;
    (*adv)->min = x[0];
    (*adv)->max = x[sizex-1];
    (*adv)->N = sizex-1;
    (*adv)->v = val;
    (*adv)->lag = malloc((2*d+2)*sizeof(double));
    (*adv)->buf = malloc((sizex)*sizeof(double));
    if (mpi_rank == 0) {
        printf("#adv1d_periodic_lag: d=%d min=%1.20lg max=%1.20lg N=%d v=%1.20lg\n",(*adv)->d,
            (*adv)->min,(*adv)->max,(*adv)->N,(*adv)->v);
    }
}

void adv1d_periodic_lag_semi_lag_advect_classical(double* buf, int N, int i0, double* lag, int d, double* f) {
    int i, j;
    for (i = 0; i < N+1; i++) {
        buf[i] = f[i];
    }    
    for (i = 0; i < N+1; i++) {
        f[i] = 0.;
        for (j = -d; j <= d+1; j++) {
            f[i] += lag[j+d] * buf[(i+j+i0+N)%N];
        }
    }
}

void adv1d_periodic_lag_compute(adv1d_periodic_lag_t* adv, double* f, double dt){
    int d;
    double min;
    double max;
    int N;
    double v;
    int i0;
    double alpha;
    
    d = adv->d;
    N = adv->N;
    min = adv->min;
    max = adv->max;
    v = adv->v;
    //printf("Lag d=%d  N=%d v=%1.20lg\n",adv->d,adv->N,adv->v);
    
    // adv1d_periodic_lag_compute_i0_and_alpha(v, dt, min, max, N, &i0, &alpha);
    // adv1d_periodic_lag_compute_lag(alpha, d, adv->lag);
    adv1d_compute_i0_and_alpha(v, dt, min, max, N, &i0, &alpha);
    adv1d_compute_lag(alpha, d, adv->lag);
    adv1d_periodic_lag_semi_lag_advect_classical(adv->buf, N, i0, adv->lag, d, f);
}

/////////////////////////////////////////////////////////////////
// Non-periodic world
/////////////////////////////////////////////////////////////////

typedef struct adv1d_non_periodic_lag_t {
    int d; // interpolation stencil is of width 2d+2
    int N; // number of intervals of the 1D mesh
    int kb; // order of the polynomial extrapolation (if needed)
    double min; // lower bound of the mesh (included)
    double max; // upper bound of the mesh (included)
    double v; // speed factor
    double* lag; // Lagrange coefficients
    double* buf; // Stores f + d points on the left + d points on the right
                 // for butterfly method (extrapolation)
    double* extrap; // coefficients of the extrapolation (size kb+1) 
                    // u_{J+1} := sum_{i=0}^{kb} extrap[i] * u_{J-kb+i}
} adv1d_non_periodic_lag_t;

/*
 * Initializes the coefficients for polynomial extrapolation.
 *
 * @param[in] kb : order of the polynomial
 * @param[out] extrap : array OF SIZE kb+1, ALREADY ALLOCATED
 */
void init_extrap (int kb, double* extrap) {
    int i,j; 
    // computation of the extrapolation coefficients
    // example with kb=3 : 
    //           0    1    2    3   <- content of extrap
    //   init    x    x    x   1.0    no call to j
    //    i=2    x    x   1.0  2.0    no call to j
    //    i=1    x   1.0  3.0  3.0    j from 2 to 2
    //    i=0   1.0  4.0  6.0  4.0    j from 1 to 2
    // then add sign to have -1.0 4.0 -6.0 4.0 (last coeff always positive).
    extrap[kb] = 1.0;
    for (i=kb-1; i>=0; --i) { // Pascal triangle 
        extrap[i] = 1.0; // adding new 1.0
        for (j=i+1; j<kb; ++j) { // internal Pascal coeff
            extrap[j] += extrap[j+1];
        }
        extrap[kb] += 1.0; // last Pascal coeff ("external 1")
        // printf("i=%d : ", i);
        // for (j=0; j<=kb; ++j) {
        //     printf("%f ", extrap[j]);
        // }
        // printf("\n");
    }
    for (i=1; i<=kb; i+=2) { // add sign.
        extrap[kb-i] *= -1.0; 
    }

    // printf("end : ");
    // for (j=0; j<=kb; ++j) {
    //     printf("%f ", extrap[j]);
    // }
    // printf("\n");
}

void adv1d_non_periodic_lag_init(adv1d_non_periodic_lag_t* *adv, PC_tree_t conf, double* x, int sizex, int mpi_rank) {
    long tmp_d, tmp_kb;
    double val;
    int i;
    if (PC_get(conf, ".adv1d_non_periodic_lag")) {
        if (PC_get(PC_get(conf,".adv1d_non_periodic_lag"), ".d")) {
            PC_int(PC_get(PC_get(conf,".adv1d_non_periodic_lag"), ".d"), &tmp_d);
        } else {
            ERROR_MESSAGE("#Error in advection %s: missing d.\n", conf->key);
        }
        if (PC_get(PC_get(conf,".adv1d_non_periodic_lag"), ".kb")) {
            PC_int(PC_get(PC_get(conf,".adv1d_non_periodic_lag"), ".kb"), &tmp_kb);
        } else {
            ERROR_MESSAGE("#Error in advection %s: missing kb.\n", conf->key);
        }
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing adv1d_non_periodic_lag.\n", conf->key);
    }
    if (PC_get(conf, ".v")) {
        PC_double(PC_get(conf,".v"), &val);
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing v.\n", conf->key);
    }
    *adv = malloc(sizeof(adv1d_non_periodic_lag_t));
    // int d = (int)tmp_d;
    // (*adv)->d = d;
    (*adv)->d = (int)tmp_d;
    (*adv)->kb = (int)tmp_kb;
    (*adv)->min = x[0];
    (*adv)->max = x[sizex-1];
    (*adv)->N = sizex-1;
    (*adv)->v = val;
    (*adv)->lag = malloc((2*(*adv)->d+2)*sizeof(double));
    (*adv)->buf = malloc((sizex+1+2*(*adv)->d)*sizeof(double)); // at least one outside value
    (*adv)->extrap = malloc(((*adv)->kb+1)*sizeof(double));
    init_extrap ((*adv)->kb, (*adv)->extrap);
    if (mpi_rank == 0) {
        printf("#adv1d_non_periodic_lag:d=%d min=%1.20lg max=%1.20lg N=%d v=%1.20lg kb=%d\n",(*adv)->d,
            (*adv)->min,(*adv)->max,(*adv)->N,(*adv)->v, (*adv)->kb);
        // printf("#extrap : ");
        for (i=0; i<=(*adv)->kb; ++i) {
            printf("%f ", (*adv)->extrap[i]);
        }
        printf("\n");
    }
}

/*
 * Computes the advection with polynomial extension of 
 * the boundary condition. 
 *
 * @param[in] coeff the speed of the advection
 * @param[in] N the number of intervals
 * @param[in] i0 (int) : index displacement of the advection
 * @param[in] lag (double*) : array of size 2d+2, coefficients of Lagrange polynomials
 * @param[in] d (int) : stencil width of the interpolation (2d+2 points)
 * @param[in] kb (int) : order of polynomial extension (0 constant, 1 linear...)
 * @param[in] extrap (double*) : size kb+1, recursion coefficients (kb=0 -> [1.0], kb=1 -> [-1.,2.], ...)
 * @param[out] buf : already-initialized temporary variable, same size as f_in_and_out
 * @param[out] f_in_and_out : array of size N+2d
 *
 */
void adv1d_non_periodic_lag_semi_lag_advect_classical (
        double coeff, int N, int i0, double* lag, int d, 
        int kb, double* extrap,
        double* buf, double* f_in_and_out) {
    int i=0, j=0, maxi=0, mini=0;

    // printf("i0 : %d, N : %d, d=%d, coeff=%f\n", i0, N, d, coeff);

    /*
        f_i&o  is of size N+1, index from [0 to N]
        Buffer is of size N+2+2d, index from [0 to N+2d+1]
     */

    // Shifted buffer. i runs along buf[.]
    if (coeff >= 0) { // advection towards the right : inflow at left, i0 negative
        mini = 0;
        maxi = imin(N+2*d+1,d-i0-1);
        // printf("inflow : going from %d to %d\n", mini, maxi);
        for (i=mini; i<=maxi; i++) { // filling the inflow condition (0.0)
            buf[i] = 0.0;
        }
        mini = maxi+1;
        maxi = N+d+imin(-i0,d+1); // respects the size of f and buf. min(N+d-i0, N+1+2d)
        // printf("middle : going from %d to %d,  f from %d to %d\n", mini, maxi, mini-d+i0, maxi-d+i0);
        for (i=mini; i<=maxi; i++) { // values of f^{n}
            buf[i] = f_in_and_out[i-(d-i0)];
        }
        mini = maxi+1;
        maxi = N+2*d+1;
        // printf("end : going from %d to %d\n", mini, maxi);
        for (i=mini; i<=maxi; i++) { // remaining terms of the buffer
            buf[i] = 0.0;
            for (j=0; j<kb+1; ++j) { // polynomial extrapolation of order kb 
                // printf("extrap i=%d, buf[%d]\n", i, i-kb-1+j);
                buf[i] += extrap[j] * buf[i-kb-1+j];
            }
            // printf("Setting buf[%d] to %f\n", i, buf[i]);
        }
    } else if (coeff < 0) { // advection towards the left : inflow at right, i0 positive 
        maxi = N+2*d+1; // last index of buf
        mini = imax(0,maxi-(d+i0));
        // printf("inflow : going from %d to %d\n", maxi, mini);
        for (i=maxi; i>=mini; --i) { // filling the infow condition (0.0)
            buf[i] = 0.0;
        }
        maxi = mini-1;
        mini = imax(0,d-i0); // taking as many f values as possible
        // printf("middle : going from %d to %d, f from %d to %d\n", maxi, mini,maxi-d+i0, mini-d+i0);
        for (i=maxi; i>=mini; --i) { // values of f^{n}
            buf[i] = f_in_and_out[i-d+i0]; // starts from f[N] and decreases. N+i-(maxi-(d+i0)-1)
        }
        maxi = mini-1;
        mini = 0;
        // printf("end : going from %d to %d\n", maxi, mini);
        for (i=maxi; i>=mini; --i) { // remaining terms of the buffer
            buf[i] = 0.0;
            for (j=0; j<kb+1; ++j) { // polynomial extrapolation of order kb 
                // printf("extrap i=%d, buf[%d]\n", i, i+j+1);
                buf[i] += extrap[kb-j] * buf[i+j+1];
            }
            // printf("Setting buf[%d] to %f\n", i, buf[i]);
        }
    }
    
    if (coeff>=0) {
        f_in_and_out[0] = 0.;
        mini=1; maxi=N;
    } else {
        f_in_and_out[N] = 0.;
        mini=0; maxi=N-1;
    }
    // mini=0; maxi=N;
    // advection part : Lagrange interpolation 
    // printf("f : ");
    for (i = mini; i <= maxi; i++) { // for each internal dof
        f_in_and_out[i] = 0.;
        for (j = -d; j <= d+1; j++) { // interpolation
            f_in_and_out[i] += lag[j+d] * buf[i+d+j];
        } // interpolation
        // printf("Setting i=%d with buf values from %d to %d\n", i, i, i+2*d+1);
        // printf("%5.2e\t", f_in_and_out[i]);
    } // for each internal dof
    // printf("\n");
}

/*
    TODO : 
    prendre du degré plus faible en arrivant au bord
    pour que le stencil reste à l'intérieur (choisir le max en fonction de d)
 */


/*
 * Computes the advection with polynomial extension of 
 * the boundary condition. 
 *
 * @param[in] coeff the speed of the advection
 * @param[in] N the number of intervals
 * @param[in] i0 (int) : index displacement of the advection
 * @param[in] lag (double*) : array of size 2d+2, coefficients of Lagrange polynomials
 * @param[in] d (int) : stencil width of the interpolation (2d+2 points)
 * @param[in] kb (int) : order of polynomial extension (0 constant, 1 linear...)
 * @param[in] extrap (double*) : size kb+1, recursion coefficients (kb=0 -> [1.0], kb=1 -> [-1.,2.], ...)
 * @param[out] buf : already-initialized temporary variable, same size as f_in_and_out
 * @param[out] f_in_and_out : array of size N+2d
 *
 */
void adv1d_non_periodic_lag_semi_lag_advect_classical_invert (
        double coeff, int N, int i0, double* lag, int d, 
        int kb, double* extrap,
        double* buf, double* f_in_and_out) {
    int i=0, j=0, maxi=0, mini=0;
    bool invert = (coeff < 0.0);

    double* f_right = (double*) malloc ((N+1)*sizeof(double));
    if (invert) {
        for (i=0; i<=N; ++i) {
            f_right [i] = f_in_and_out [N-i];
        }
        // i0 = -i0-1;
        // alpha = 1.0-alpha;        
    } else {
        for (i=0; i<=N; ++i) {
            f_right [i] = f_in_and_out [i];
        }
    }

    mini = 0;
    maxi = imin(N+2*d+1,d-i0-1);
    // printf("inflow : going from %d to %d\n", mini, maxi);
    for (i=mini; i<=maxi; i++) { // filling the inflow condition (0.0)
        buf[i] = 0.0;
    }
    mini = maxi+1;
    maxi = N+d+imin(-i0,d+1); // respects the size of f and buf. min(N+d-i0, N+1+2d)
    // printf("middle : going from %d to %d,  f from %d to %d\n", mini, maxi, mini-d+i0, maxi-d+i0);
    for (i=mini; i<=maxi; i++) { // values of f^{n}
        buf[i] = f_right[i-(d-i0)];
    }
    mini = maxi+1;
    maxi = N+2*d+1;
    // printf("end : going from %d to %d\n", mini, maxi);
    for (i=maxi; i<=maxi; i++) { // remaining terms of the buffer
        buf[i] = 0.0;
        for (j=0; j<kb+1; ++j) { // polynomial extrapolation of order kb 
            // printf("extrap i=%d, buf[%d]\n", i, i-kb-1+j);
            buf[i] += extrap[j] * buf[i-kb-1+j];
        }
        // printf("Setting buf[%d] to %f\n", i, buf[i]);
    }

    // advection part : Lagrange interpolation 
    for (i = 0; i <= N; i++) { // for each internal dof
        f_right[i] = 0.;
        for (j = -d; j <= d+1; j++) { // interpolation
            f_right[i] += lag[j+d] * buf[i+d+j];
        } // interpolation
        // printf("Setting i=%d with buf values from %d to %d\n", i, i, i+2*d+1);
    } // for each internal dof
    
    if (invert) {
        for (i=0; i<=N; ++i) {
            f_in_and_out [i] = f_right [N-i];
        }
    } else {
        for (i=0; i<=N; ++i) {
            f_in_and_out [i] = f_right [i];
        }
    }
    free (f_right);
}

/*
 * Elementary advection of the function f
 * following the parameters of adv on time dt
 *
 * @param[in] adv : advection parameters
 * @param[out] f_in_and_out : values to advect in-place
 * @param[out] dt : time step
 */
void adv1d_non_periodic_lag_compute(adv1d_non_periodic_lag_t* adv,
        double* f_in_and_out, double dt){
    int d;
    double min;
    double max;
    int N;
    double v;
    int i0;
    int kb;
    double alpha;
    
    d = adv->d;
    N = adv->N;
    min = adv->min;
    max = adv->max;
    v = adv->v;
    kb = adv->kb;
    adv1d_compute_i0_and_alpha(v, dt, min, max, N, &i0, &alpha);
    // printf("Lag d=%d  N=%d v=%1.20lg alpha=%f i0=%d\n",adv->d,adv->N,adv->v, alpha, i0);
    adv1d_compute_lag(alpha, d, adv->lag);
    adv1d_non_periodic_lag_semi_lag_advect_classical(v, N, i0, adv->lag, d, kb, adv->extrap, adv->buf, f_in_and_out);
}


/////////////////////////////////////////////////////////////////
// Common world
/////////////////////////////////////////////////////////////////

typedef struct adv1d_lag_t {
    adv1d_periodic_lag_t* periodic_adv;
    adv1d_non_periodic_lag_t* non_periodic_adv;
} adv1d_lag_t;

void adv1d_init(adv1d_lag_t* *adv, PC_tree_t conf, double* x, int sizex, int mpi_rank) {
    *adv = malloc(sizeof(adv1d_lag_t));
    (*adv)->periodic_adv = (void*)0;
    (*adv)->non_periodic_adv = (void*)0;
    if (PC_get(conf, ".adv1d_periodic_lag")) {
        adv1d_periodic_lag_init(&((*adv)->periodic_adv), conf, x, sizex, mpi_rank);
    } else if (PC_get(conf, ".adv1d_non_periodic_lag")) {
        adv1d_non_periodic_lag_init(&((*adv)->non_periodic_adv), conf, x, sizex, mpi_rank);
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing adv1d_periodic_lag or adv1d_non_periodic_lag.\n", conf->key);
    }
}

int advection_v(parallel_stuff* par_variables, adv1d_v_t* adv, double dt, double* E){
    int i_x;
    int i_v;
    double tmp;
    // int ihello=0;
    if (!par_variables->is_par_x) { // align correctly
        exchange_parallelizations(par_variables);
    }
    local_to_global_2d(par_variables, 0, 0); // modifies par_variables
    if (adv->periodic_adv) {
        adv1d_periodic_lag_t* chosen_adv = adv->periodic_adv;
        tmp = chosen_adv->v;
        for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
            for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++) {
                par_variables->f_1d[i_v] = par_variables->f_parallel_in_x[i_x][i_v];
            }
            chosen_adv->v = tmp*E[i_x+par_variables->global_indices[0]];
            adv1d_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
            for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++) {
                par_variables->f_parallel_in_x[i_x][i_v] = par_variables->f_1d[i_v];
            }
        }
        chosen_adv->v = tmp;
    } else {
        // printf("\t\thello %d\n", ++ihello);
        adv1d_non_periodic_lag_t* chosen_adv = adv->non_periodic_adv;
        // printf("\t\thello %d\n", ++ihello);
        tmp = chosen_adv->v; // save the speed factor
        // printf("\t\thello %d\n", ++ihello);
        for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) { // chosen subdomain in x
            // printf("\t\tix=%d, d = %d\n", i_x, chosen_adv->d);
            for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++) { // temporary container
                par_variables->f_1d[i_v] = par_variables->f_parallel_in_x[i_x][i_v];
            }
            // printf("\t\tix=%d 2, d = %d\n", i_x, chosen_adv->d);
            chosen_adv->v = tmp*E[i_x+par_variables->global_indices[0]]; // speed = electric field 
            // printf("\t\tix=%d 3, d = %d\n", i_x, chosen_adv->d);
            adv1d_non_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt); // advection along one dimension in v
            // printf("\t\tix=%d 4, d = %d\n", i_x, chosen_adv->d);
            for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++) { // retrieve the solution
                par_variables->f_parallel_in_x[i_x][i_v] = par_variables->f_1d[i_v];
            }
            // printf("\t\tix=%d 5, d = %d\n", i_x, chosen_adv->d);
        }
        chosen_adv->v = tmp; // reset the speed factor
    }
    return 0;
}

void advection_x(parallel_stuff* par_variables, adv1d_x_t* adv, double dt, double* v){
    int i_x;
    int i_v;
    double tmp;
    if (par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    local_to_global_2d(par_variables, 0, 0);
    if (adv->periodic_adv) {
        adv1d_periodic_lag_t* chosen_adv = adv->periodic_adv;
        tmp = chosen_adv->v;
        for (i_v = 0; i_v < par_variables->size_v_par_v; i_v++){
            for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++) {
                par_variables->f_1d[i_x] = par_variables->f_parallel_in_v[i_x][i_v];
            }
            chosen_adv->v = tmp*v[i_v+par_variables->global_indices[1]];
            adv1d_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
            for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++) {
                par_variables->f_parallel_in_v[i_x][i_v] = par_variables->f_1d[i_x];
            }
        chosen_adv->v = tmp;
        }
    } else {
        adv1d_non_periodic_lag_t* chosen_adv = adv->non_periodic_adv;
        tmp = chosen_adv->v;
        for (i_v = 0; i_v < par_variables->size_v_par_v; i_v++){
            for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++) {
                par_variables->f_1d[i_x] = par_variables->f_parallel_in_v[i_x][i_v];
            }
            chosen_adv->v = tmp*v[i_v+par_variables->global_indices[1]];
            adv1d_non_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
            for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++) {
                par_variables->f_parallel_in_v[i_x][i_v] = par_variables->f_1d[i_x];
            }
        }
        chosen_adv->v = tmp;
    }
}

#endif // ifndef SELA_VP_1D1V_CART_ADVECT
