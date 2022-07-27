#ifndef SELA_ADV1D_NON_PERIODIC_LAG
#define SELA_ADV1D_NON_PERIODIC_LAG

#include <math.h>             // function  floor
#include "parameter_reader.h" // type      PC_tree_t
                              // functions PC_get, PC_double, PC_int
#include "string_helpers.h"   // macro     ERROR_MESSAGE

typedef struct adv1d_non_periodic_lag_t {
    int d;
    double min;
    double max;
    int N;
    double v;
    double* lag; // Lagrange coefficients
    double* buf; // Stores f + d points on the left + d points on the right
                 // for butterfly method (extrapolation)
} adv1d_non_periodic_lag_t;

void adv1d_non_periodic_lag_init(adv1d_non_periodic_lag_t* *adv, PC_tree_t conf, double* x, int sizex, int mpi_rank) {
    long tmp;
    double val;
    if (PC_get(conf, ".adv1d_lag")) {
        if (PC_get(PC_get(conf,".adv1d_lag"), ".d")) {
            PC_int(PC_get(PC_get(conf,".adv1d_lag"), ".d"), &tmp);
        } else {
            ERROR_MESSAGE("#Error in advection %s: missing d.\n", conf->key);
        }
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing adv1d_lag.\n", conf->key);
    }
    if (PC_get(conf, ".v")) {
        PC_double(PC_get(conf,".v"), &val);
    } else {
        ERROR_MESSAGE("#Error in advection %s: missing v.\n", conf->key);
    }
    *adv = malloc(sizeof(adv1d_non_periodic_lag_t));
    int d = (int)tmp;
    (*adv)->d = d;
    (*adv)->min = x[0];
    (*adv)->max = x[sizex-1];
    (*adv)->N = sizex-1;
    (*adv)->v = val;
    (*adv)->lag = malloc((2*d+2)*sizeof(double));
    (*adv)->buf = malloc((sizex+2*d)*sizeof(double));
    if (mpi_rank == 0) {
        printf("#adv1d_lag:d=%d min=%1.20lg max=%1.20lg N=%d\n",(*adv)->d,
            (*adv)->min,(*adv)->max,(*adv)->N);
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
void adv1d_non_periodic_lag_compute_i0_and_alpha(double coeff, double dt, double xmin,
        double xmax, int N, int* i0, double* alpha) {
    double dx = (xmax - xmin) / ((double)N);
    double displacement = -coeff*dt/dx;
    *i0 = (int)floor(displacement);
    *alpha = displacement - (double)(*i0);
}

/*
 * Computes Lagrange coefficients.
 *
 * @param[in]  x displacement
 * @param[in]  d degree of the interpolator is 2d+1
 * @param[out] lag coefficients of the interpolator
 */
void adv1d_non_periodic_lag_compute_lag(double x, int d, double* lag) {
    int i;
    double a;
    
    if (d>0) {
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
    } else {
        lag[0] = 1. - x;
        lag[1] = x;
    }
}

void adv1d_non_periodic_lag_semi_lag_advect_classical(
        int N, int i0, double* lag, int d, double* buf, double* f_in_and_out) {
    int i, j, index;
    
    for (i = 0; i < N+1; i++)
        buf[i + d] = f_in_and_out[i];
    // Fill buffer on the left
    // TODO: butterfly
    for (i = 0; i < d; i++) {
        buf[i] = f_in_and_out[0];
    }
    // Fill buffer on the right
    // TODO: butterfly
    for (i = 0; i < d; i++) {
        buf[i + N+1+d] = f_in_and_out[N];
    }
    
    for (i = 0; i < N+1; i++) {
        f_in_and_out[i] = 0.;
        for (j = -d; j <= d+1; j++) {
            index = i+j+i0 + d;
            if (index < 0)
                index = 0;
            if (index > N+2*d)
                index = N+2*d;
            f_in_and_out[i] += lag[j+d] * buf[index];
        }
    }
}


void adv1d_non_periodic_lag_compute(adv1d_non_periodic_lag_t* adv,
        double* useless_parameter, double* f_in_and_out, double dt){
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
    
    adv1d_non_periodic_lag_compute_i0_and_alpha(v, dt, min, max, N, &i0, &alpha);
    adv1d_non_periodic_lag_compute_lag(alpha, d, adv->lag);
    adv1d_non_periodic_lag_semi_lag_advect_classical(N, i0, adv->lag, d, adv->buf, f_in_and_out);
}


#endif // ifndef SELA_ADV1D_NON_PERIODIC_LAG
