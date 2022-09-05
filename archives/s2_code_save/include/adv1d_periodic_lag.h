#ifndef SELA_ADV1D_PERIODIC_LAG
#define SELA_ADV1D_PERIODIC_LAG

#include <math.h>             // function  floor
#include "parameter_reader.h" // type      PC_tree_t
                              // functions PC_get, PC_double, PC_int
#include "string_helpers.h"   // macro     ERROR_MESSAGE

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
void adv1d_periodic_lag_compute_i0_and_alpha(double coeff, double dt, double xmin,
        double xmax, int N, int* i0, double* alpha) {
    
    *alpha = -coeff*dt/(xmax-xmin);
    *alpha = *alpha-floor(*alpha);
    *alpha *= (double)N;
    *i0 = (int)floor(*alpha);
    if (*i0 == N) {
        *alpha = 0.;
        *i0 = 0;
    }
    *alpha = *alpha-((double)*i0);
}

void adv1d_periodic_lag_semi_lag_advect_classical(double* buf, int N, int i0, double* lag, int d, double* f) {
    int i, j;
    
    for (i = 0; i < N+1; i++)
        buf[i] = f[i];
    
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
    
    adv1d_periodic_lag_compute_i0_and_alpha(v, dt, min, max, N, &i0, &alpha);
    adv1d_periodic_lag_compute_lag(alpha, d, adv->lag);
    adv1d_periodic_lag_semi_lag_advect_classical(adv->buf, N, i0, adv->lag, d, f);
}


#endif // ifndef SELA_ADV1D_PERIODIC_LAG
