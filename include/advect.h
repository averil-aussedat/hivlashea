#ifndef SELA_VP_1D1V_CART_ADVECT
#define SELA_VP_1D1V_CART_ADVECT

#include "adv1d_non_periodic_lag.h" // type      adv1d_non_periodic_lag_t
                                    // functions adv1d_non_periodic_lag_compute, 
                                    //           adv1d_non_periodic_lag_init
                                    // (used not directly in this file but
                                    //  used by simulations that include this file)
#include "adv1d_periodic_lag.h" 	// type      adv1d_periodic_lag_t
                                    // functions adv1d_periodic_lag_compute, 
                                    //			 adv1d_periodic_init
                                    // (used not directly in this file but
                                    //  used by simulations that include this file)
#include "remap.h"                  // type      parallel_stuff
                                    // functions exchange_parallelizations, local_to_global_2d

#define adv1d_x_t adv1d_lag_t
#define adv1d_x_init adv1d_init

#define adv1d_v_t adv1d_lag_t
#define adv1d_v_init adv1d_init

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
    if (!par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    local_to_global_2d(par_variables, 0, 0);
if (adv->periodic_adv) {
    adv1d_periodic_lag_t* chosen_adv = adv->periodic_adv;
    tmp = chosen_adv->v;
    for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++)
            par_variables->f_1d[i_v] = par_variables->f_parallel_in_x[i_x][i_v];
        chosen_adv->v = tmp*E[i_x+par_variables->global_indices[0]];
        adv1d_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++)
            par_variables->f_parallel_in_x[i_x][i_v] = par_variables->f_1d[i_v];
    }
    chosen_adv->v = tmp;
} else {
    adv1d_non_periodic_lag_t* chosen_adv = adv->non_periodic_adv;
    tmp = chosen_adv->v;
    for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++)
            par_variables->f_1d[i_v] = par_variables->f_parallel_in_x[i_x][i_v];
        chosen_adv->v = tmp*E[i_x+par_variables->global_indices[0]];
        adv1d_non_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++)
            par_variables->f_parallel_in_x[i_x][i_v] = par_variables->f_1d[i_v];
    }
    chosen_adv->v = tmp;
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
        for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++)
            par_variables->f_1d[i_x] = par_variables->f_parallel_in_v[i_x][i_v];
        chosen_adv->v = tmp*v[i_v+par_variables->global_indices[1]];
        adv1d_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
        for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++)
            par_variables->f_parallel_in_v[i_x][i_v] = par_variables->f_1d[i_x];
    chosen_adv->v = tmp;
    }
} else {
    adv1d_non_periodic_lag_t* chosen_adv = adv->non_periodic_adv;
    tmp = chosen_adv->v;
    for (i_v = 0; i_v < par_variables->size_v_par_v; i_v++){
        for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++)
            par_variables->f_1d[i_x] = par_variables->f_parallel_in_v[i_x][i_v];
        chosen_adv->v = tmp*v[i_v+par_variables->global_indices[1]];
        adv1d_non_periodic_lag_compute(chosen_adv, par_variables->f_1d, dt);
        for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++)
            par_variables->f_parallel_in_v[i_x][i_v] = par_variables->f_1d[i_x];
    }
    chosen_adv->v = tmp;
}
}


#endif // ifndef SELA_VP_1D1V_CART_ADVECT
