#ifndef SELA_VP_1D1V_CART_ADVECT
#define SELA_VP_1D1V_CART_ADVECT

#include "adv1d_non_periodic_lag.h" // type      adv1d_t
                                    // functions adv1d_compute, adv1d_init
                                    // (used not directly in this file but
                                    //  used by simulations that include this file)
#include "remap.h"                  // type      parallel_stuff
                                    // functions exchange_parallelizations, local_to_global_2d

// If you wish to change to another advection, you can include another advection file and
// change the following lines:
#define adv1d_t adv1d_non_periodic_lag_t
#define adv1d_compute adv1d_non_periodic_lag_compute
#define adv1d_init adv1d_non_periodic_lag_init

void advection_v(parallel_stuff* par_variables, adv1d_t* adv, double dt, double* E){
    int i_x;
    int i_v;
    double tmp;
    tmp = adv->v;
    if (!par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    local_to_global_2d(par_variables, 0, 0);
    for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++)
            par_variables->f_1d_in[i_v] = par_variables->f_parallel_in_x[i_x][i_v];
        adv->v = tmp*E[i_x+par_variables->global_indices[0]];
        adv1d_compute(adv, par_variables->f_1d_in, par_variables->f_1d_out, dt);
        //interpol_1d_semi_lag_advect(&par_variables->lag_v, par_variables->f_1d);
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++)
            par_variables->f_parallel_in_x[i_x][i_v] = par_variables->f_1d_out[i_v];
    }
    adv->v = tmp;
}

void advection_x(parallel_stuff* par_variables, adv1d_t* adv, double dt, double* v){
    int i_x;
    int i_v;
    double tmp;
    tmp = adv->v;
    if (par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    local_to_global_2d(par_variables, 0, 0);
    for (i_v = 0; i_v < par_variables->size_v_par_v; i_v++){
        for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++)
            par_variables->f_1d_in[i_x] = par_variables->f_parallel_in_v[i_x][i_v];
        adv->v = tmp*v[i_v+par_variables->global_indices[1]];
        adv1d_compute(adv, par_variables->f_1d_in, par_variables->f_1d_out, dt);
        //printf("par_variables->size_x_par_v=%d\n",par_variables->size_x_par_v);
        //interpol_1d_semi_lag_advect(&par_variables->lag_x, par_variables->f_1d);
        for (i_x = 0; i_x < par_variables->size_x_par_v; i_x++)
            par_variables->f_parallel_in_v[i_x][i_v] = par_variables->f_1d_out[i_x];
    }
    adv->v = tmp;
}





#endif // ifndef SELA_VP_1D1V_CART_ADVECT
