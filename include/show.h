#ifndef SELA_SHOW
#define SELA_SHOW

#include <stdio.h>                  // functions printf, fprintf
#include "mesh_1d.h"                // type      mesh_1d
#include "advect.h"                 // types     adv1d_x_t, adv1d_v_t

/***********************************************************************************
 *                Bank of printing functions --- CEMRACS 2022 project              *
 ***********************************************************************************/

void show_mesh1d (mesh_1d *mesh, const char* name) {
    printf("%8s : %4d points in [%6.2f, %6.2f].\n", name, mesh->size,  mesh->min,  mesh->max);
} 

void show_adv (adv1d_non_periodic_lag_t *adv, const char* name) {
    int k=0;
    printf("%8s : d=%2d, N=%4d, kb=%2d, interval [%6.2f, %6.2f], v=%6.2f", 
        name, adv->d, adv->N, adv->kb, adv->min, adv->max, adv->v);

    // printf(", lag = [");
    // for (k=0; k<=2*adv->d+1; ++k) {
    //     printf("%f ", adv->lag[k]);
    // }
    // printf("]");
    printf(", extrap = [");
    for (k=0; k<=adv->kb; ++k) {
        printf("%6.2f ", adv->extrap[k]);
    }
    printf("]");
    printf(".\n");
}

 #endif // ifndef SELA_SHOW