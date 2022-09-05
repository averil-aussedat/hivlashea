// Bank of unitary tests.
//

// todo : linear advection + poisson solver


#include <math.h>                   // function  floor
#include <stdbool.h>                // type      bool
#include <stdio.h>                  // functions printf, fprintf
#include <stdlib.h>                 // type      size_t
// #include "adv1d_non_periodic_lag.h" // function no_boundary_lag_compute_index_and_alpha
#include "advect.h"                 // types     adv1d_x_t, adv1d_v_t, adv1d_[non_]periodic_lag_t
                                    // functions adv1d_x_init, adv1d_y_init, advection_x, advection_v, init_extrap
                                    // functions adv1d_compute_lag, adv1d_compute_i0_and_alpha
#include "diag.h"                   // functions diag_energy, diag_f, diag_1d
#include "mesh_1d.h"                // type      mesh_1d
                                    // function  mesh_1d_create
#include "poisson_direct.h"         // type      poisson_solver_direct
                                    // function  new_poisson_solver_direct
//#include "poisson_periodic_fft.h" // type      poisson_solver_periodic_fft
                                    // function  new_poisson_solver_periodic_fft
                                    // COMMONFUN update_E_from_rho_and_current_1d
#include "rho.h"                    // functions update_spatial_density, update_current, compute_mass
                                    // update_rho, update_rho_and_current
#include "string_helpers.h"         // macro     ERROR_MESSAGE

// MACROS
#define VERBOSE True // comment to hide terminal output

int advection1D_test () {
    int N = 99; // number of intervals/cells
    int i=0; // index
    adv1d_non_periodic_lag_t* adv = (adv1d_non_periodic_lag_t*)malloc(sizeof(adv1d_non_periodic_lag_t));

    adv->d = 5;
    adv->kb = 3;
    adv->min = -1.0;
    adv->max = 3.0;
    adv->N = N; // number of points - 1
    adv->v = -2.0;
    adv->lag = malloc((2*adv->d+2)*sizeof(double));
    adv->buf = malloc((N+1+2*adv->d)*sizeof(double));
    adv->extrap = malloc((adv->kb+1)*sizeof(double));
    init_extrap (adv->kb, adv->extrap);

    printf("extrap : ");
    for (i=0; i<=adv->kb; ++i) {
        printf("%f ", adv->extrap[i]);
    }
    printf("\n");

    double dx = (adv->max - adv->min) / N; 
    double* xx = (double*) malloc((N+1)*sizeof(double));
    double* ff = (double*) malloc((N+1)*sizeof(double));
    for (i=0; i<=N; ++i) { // do not forget to change error computation
        xx[i] = adv->min + i * dx;
        // ff[i] = cos(xx[i])-cos(-1.0);
        ff[i] = cos(xx[i])-cos(3.0);
        // ff[i] = exp(-xx[i]*xx[i]/0.1);
        // ff[i] = 1.0;
        // ff[i] = i * dx; 
    }

    // (double *E, double *x, int sizex, char* array_name, int iplot) 
    diag_1d (ff, xx, N+1, "test_adv", "", 0, 0.0);

    double dt = 0.01;
    int Nt=50;
    double T = Nt * dt;
    for (i=1; i<=Nt; ++i) {
        adv1d_non_periodic_lag_compute(adv, ff, dt);
    }
    printf("T = %f, nu : %f\n", T, adv->v * dt / dx);

    double error_max = -1.0, sol=0.0, yy=0.0;
    for (i=0; i<=N; ++i) { // do not forget to change initial condition
        yy = xx[i] - adv->v * T;
        sol = (yy>3.0 ? 0.0 :cos(yy)-cos(3.0)); // truncate to have 0 boundary condition (! negative speed only)
        // error_max = fmax(error_max, fabs(ff[i] - cos(yy)+cos(-1.0)));
        // error_max = fmax(error_max, fabs(ff[i] - exp(-yy*yy/0.1)));
        // error_max = fmax(error_max, fabs(ff[i] - 1.0));
        // error_max = fmax(error_max, fabs(ff[i] - i*dx));
        error_max = fmax(error_max, fabs(ff[i] - sol));
    }
    printf("Erreur maximale : %f.\n", error_max);
    diag_1d (ff, xx, N+1, "test_adv", "", 1, 0.0);

    free (xx); 
    free (ff); 
    free (adv->lag);
    free (adv->buf);
    free (adv->extrap);
    free (adv);

    return 0;
}

/*
 * Inputs : 
 *     * number of mesh points (intervals+1)
 * Output : 
 *     * error at final time
 *
 */
double solverE_test (int sizex) {
    int i=0;
    double lambda = 1.0, nu=1.0, dx = (1.0 - (-1.0))/(sizex-1);
    double * x = (double*)malloc(sizex*sizeof(double));

    for (i=0; i<sizex; ++i) { // uniform mesh from -1. to 1. with sizex points, i.e sizex-1 cells 
        x[i] = -1.0 + dx * i;
    }
    poisson_solver_direct psolv = new_poisson_solver_direct (x, sizex, lambda, nu);

    double dt=1.0, Mass_e=0.0, error=0.0;
    double* current = (double*)malloc(sizex*sizeof(double)); // J 
    double* rho = (double*)malloc(sizex*sizeof(double)); // what we integrate
    double* E = (double*)malloc(sizex*sizeof(double)); // in and out variable
    double* realE = (double*)malloc(sizex*sizeof(double)); // exact solution

    int thecase = 0; // sinus, no update of E at boundaries
    // initialization
    switch (thecase) {
        case 0: {
            Mass_e = 0.0;
            for (i=0; i<sizex; ++i) {
                current[i] = 0.0; // conjugated with Mass_e=0, yields no update of E in +-1            
                rho[i] = cos(x[i]); // pair function
                E[i] = 0.0; // initially null !!!! set boundary conditions accordingly
                realE[i] = sin(x[i]); // analytical solution
            }
            E[sizex-1] = sin(1.0); // int_0^{1} rho(y) dy
            E[0] = - E[sizex-1];
            break;
        }
        default: {
            printf("The case %d unknown. All to 0 by default\n", thecase);            
            for (i=0; i<sizex; ++i) {
                current[i] = 0.0; // conjugated with Mass_e=0, yields no update of E in +-1            
                rho[i] = 0.0; // pair function
                E[i] = 0.0; // initially null !!!! set boundary conditions accordingly
            }
        }
    }

    diag_1d (E, x, sizex, "test_E_E", "", 0, 0.0);
    diag_1d (rho, x, sizex, "test_E_rho", "", 0, 0.0);

    update_E_from_rho_and_current_1d (psolv, dt, Mass_e, rho, current, E); // modifies E

    error = -1.0;
    for (i=0; i<sizex; ++i) {
        error = fmax(error, fabs(E[i] - realE[i])); 
    }

    diag_1d (E, x, sizex, "test_E_E", "", 1, 0.0);

    free (x);
    free (current);
    free (rho);
    free (E);
    free (realE);
    return error;
}


int main(int argc, char *argv[]) {
    #ifdef VERBOSE
        printf("Welcome in unitary_tests.\n");
    #endif

    FILE* errfile = fopen("test_errors.dat", "w");
    for (int sizex=100; sizex<=1000; sizex+=100) {
        fprintf(errfile, "%f\t%.10f\n", 2.0/(sizex-1), solverE_test (sizex));
    }

    // advection1D_test ();  
    // solverE_test (100);  
    fclose (errfile);

    #ifdef VERBOSE
        printf("Bye\n");
    #endif
    return 0;
}

