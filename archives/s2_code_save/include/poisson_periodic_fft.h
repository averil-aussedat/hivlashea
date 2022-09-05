#ifndef SELA_VP_1D1V_CART_POISSON
#define SELA_VP_1D1V_CART_POISSON

/*****************************************************************************
 *                             Poisson solver                                *
 *                                                                           *
 * selalib/src/field_solvers/poisson_solvers/sll_m_poisson_2d_periodic.F90   *
 *****************************************************************************/

#include <complex.h> // see http://www.fftw.org/fftw3_doc/Complex-numbers.html
#include <fftw3.h>   // type      fftw_plan
                     // functions fftw_plan_r2r_1d, fftw_execute
                     // constants FFTW_R2HC, FFTW_ESTIMATE, FFTW_HC2R
#include <math.h>    // constant  M_PI
#include <stdlib.h>  // function  malloc

typedef struct poisson_solver_periodic_fft poisson_solver_periodic_fft;
struct poisson_solver_periodic_fft {
    int nc;        // cells number
    double L;
    double *in;
    double *out;
    fftw_plan p1;
    fftw_plan p2;
};

/*
 * @param[in] mesh, the mesh on which we're working.
 * @return    a poisson solver that uses the library FFTW.
 */
poisson_solver_periodic_fft new_poisson_solver_periodic_fft(double *x, int sizex) {
    poisson_solver_periodic_fft p;
    p.nc = sizex - 1;
    p.L = x[sizex - 1] - x[0];
    p.in  = malloc(p.nc * sizeof(double));
    p.out = malloc(p.nc * sizeof(double));
    p.p1 = fftw_plan_r2r_1d(p.nc, p.in, p.out, FFTW_R2HC, FFTW_ESTIMATE);
    p.p2 = fftw_plan_r2r_1d(p.nc, p.in, p.out, FFTW_HC2R, FFTW_ESTIMATE);
    return p;
}

void poisson1drealftw(double *in, double *out, fftw_plan p1, fftw_plan p2, double *E, double L, int N){
    // L = xmax - xmin
    // we suppose that E = rho - 1 at the beginning
    // corresponds to epsilon = -1 in the pdf file
    int i;
    double tmp = 0.5 * (L / M_PI) / (double)N;  
    double re, im;
    
    //realftw
    for (i = 0; i < N; i++)
        in[i] = E[i];
    fftw_execute(p1);
    for (i = 0; i < N; i++)
        E[i] = out[i];
    
    //Re(p_0),...,Re(p_{N/2}),Im(p_{(N+1)/2-1}),...,Im(p_{1})  
    E[0] = 0.;
    if (N%2 == 0)
        E[N/2] = 0.;
    for (i = 1; i < (N+1)/2; i++) {
        re = E[i];
        im = E[N-i];
        E[i] = (tmp/(double)i)*im;
        E[N-i] = -(tmp/(double)i)*re;
    }
    
    //realftw inverse    
    for (i = 0; i < N; i++)
        in[i] = E[i];
    fftw_execute(p2);
    for (i = 0; i < N; i++)
        E[i] = out[i];
}

void update_E_from_rho_and_current_1d(poisson_solver_periodic_fft p,
        double delta_t, double Mass_e, double* rho, double* current, double* E) {
    for (int j = 0; j < p.nc; j++)
        E[j] = rho[j] - 1.;
    poisson1drealftw(p.in, p.out, p.p1, p.p2, E, p.L, p.nc);
    E[p.nc] = E[0];
}

#endif // ifndef SELA_VP_1D1V_CART_POISSON
