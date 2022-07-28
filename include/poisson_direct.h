#ifndef SELA_VP_1D1V_CART_POISSON_DIRECT
#define SELA_VP_1D1V_CART_POISSON_DIRECT

#include "rho.h"              // functions compute_mass, update_current
#include "string_helpers.h"   // macro     ERROR_MESSAGE

/***********************************************************************************
 *                Poisson solver --- Only for CEMRACS 2022 project                 *
 * E(x)=1/(2*lambda^2)*(int_{-1}^x rho(x) dx - int_x^1 rho(x) dx) + (E(-1)+E(1))/2 *
 ***********************************************************************************/

typedef struct poisson_solver_direct poisson_solver_direct;
struct poisson_solver_direct {
    double* x;      // mesh
    int sizex;      // cells number
    double delta_t; // simulation parameter
    double lambda;  // simulation parameter
    double nu;      // simulation parameter
};

/*
 * @param[in] mesh, the mesh on which we're working.
 * @return    a poisson solver that uses integration.
 */
poisson_solver_direct new_poisson_solver_direct(double *x, int sizex,
        double delta_t, double lambda, double nu) {
    poisson_solver_direct p;
    p.x = x;
    p.sizex = sizex;
    p.delta_t = delta_t;
    p.lambda = lambda;
    p.nu = nu;
    return p;
}

void compute_E_from_rho_1d(poisson_solver_direct p, double* rho, double* current, double* E) {
    // There is one less interval than the number of points
    double delta_x = (p.x[p.sizex - 1] - p.x[0]) / (p.sizex - 1);
    
    // The domain is [-1..1]
    double E_at_one = E[p.sizex - 1];
    double E_at_minus_one = E[0];
    double mass = compute_mass(p.x, p.sizex, rho);
    
    //printf("E_at_one=%1.20g\n",E_at_one);
    
    
    
    // Formula: lambda^2 * d_t E(t,1) + J(t,1) = nu/2 mass(t)
    double E_at_one_new = (p.nu / 2. * mass - current[p.sizex-1]) * p.delta_t / (p.lambda * p.lambda) + E_at_one;
    // Formula: lambda^2 * d_t E(t,-1) + J(t,1) = -nu/2 mass(t)
    double E_at_minus_one_new = (-p.nu / 2. * mass - current[p.sizex-1]) * p.delta_t / (p.lambda * p.lambda) + E_at_minus_one;
    double should_be_zero = (E_at_minus_one_new + E_at_one_new) / 2.;
    
    for (int i = 0; i < p.sizex; i++) {
        int j;
        // Sub-formula: int_{-1}^x rho(x) dx
        // Here, uses left rectangles.
        double first_integral = 0.;
        for (j = 0; j < i; j++) {
            first_integral += rho[j];
        }
        first_integral *= delta_x;
        // Sub-formula: int_x^1 rho(x) dx
        // Here, uses left rectangles.
        double second_integral = 0.;
        for (j = i; j < p.sizex; j++) {
            second_integral += rho[j];
        }
        second_integral *= delta_x;
        // Formula: E(x)=1/(2*lambda^2)*(int_{-1}^x rho(x) dx - int_x^1 rho(x) dx) + (E(-1)+E(1))/2
        E[i] = 1. / (2. * p.lambda * p.lambda) * (first_integral - second_integral) + should_be_zero;
    }
}

#endif // ifndef SELA_VP_1D1V_CART_POISSON_DIRECT
