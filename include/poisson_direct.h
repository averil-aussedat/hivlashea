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
    double lambda;  // simulation parameter
    double nu;      // simulation parameter
};

/*
 * @param[in] mesh, the mesh on which we're working.
 * @return    a poisson solver that uses integration.
 */
poisson_solver_direct new_poisson_solver_direct(double *x, int sizex,
        double lambda, double nu) {
    poisson_solver_direct p;
    p.x = x;
    p.sizex = sizex;
    p.lambda = lambda;
    p.nu = nu;
    return p;
}

void update_E_from_rho_and_current_1d_OLD(poisson_solver_direct p,
        double delta_t,
        double* rho, double* current, double* E) {
    // There is one less interval than the number of points
    double delta_x = (p.x[p.sizex - 1] - p.x[0]) / (p.sizex - 1);
    
    // The domain is [-1..1]
    double E_at_one = E[p.sizex - 1];
    double E_at_minus_one = E[0];
    double mass = compute_mass(p.x, p.sizex, rho);
    
    //printf("E_at_one=%1.20g\n",E_at_one);
    
    // Formula: lambda^2 * d_t E(t,1) + J(t,1) = nu/2 mass(t)
    double E_at_one_new = (p.nu / 2. * mass - current[p.sizex-1]) * delta_t / (p.lambda * p.lambda) + E_at_one;
    // Formula: lambda^2 * d_t E(t,-1) + J(t,1) = -nu/2 mass(t)
    double E_at_minus_one_new = (-p.nu / 2. * mass - current[p.sizex-1]) * delta_t / (p.lambda * p.lambda) + E_at_minus_one;
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

/*
 * Copied / pasted from Mehdi BADSI.
 */
void update_E_from_rho_and_current_1d(poisson_solver_direct p,
        double dt,
        double Mass_e, double* rho, double* current, double* E) {
    double E_xmax, E_xmin;
    int Nx = p.sizex - 1;
    double dx = (p.x[Nx] - p.x[0]) / Nx;
    double nu = p.nu;
    double lambda = p.lambda;
    double J_right = current[Nx];
    double J_left = current[0];
    /*
    printf("J_right = %f; J_left = %f.\n", J_right, J_left);
    printf("nu = %f, lambda = %f\n", nu, lambda);
    for (int index = 0; index <= Nx; index++) {
        printf("E[%d] = %f ", index, E[index]);
    }
    printf("\n");
    for (int index = 0; index <= Nx; index++) {
        printf("rho[%d] = %f ", index, rho[index]);
    }
    printf("\n");
    */
    
    // Compute the new boundary conditions on the electric_field;
    E_xmax = E[Nx] + 0.5 * (nu * dt/(lambda*lambda)) * Mass_e - (dt/(lambda*lambda)) * J_right;
    E_xmin = E[0]  - 0.5 * (nu * dt/(lambda*lambda)) * Mass_e - (dt/(lambda*lambda)) * J_left;
    
    double mass_rho_minus; // Quadrature for int_{-1}^{x_i} rho(x)dx.
    double mass_rho_plus;  // Quadrature for int_{x_i}^{1} rho(x)dx
    for (int i = 0; i < Nx+1; i++) {
        mass_rho_plus  = 0;
        mass_rho_minus = 0;
        for(int k = 0; k <= i-1; k++) {
            mass_rho_minus += 0.5 * dx * (rho[k+1] + rho[k]);
        }
        for(int k = i; k <= Nx-1; k++) {
            mass_rho_plus += 0.5 * dx * (rho[k+1] + rho[k]);
        }
        E[i] = (1./(2.*lambda*lambda)) * (mass_rho_minus - mass_rho_plus) + 0.5 * (E_xmax + E_xmin);
    }
    /*
    for (int index = 0; index <= Nx; index++) {
        printf("E[%d] = %f ", index, E[index]);
    }
    printf("\n");
    exit(1);
    */
}

#endif // ifndef SELA_VP_1D1V_CART_POISSON_DIRECT
