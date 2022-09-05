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
    int sizex;      // mesh points number (so nummber of cells/intervals + 1)
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

void update_E_from_rho_and_current_1d_OLDOLD(poisson_solver_direct p,
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
void update_E_from_rho_and_current_1dOLD (poisson_solver_direct p,
        double dt,
        double Mass_e, double* rho, double* current, double* E) {
    double E_xmin, E_xmax; // prescribed values of E at xmin and xmax
    int Nx = p.sizex - 1; // number of intervals. Index from [0 to Nx]
    int i; // running index
    double dx = (p.x[Nx] - p.x[0]) / Nx; 
    double nu = p.nu;
    double lambda = p.lambda;
    double J_left  = current[0];
    double J_right = current[Nx];
    double mass_rho_minus=0.0; // Quadrature for  int_{-1}^{x_i} rho(x)dx.
    double mass_rho_plus =0.0; // Quadrature for -int_{x_i}^{ 1} rho(x)dx
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
    // E_xmax = E[Nx] + 0.5 * (nu * dt/(lambda*lambda)) * Mass_e - (dt/(lambda*lambda)) * J_right;
    // E_xmin = E[0]  - 0.5 * (nu * dt/(lambda*lambda)) * Mass_e - (dt/(lambda*lambda)) * J_left;
    E_xmin = E[0]  + dt/(lambda*lambda) * (-0.5 * nu * Mass_e - J_left);
    E_xmax = E[Nx] + dt/(lambda*lambda) * ( 0.5 * nu * Mass_e - J_right);
    for (i = 0; i < Nx+1; i++) {
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


/*
 * Double-integration "solver" with trapezes formula. Solves
 *  E(x) = 1/(2*lambda^2) [\int_{-1}^x rho(y)dy - \int_{x}^1 rho(y)dy] + 1/2(E[-1]+E[1])
 *
 *  NO PROPER BOUNDARY CONDITIONS ?
 *
 * @param[in] p : poisson solver data
 * @param[in] dt : time step
 * @param[in] Mass_e : triple integral of fe - fi on x,v and t\in[0,tn]
 * @param[out] rho :     int_{v}   (fe(t,x,v) - fi(t,x,v)) dv   for each xi
 * @param[out] current : int_{v} v (fe(t,x,v) - fi(t,x,v)) dv   for each xi
 * @param[out] E : electric field for each xi
 */
void update_E_from_rho_and_current_1d_nbc (poisson_solver_direct p, double dt,
        double Mass_e, double* rho, double* current, double* E) {

    // printf("Enters Poisson solver\n");
    double E_xmin, E_xmax; // prescribed values of E at xmin and xmax
    int Nx = p.sizex - 1; // number of intervals. Index from [0 to Nx]
    int i; // running index
    double dx = (p.x[Nx] - p.x[0]) / (double)Nx; 
    double nu = p.nu;
    double lambda = p.lambda;
    double J_left  = current[0];
    double J_right = current[Nx];
    double factor = 0.5 * 0.5 * dx / (lambda * lambda); // 0.5 from the expression of E and from trapezes
    double mass_rho_minus=0.0; // Quadrature for  int_{-1}^{x_i} rho(x)dx.
    double mass_rho_plus =0.0; // Quadrature for -int_{x_i}^{ 1} rho(x)dx

    // Boundary conditions
    E_xmin = E[0]  + dt/(lambda*lambda) * (-0.5 * nu * Mass_e - J_left);
    E_xmax = E[Nx] + dt/(lambda*lambda) * ( 0.5 * nu * Mass_e - J_right);
    E[0] = E_xmin; E[Nx] = E_xmax;  
    // printf("[update_E_from_rho_and_current_1d] |E(1) + E(-1)| : %6.9f\n", fabs(E_xmin+E_xmax));
    for (i=1; i<Nx; ++i) { E[i] = 0.0; } // reset the electric field
    // Linear-time integration with trapezes method (order 2)
    for (i=1; i<Nx; ++i) { // interior points only
        mass_rho_minus += factor * (rho[i-1] + rho[i]);       // forward  integration
        mass_rho_plus  -= factor * (rho[Nx-i+1] + rho[Nx-i]); // backward integration
        E[i]    += mass_rho_minus + E_xmin * 0.5;
        E[Nx-i] += mass_rho_plus  + E_xmax * 0.5;
    }
    // printf("Leaves Poisson solver\n");
}

/*
 * "solver" with trapezes formula. Solves
 *  E(x) = E(-1) - 1/(lambda*lambda) \int_{x}^{ 1} rho(y) dy
 *  E(x) = E( 1) + 1/(lambda*lambda) \int_{-1}^{x} rho(y) dy
 *  and takes a convex combination of both.
 *
 * @param[in] p : poisson solver data
 * @param[in] dt : time step
 * @param[in] Mass_e : triple integral of fe - fi on x,v and t\in[0,tn]
 * @param[out] rho :     int_{v}   (fe(t,x,v) - fi(t,x,v)) dv   for each xi
 * @param[out] current : int_{v} v (fe(t,x,v) - fi(t,x,v)) dv   for each xi
 * @param[out] E : electric field for each xi
 */
void update_E_from_rho_and_current_1d_cb (poisson_solver_direct p, double dt,
        double Mass_e, double* rho, double* current, double* E) {

    // printf("Enters Poisson solver\n");
    double E_xmin=0.0, E_xmax=0.0; // prescribed values of E at xmin and xmax
    int Nx = p.sizex - 1; // number of intervals. Index from [0 to Nx]
    int i; // running index
    double dx = (p.x[Nx] - p.x[0]) / (double)Nx; 
    // printf("dx : %f\n", dx);
    double nu = p.nu;
    double lambda = p.lambda;
    double J_left  = current[0];
    double J_right = current[Nx];
    double factor = 0.5 * dx / (lambda * lambda); // 0.5 from the expression of E and from trapezes
    double convex_factor = 0.0; 
    double mass_rho_minus=0.0; // Quadrature for   1/lambda^2 int_{-1}^{x_i} rho(x)dx.
    double mass_rho_plus =0.0; // Quadrature for - 1/lambda^2 int_{x_i}^{ 1} rho(x)dx

    double* Eforward = (double*)malloc (p.sizex*sizeof(double));
    double* Ebackward = (double*)malloc (p.sizex*sizeof(double));

    // Boundary conditions
    E_xmin = E[0]  + dt/(lambda*lambda) * (-0.5 * nu * Mass_e - J_left);
    E_xmax = E[Nx] + dt/(lambda*lambda) * ( 0.5 * nu * Mass_e - J_right);

    E_xmin = 0.0;
    E_xmax = 0.0;

    E[0] = E_xmin; E[Nx] = E_xmax;  

    // printf("[update_E_from_rho_and_current_1d] |E(1) + E(-1)| : %6.9f\n", fabs(E_xmin+E_xmax));
    for (i=1; i<Nx; ++i) { E[i] = 0.0; } // reset the electric field
    // Linear-time integration with trapezes method (order 2)
    for (i=1; i<Nx; ++i) { // interior points only
        mass_rho_minus += factor * (rho[i-1] + rho[i]);       // forward  integration
        mass_rho_plus  -= factor * (rho[Nx-i+1] + rho[Nx-i]); // backward integration
        convex_factor = 1.0 - (double)i/(double)Nx; // both boundary conditions
        // convex_factor = 0.0; // left-sided
        // convex_factor = 0.5; // middle
        // convex_factor = 1.0; // right-sided
        E[i]    += (mass_rho_minus + E_xmin) * convex_factor;
        E[Nx-i] += (mass_rho_plus  + E_xmax) * convex_factor;
        // printf("%d and %d : %f\n", i, Nx-i, convex_factor);
    }

    Eforward[0] = E_xmin;
    for (int ii=1; ii<p.sizex; ++ii) {
        Eforward[ii] = Eforward[ii-1] + factor * (rho[ii-1] + rho[ii]);
    }

    Ebackward[p.sizex-1] = E_xmax;
    for (int ii=p.sizex-2; ii>=0; --ii) {
        Ebackward[ii] = Ebackward[ii+1] - factor * (rho[ii+1] + rho[ii]);
    }

    double errf=0.0, errb=0.0, errc=0.0, mmu = 0.0;
    for (int ii=0; ii<p.sizex; ++ii) {
        errf = fmax(errf, fabs(Eforward[ii] - E[ii]));
        errb = fmax(errb, fabs(Ebackward[ii] - E[ii]));
        // mmu = 1.0 - (double)ii/(double)(p.sizex-1);
        mmu = (double)ii/(double)(p.sizex-1);
        errc = fmax(errc, fabs(((1.0 - mmu) * Eforward[ii] + mmu * Ebackward[ii]) - E[ii]));

        E[ii] = (1.0 - mmu) * Eforward[ii] + mmu * Ebackward[ii];
    }
    // printf("errf : %e, errb : %e, errc : %e\n", errf, errb);

    // double err_dxE = 0.0;
    // for (i=1; i<Nx; ++i) {
    //     printf("Ei : %e, rho_i : %e, err : %e \n", E[i], rho[i], lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
    //     err_dxE = fmax(err_dxE, lambda*lambda*(E[i+1] - E[i-1])/(2*dx) - rho[i]);
    // }
    // printf("Err_dxE : %e\n", err_dxE);

    // printf("Mass e : %e, Jleft : %e, Jright : %e, evolution : (%e,%e), Exmin/max : (%e, %e).\n", Mass_e, J_left, J_right, 
    //     dt/(lambda*lambda) * (-0.5 * nu * Mass_e - J_left), dt/(lambda*lambda) * ( 0.5 * nu * Mass_e - J_right), E_xmin, E_xmax);

    // // juste pour tester
    // printf("Attention test dans poisson direct\n");
    // E[0] = mass_rho_plus; E[Nx] = mass_rho_minus;

    // printf("Leaves Poisson solver\n");

    free(Eforward);
    free(Ebackward);
}

/*
 * "solver" with trapezes formula. Solves
 *  E(x) = E(-1) - 1/(lambda*lambda) \int_{x}^{ 1} rho(y) dy
 *  E(x) = E( 1) + 1/(lambda*lambda) \int_{-1}^{x} rho(y) dy
 *  and takes a convex combination of both.
 *
 * @param[in] N : number of intervals of the space mesh
 * @param[in] dx : space step
 * @param[in] lambda : Debye length, real
 * @param[out] rho :  int_{v}   (fe(t,x,v) - fi(t,x,v)) dv   for each xi
 * @param[out] E : electric field for each xi
 */
void update_E_from_rho_and_current_1d (int N, double dx, double lambda, double* rho, double* E) {

    int i=0, i0 = (int)(N/2); // middle of the domain 
    double factor = 0.5 * dx / (lambda*lambda); // 0.5 from trapezes

    if (2*i0 != N) {
        printf("\n\n[update_E_from_rho_and_current_1d] WARNING : 0 n'appartient pas au maillage. Milieu i0 = %d, N=%d.\n\n", i0, N);
    }

    // Boundary conditions
    E[i0] = 0.0; // not even needed since it is never updated and E[0]=0 at t=0, but more safe
    for (i=i0+1; i<=N; ++i) {
        E[i] = E[i-1] + factor * (rho[i-1] + rho[i]);
    }
    for (i=i0-1; i>=0; --i) {
        E[i] = E[i+1] - factor * (rho[i+1] + rho[i]);
    }
}

#endif // ifndef SELA_VP_1D1V_CART_POISSON_DIRECT
