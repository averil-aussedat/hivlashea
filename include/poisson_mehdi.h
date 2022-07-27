#ifndef SELA_VP_1D1V_CART_POISSON_MEHDI
#define SELA_VP_1D1V_CART_POISSON_MEHDI

#include <stdlib.h>         // functions malloc, free ((de)allocate memory)
#include "array_helpers.h"  // functions allocate_2d_array, deallocate_2d_array
#include "rho.h"            // functions compute_mass, update_current
#include "string_helpers.h" // macro     ERROR_MESSAGE

/*****************************************************************************
 *                             Poisson solver                                *
 * Author: Mehdi Badsi (transcription from C++: Yann Barsamian)              *
 *****************************************************************************/

typedef struct poisson_solver_mehdi poisson_solver_mehdi;
struct poisson_solver_mehdi {
    double* x;       // mesh
    int sizex;       // cells number
    // The full domain is on [0 .. sizex - 1] which is, with other notations,
    // [-Nx .. Nx] (thus, p.size is odd: sizex = Nx + 1).
    // As a consequence, to access array[i_bis] with i_bis in [-Nx .. Nx], we
    // look for array[i_bis - Nx + sizex - 1].
    int Nx;
    double** Lap;    // Laplacian matrix
    double* sigma;   // permutation
    double* phi;     // 
    double* rhs;     // right-hand side for the A * X = B equation
    double delta_t;  // simulation parameter
    double lambda;   // simulation parameter
    double nu;       // simulation parameter
};

/*
 * @param[in] mesh, the mesh on which we're working.
 * @return    a poisson solver that uses the Gauss elimination.
 */
poisson_solver_mehdi new_poisson_solver_mehdi(double *x, int sizex,
        double delta_t, double lambda, double nu) {
    poisson_solver_mehdi p;
    p.x = x;
    p.sizex = sizex;
    p.Nx = (sizex - 1) / 2;
    p.Lap = allocate_2d_array(p.Nx, p.Nx);
    p.sigma = malloc(p.Nx * sizeof(double));
    p.phi = malloc(p.Nx * sizeof(double));
    p.rhs = malloc(p.Nx * sizeof(double));
    p.delta_t = delta_t;
    p.lambda = lambda;
    p.nu = nu;
    int i;
    
    // Define the matrix of the Laplacian with s.o consistent Neuman b.c
    for (i = 0; i < p.Nx; i++) {
        p.Lap[i][i] = -2.0;
    }
    for (i = 0 ; i < p.Nx - 1; i++) {
        p.Lap[i][i+1] = 1.0;
        p.Lap[i+1][i] = 1.0;
    }
    p.Lap[Nx-1][Nx-2] = 2.0;
    
    // LU Factorization of the Laplacian matrix
    for (i = 0; i < p.Nx; i++)
        p.sigma[i] = i;
    Gauss(p.Lap, p.sigma);
    return p;
}

#define EPSILON 1e-12

/*
 * Gauss elimination of matrix A (stores the permutation in sigma).
 */
void Gauss(double** A, double* sigma, int size) {
    double max, C;
    int stock, line;
    int i;
    for (int j = 0; j < size; j++) {
        // Test sur le pivot
        if (fabs(A[sigma[j]][j]) < EPSILON) {
            // Ligne du pivot nul
            line = sigma[j];
            // Initialisation du max
            max = A[sigma[j]][j];
            for (i = j; i < size; i++) {
                if (fabs(max) < fabs(A[sigma[i]][j])) {
                    max = fabs(A[sigma[i]][j]);
                    line = sigma[i];
                }
            }
            if (fabs(max) < EPSILON) {
                ERROR_MESSAGE("Bad matrix for the Gauss elimination: singular matrix.\n");
            }
            stock = sigma[j];
            sigma[j] = sigma[line];
            sigma[line] = stock;
        }
        for (i = j+1; i < size; i++) {
            C = A[sigma[i]][j] / A[sigma[j]][j];
            A[sigma[i]][j] = C;
            for (int k = j+1; k < size; k++) {
                A[sigma[i]][k] = A[sigma[i]][k] - C * A[sigma[j]][k];
            }
        }
    }
}

/*
 * Solves A * X = B.
 */
void Solve(double** A, double* sigma, double* B, int size, double* X) {
    int i, k;
    for (i = 0; i < size i++) {
        X[i] = B[sigma[i]];
        for (k = 0; k < i; k++) {
            X[i] = X[i] - A[sigma[i]][k] * X[k];
        }
    }
    for (i = size-1; i >= 0; i--) {
        for (k = i+1; k < size; k++) {
            X[i] = X[i] - A[sigma[i]][k] * X[k];
        }
        X[i] = X[i] / A[sigma[i]][i];
    }
}

void compute_E_from_rho_1d(poisson_solver p, double* rho, double* current, double* E) {
    double delta_x = (x[p.sizex - 1] - x[0]) / (p.sizex - 1);
    double mass = compute_mass(p.x, p.sizex, rho);
    double b_Nx = p.nu / 2 * mass - current[p.sizex];
    double f_neum = (delta_t / (p.lambda * p.lambda)) * b_Nx + E[p.sizex];
    
    // Define the right hand side with s.o consistent Neuman b.c
    for (int i = 1; i <= p.Nx; i++) {
        rhs[i - 1] = -(delta_x / p.lambda) * (delta_x / p.lambda) * rho[i - p.Nx + p.sizex - 1];
    }
    p.rhs[p.Nx - 1] += 2 * delta_x * f_neum;
    Solve(p.Lap, p.sigma, p.rhs, p.sizex, p.phi);
    // TODO: here, we have phi on [1 .. Nx]. We have to convert it to phi[-Nx .. Nx].
    // TODO: Convert phi to E and remove the error message.
    ERROR_MESSAGE("Don't use this Poisson solver, it is not finished.\n");
}

#endif // ifndef SELA_VP_1D1V_CART_POISSON_MEHDI
