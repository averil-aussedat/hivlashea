#ifndef SELA_VP_1D1V_CART_RHO
#define SELA_VP_1D1V_CART_RHO

#include <mpi.h>     // function  MPI_Allgatherv
                     // constants MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
#include <stdbool.h> // type      bool
#include <stdlib.h>  // type      size_t
#include "remap.h"   // type      parallel_stuff
                     // function  exchange_parallelizations

/*****************************************************************************
 *                                Update rho                                 *
 *****************************************************************************/

/*
 * Computes rho = int_{v=v_min}^{v_max} f(x,v) dv. If the current parallelization is on v,
 * starts by changing it to parallel on x.
 *
 * @param[in, out] par_variables
 * @param[in]      spatial_mesh (x[0..size_x-1], size_x) : mesh on x
 * @param[in]      velocity_mesh (v[0..size_v-1], size_v) : mesh on v
 * @param[out]     rho : charge density (integral of f on v)
 */
void update_spatial_density(parallel_stuff* par_variables, double *x, int sizex,
        double *v, int sizev, double* rho, bool is_periodic) {
    int pos;
    
    double delta_v = (v[sizev - 1] - v[0]) / (double)(sizev - 1);
    
    if (!par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    
    // Integration of f_local in v to generate linearised rho_local
    pos = 0;
    for (size_t i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        par_variables->send_buf[pos] = 0.;
        // // Integration: left rectangles method
        // for (size_t i_v = 0; i_v < par_variables->size_v_par_x-1; i_v++) {
        //     par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][i_v];
        // }
        // Integration: trapezes method
        par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][0]*0.5;
        for (size_t i_v = 1; i_v < par_variables->size_v_par_x-1; i_v++) {
            par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][i_v];
        }
        par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][par_variables->size_v_par_x-1]*0.5;
        par_variables->send_buf[pos] *= delta_v;
    	pos++;
    }
    // Gather the concatenation of all the linearised rho_locals
    MPI_Allgatherv(par_variables->send_buf, 
    	par_variables->size_x_par_x, MPI_DOUBLE_PRECISION,
        par_variables->recv_buf, par_variables->recv_counts, 
        par_variables->displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
    
    // Rebuild rho from the concatenation of the linearised rho_locals
    pos = 0;
    for (int process = 0; process < par_variables->mpi_world_size; process++) {
        for (size_t i = par_variables->layout_par_x.boxes[process].i_min; i <= par_variables->layout_par_x.boxes[process].i_max; i++)
            rho[i] = par_variables->recv_buf[pos++];
//        printf("Last i value is: %d.\n", par_variables->layout_par_x.boxes[process].i_max);
//        printf("sizex - 1 value is: %d.\n", sizex - 1);
    }
    
    // For periodic case only
    if (is_periodic) {
        rho[sizex - 1] = rho[0];
    } else {
        // Symmetrization of rho
//         for (int i=0; i<(int)(sizex-1)/2; ++i) { 
//             rho[i] = (rho[i] + rho[sizex-1-i])*0.5;
//             rho[sizex-1-i] = rho[i];
//         }
    }
}

/*
 * Computes mass = int_{x=x_min}^{x_max} rho(x) dx.
 *
 * @param[in] spatial_mesh (x[0..size_x-1], size_x) : mesh on x
 * @param[in] rho : charge density
 */
double compute_mass(double *x, int sizex, double* rho) {
    double delta_x = (x[sizex - 1] - x[0]) / (double)(sizex - 1);
    double mass = 0.0;
    // Integration: trapezius method
    mass += rho[0] / 2.;
    for (size_t i_x = 1; i_x < sizex - 1; i_x++) {
        mass += rho[i_x];
    }
    mass += rho[sizex - 1] / 2.;
    mass *= delta_x;
    return mass;
}

/*
 * Computes currentrho = int_{v=v_min}^{v_max} f(x,v) v dv. If the current parallelization is on v,
 * starts by changing it to parallel on x.
 *
 * @param[in, out] par_variables
 * @param[in]      spatial_mesh (x[0..size_x-1], size_x) : mesh on x
 * @param[in]      velocity_mesh (v[0..size_v-1], size_v) : mesh on v
 * @param[out]     current : integral of f*v on v
 */
void update_current(parallel_stuff* par_variables, double *x, int sizex,
        double *v, int sizev, double* current, bool is_periodic) {
    int pos;
    
    double delta_v = (v[sizev - 1] - v[0]) / (double)(sizev - 1);
    
    if (!par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    
    // Integration of f_local*v in v to generate linearised current_local
    pos = 0;
    for (size_t i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        par_variables->send_buf[pos] = 0.;
        // // Integration: left rectangles method
        // for (size_t i_v = 0; i_v < par_variables->size_v_par_x-1; i_v++)
        //     par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][i_v] * v[i_v];
        // par_variables->send_buf[pos] *= delta_v;
        // Integration: trapezes method
        par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][0]*0.5*v[0];
        for (size_t i_v = 1; i_v < par_variables->size_v_par_x-1; i_v++) {
            par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][i_v] * v[i_v];
        }
        par_variables->send_buf[pos] += par_variables->f_parallel_in_x[i_x][par_variables->size_v_par_x-1]*0.5*v[par_variables->size_v_par_x-1];
        par_variables->send_buf[pos] *= delta_v;
    	pos++;
    }
    // Gather the concatenation of all the linearised current_local
    MPI_Allgatherv(par_variables->send_buf, 
    	par_variables->size_x_par_x, MPI_DOUBLE_PRECISION,
        par_variables->recv_buf, par_variables->recv_counts, 
        par_variables->displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
    
    // Rebuild current from the concatenation of the linearised current_local
    pos = 0;
    for (int process = 0; process < par_variables->mpi_world_size; process++)
        for (size_t i = par_variables->layout_par_x.boxes[process].i_min; i <= par_variables->layout_par_x.boxes[process].i_max; i++)
            current[i] = par_variables->recv_buf[pos++];
    
    // For periodic case only
    if (is_periodic) {
        current[sizex - 1] = current[0];
    }
}

/*
 * Creates diagnostics of a function f.
 */
void compute_diag_f(parallel_stuff* par_variables, double *x, int sizex,
        double *v, int sizev, double* diag) {
    size_t i, i_x, i_v;
    int process, pos;
    double tmp;
    
    double delta_v = (v[sizev - 1] - v[0]) / (double)(sizev - 1);
    double delta_x = (x[sizex - 1] - x[0]) / (double)(sizex - 1);
    
    //printf("par_variables->size_v_par_x=%d\n",par_variables->size_v_par_x);
    if (!par_variables->is_par_x)
        exchange_parallelizations(par_variables);
    
    /**************************************************************************
     *                              Diag 0
     *************************************************************************/
    // Integration of [f_local XXX] to generate linearised [XXX_local]
    pos = 0;
    for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        par_variables->send_buf[pos] = 0.;
        //for (i_v = 0; i_v < par_variables->size_v_par_x-1; i_v++){
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++){
        	tmp = par_variables->f_parallel_in_x[i_x][i_v];
            par_variables->send_buf[pos] += tmp*tmp;
        }
        par_variables->send_buf[pos] *= delta_v;
    	pos++;
    }
    // Gather the concatenation of all the linearised [XXX_local]
    MPI_Allgatherv(par_variables->send_buf, 
    	par_variables->size_x_par_x, MPI_DOUBLE_PRECISION,
        par_variables->recv_buf, par_variables->recv_counts, 
        par_variables->displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
    
    // Rebuild [XXX] from the concatenation of the linearised [XXX_local]
    diag[0] = 0.;
    pos = 0;
    for (process = 0; process < par_variables->mpi_world_size; process++)
        for (i = par_variables->layout_par_x.boxes[process].i_min; i <= par_variables->layout_par_x.boxes[process].i_max; i++)
            diag[0] += par_variables->recv_buf[pos++];
    diag[0] *= delta_x;
    diag[0] = sqrt(diag[0]);
    
    /**************************************************************************
     *                              Diag 1
     *************************************************************************/
    // Integration of [f_local XXX] to generate linearised [XXX_local]
    pos = 0;
    for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        par_variables->send_buf[pos] = 0.;
        //for (i_v = 0; i_v < par_variables->size_v_par_x-1; i_v++){
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++){
        	tmp = par_variables->f_parallel_in_x[i_x][i_v];
            par_variables->send_buf[pos] += tmp*v[i_v]*v[i_v];
        }
        par_variables->send_buf[pos] *= delta_v;
    	pos++;
    }
    // Gather the concatenation of all the linearised [XXX_local]
    MPI_Allgatherv(par_variables->send_buf, 
    	par_variables->size_x_par_x, MPI_DOUBLE_PRECISION,
        par_variables->recv_buf, par_variables->recv_counts, 
        par_variables->displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
    
    // Rebuild [XXX] from the concatenation of the linearised [XXX_local]
    diag[1] = 0.;
    pos = 0;
    for (process = 0; process < par_variables->mpi_world_size; process++)
        for (i = par_variables->layout_par_x.boxes[process].i_min; i <= par_variables->layout_par_x.boxes[process].i_max; i++)
            diag[1] += par_variables->recv_buf[pos++];
    diag[1] *= delta_x;
    
    /**************************************************************************
     *                              Diag 2
     *************************************************************************/
    // Integration of [f_local XXX] to generate linearised [XXX_local]
    pos = 0;
    for (i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        par_variables->send_buf[pos] = 0.;
        //for (i_v = 0; i_v < par_variables->size_v_par_x-1; i_v++){
        for (i_v = 0; i_v < par_variables->size_v_par_x; i_v++){
        	tmp = par_variables->f_parallel_in_x[i_x][i_v];
            par_variables->send_buf[pos] += fabs(tmp);
        }
        par_variables->send_buf[pos] *= delta_v;
    	pos++;
    }
    // Gather the concatenation of all the linearised [XXX_local]
    MPI_Allgatherv(par_variables->send_buf, 
    	par_variables->size_x_par_x, MPI_DOUBLE_PRECISION,
        par_variables->recv_buf, par_variables->recv_counts, 
        par_variables->displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
    
    // Rebuild [XXX] from the concatenation of the linearised [XXX_local]
    diag[2] = 0.;
    pos = 0;
    for (process = 0; process < par_variables->mpi_world_size; process++)
        for (i = par_variables->layout_par_x.boxes[process].i_min; i <= par_variables->layout_par_x.boxes[process].i_max; i++)
            diag[2] += par_variables->recv_buf[pos++];
    diag[2] *= delta_x;
}

/*
 * Given the density functions fi(x, v) and fe(x, v) of ions and electrons,
 * compute the charge density rho = rhoi - rhoe = int_v (fi(x,v) - fe(x,v)) dv.
 * NB.: rho(e/i) could be allocated and destroyed by the function, they
 *      are temporary, but we give them to avoid frequent malloc / free of arrays.
 */
void update_rho(mesh_1d* meshx, mesh_1d* meshve, mesh_1d* meshvi,
        parallel_stuff* electrons, parallel_stuff* ions,
        double* rhoe, double* rhoi, double* rho,
	    bool is_periodic) {
	update_spatial_density(electrons, meshx->array, meshx->size, meshve->array, meshve->size, rhoe, is_periodic);
	update_spatial_density(ions,      meshx->array, meshx->size, meshvi->array, meshvi->size, rhoi, is_periodic);
	for (size_t i = 0; i < meshx->size; i++) {
        rho[i] = rhoi[i] - rhoe[i];
	}
}

/*
 * Given the density functions fi(x, v) and fe(x, v) of ions and electrons,
 * compute:
 *     >>> the charge density rho = rhoi - rhoe = int_v (fi(x,v) - fe(x,v)) dv.
 *     >>> the current current = currenti - currente = int_v (fi(x,v) - fe(x,v)) v dv.
 * NB.: rho(e/i) and current(e/i) could be allocated and destroyed by the function, they
 *      are temporary, but we give them to avoid frequent malloc / free of arrays.
 */
void update_rho_and_current(mesh_1d* meshx, mesh_1d* meshve, mesh_1d* meshvi,
        parallel_stuff* electrons, parallel_stuff* ions,
        double* Mass_e,
        double* rhoe, double* rhoi, double* rho,
        double* currente, double* currenti, double* current,
	    bool is_periodic) {
    update_rho(meshx, meshve, meshvi, electrons, ions, rhoe, rhoi, rho, is_periodic);
    *Mass_e = compute_mass(meshx->array, meshx->size, rhoe);
	update_current(electrons, meshx->array, meshx->size, meshve->array, meshve->size, currente, is_periodic);
	update_current(ions,      meshx->array, meshx->size, meshvi->array, meshvi->size, currenti, is_periodic);
	for (size_t i = 0; i < meshx->size; i++) {
        current[i] = currenti[i] - currente[i];
	}
}


#endif // ifndef SELA_VP_1D1V_CART_RHO
