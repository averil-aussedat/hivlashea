//#hivlashea
//
// Compilation and run:
// ./compile_hivlashea.sh
// cd build && ./hivlashea.out ../yaml/cemracs2022_2sp.yaml

#include <mpi.h>                // constants MPI_THREAD_FUNNELED, MPI_COMM_WORLD
                                // functions MPI_Init_thread, MPI_Comm_size, MPI_Comm_rank
#include <stdbool.h>            // type      bool
#include <stdio.h>              // functions printf, fprintf
#include <stdlib.h>             // type      size_t
#include "adv1d_periodic_lag.h" // function  adv1d_periodic_lag_compute_lag
#include "advect.h"             // type      adv1d_t
                                // functions adv1d_init, advection_x, advection_v
#include "diag.h"               // functions diag_energy, diag_f
#include "init_any_expr.h"      // function  fun_1d1v
#include "mesh_1d.h"            // type      mesh_1d
                                // function  mesh_1d_create
#include "parameter_reader.h"   // type      PC_tree_t
                                // functions PC_get, PC_parse_path
#include "poisson_direct.h"     // type      poisson_solver_direct
                                // functions new_poisson_solver_direct, compute_E_from_rho_1d
#include "remap.h"              // type      parallel_stuff
                                // functions exchange_parallelizations, local_to_global_2d
#include "rho.h"                // function  update_spatial_density
#include "split.h"              // type      sll_t_splitting_coeff
#include "string_helpers.h"     // macro     ERROR_MESSAGE

/*
 * Read the simulation-specific parameters values from the yaml file.
 */
void read_simulation_parameters(PC_tree_t conf, int mpi_rank, double* lambda, double* nu) {
    if (PC_get(conf, ".lambda")) {
        PC_double(PC_get(conf, ".lambda"), lambda);
    } else {
        ERROR_MESSAGE("#Error in %s: missing the lambda value.\n", conf->key);
    }
    if (PC_get(conf, ".nu")) {
        PC_double(PC_get(conf, ".nu"), nu);
    } else {
        ERROR_MESSAGE("#Error in %s: missing the nu value.\n", conf->key);
    }
    
    // Print the values.
    if (mpi_rank == 0) {
        printf("#Reading the simulation parameters.\n");
        printf("#lambda = %1.20g, nu = %1.20g.\n", *lambda, *nu);
    }
}

/*
 * Given a value target inside (x[0], x[N-1]), computes alpha and io such that
 * x[i0] + alpha = target and 0 <= alpha < delta_x.
 * TODO: this function is called multiple times with increasing values of
 *       "target", thus time could be saved because we could continue to
 *       increase values of i to test, and not start to 0 each time.
 *
 * @param[in] target : target value to reach.
 * @param[in] xmin   : minimum value of the mesh.
 * @param[in] xmax   : maximum value of the mesh.
 * @param[in] N      : number of cells in the mesh.
 * @param[out] i0    : integer part.
 * @param[out] alpha : fractional part.
 */
void no_boundary_lag_compute_i0_and_alpha(double target, mesh_1d* meshx, int* i0, double* alpha) {
    int i;
    i = 0;
    while (i < meshx->size && meshx->array[i] > target) {
        i++;
    }
    if (i == meshx->size) {
        ERROR_MESSAGE("#The target value is outside the mesh.\n");
    }
    *i0 = i;
    
    double delta_x = (meshx->array[meshx->size - 1] - meshx->array[0]) / (double)(meshx->size - 1);
    double a = (target - meshx->array[i]) / delta_x;
    *alpha = a;
}

/*
 * Solve the equation:
 *     f_i^{n+1} = f_i^n + nu * delta_t/2 * f_e
 * It uses Lagrange interpolation with d=1 (4 points, degree 3).
 */
void source_term(parallel_stuff* electrons, parallel_stuff* ions,
        mesh_1d* meshve, mesh_1d* meshvi,
        double nu, double dt, double* rho) {
    // Lagrange interpolation to interpolate values of the electron mesh
    // on the ion mesh. The electron velocity mesh is on [-500; 500] whereas the
    // ion velocity mesh is on [-10; 10] and have the same number of points, thus
    // we need to interpolate values of fe on the ion velocity mesh.
    int d = 1;
    double* lag = malloc((2*d+2)*sizeof(double));
    int i0;
    double alpha;
    
    if (!electrons->is_par_x) {
        exchange_parallelizations(electrons);
    }
    if (!ions->is_par_x) {
        exchange_parallelizations(ions);
    }
    local_to_global_2d(electrons, 0, 0);
    local_to_global_2d(ions, 0, 0);
    for (int i_x = 0; i_x < ions->size_x_par_x; i_x++) {
        for (int i_v = 0; i_v < ions->size_v_par_x; i_v++) {
            double target = meshvi->array[i_v];
            no_boundary_lag_compute_i0_and_alpha(target, meshve, &i0, &alpha);
            adv1d_periodic_lag_compute_lag(alpha, d, lag);
            for (int j = -d; j <= d+1; j++) {
                ions->f_parallel_in_x[i_x][i_v] += lag[j+d] * electrons->f_parallel_in_x[i_x][i0+j];
            }
        }
    }
    free(lag);
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
        double* rhoe, double* rhoi, double* rho,
        double* currente, double* currenti, double* current) {
	update_spatial_density(electrons, meshx->array, meshx->size, meshve->array, meshve->size, rhoe);
	update_spatial_density(ions,      meshx->array, meshx->size, meshvi->array, meshvi->size, rhoi);
	for (size_t i = 0; i < meshx->size; i++) {
        rho[i] = rhoi[i] - rhoe[i];
	}
	update_current(electrons, meshx->array, meshx->size, meshve->array, meshve->size, currente);
	update_current(ions,      meshx->array, meshx->size, meshvi->array, meshvi->size, currenti);
	for (size_t i = 0; i < meshx->size; i++) {
        current[i] = currenti[i] - currente[i];
	}
}

int main(int argc, char *argv[]) {
    // MPI parallelism
    int mpi_world_size, mpi_rank;
    int mpi_thread_support;
    // Electric field and charge and current
    double *rho, *rhoe, *rhoi;
    double *current, *currente, *currenti;
    double *E;
    // Splitting
    sll_t_splitting_coeff split;
    double delta_t;
    int num_iteration;
    // Advection
    adv1d_t *adv_xe,*adv_xi;
    adv1d_t *adv_ve,*adv_vi;
    // Other simulation parameters
    double lambda, nu;
    
    // Electric energy computations
    double ee;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_support);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    //printf("#mpi_rank=%d mpi_world_size=%d\n", mpi_rank, mpi_world_size);
    
    // Check that we correctly provide the yaml file
    if (argc != 2) {
        ERROR_MESSAGE("#Usage: %s file.yaml\n", argv[0]);
    }
    PC_tree_t conf = PC_parse_path(argv[1]);
    
    // Read the mesh
    mesh_1d meshx, meshve, meshvi;
    if (PC_get(conf, ".meshx")) {
        meshx = mesh_1d_create(PC_get(conf, ".meshx"), mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing meshx in %s\n", argv[1]);
    }
    if (PC_get(conf, ".meshve")) {
        meshve = mesh_1d_create(PC_get(conf, ".meshve"), mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing meshve in %s\n", argv[1]);
    }
    if (PC_get(conf, ".meshvi")) {
        meshvi = mesh_1d_create(PC_get(conf, ".meshvi"), mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing meshvi in %s\n", argv[1]);
    }
    
    // Read the initial condition for electrons
    parallel_stuff pare;
    init_par_variables(&pare, mpi_world_size, mpi_rank, meshx.size, meshve.size);
    if (PC_get(conf, ".f0e")) {
	    fun_1d1v(PC_get(conf, ".f0e"), &pare, meshx.array, meshve.array);
    } else {
        ERROR_MESSAGE("#Missing f0e in %s\n", argv[1]);
    }
    diag_f(&pare, 1, meshx, meshve, 0, "f0e");
    // Read the initial condition for ions
    parallel_stuff pari;
    init_par_variables(&pari, mpi_world_size, mpi_rank, meshx.size, meshvi.size);
    if (PC_get(conf, ".f0i")) {
	    fun_1d1v(PC_get(conf, ".f0i"), &pari, meshx.array, meshvi.array);
    } else {
        ERROR_MESSAGE("#Missing f0i in %s\n", argv[1]);
    }
    diag_f(&pari, 1, meshx, meshvi, 0, "f0i");
    // Read the simulation parameters
    if (PC_get(conf, ".simulation_parameters")) {
	    read_simulation_parameters(PC_get(conf, ".simulation_parameters"), mpi_rank, &lambda, &nu);
    } else {
        ERROR_MESSAGE("#Missing simulation_parameters in %s\n", argv[1]);
    }
    
    // Read the advection parameters
    if (PC_get(conf, ".time_parameters")) {
	    splitting(PC_get(conf, ".time_parameters"), &split);
	    delta_t = split.dt;
	    num_iteration = split.num_iteration;
    } else {
        ERROR_MESSAGE("#Missing time_parameters in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_xe")) {
	    adv1d_init(&adv_xe, PC_get(conf, ".adv_xe"), meshx.array, meshx.size, mpi_rank);
    } else if (PC_get(conf, ".adv_x")) {
	    adv1d_init(&adv_xe, PC_get(conf, ".adv_x"), meshx.array, meshx.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_x or adv_xe in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_xi")) {
	    adv1d_init(&adv_xi, PC_get(conf, ".adv_xi"), meshx.array, meshx.size, mpi_rank);
    } else if (PC_get(conf, ".adv_x")) {
	    adv1d_init(&adv_xi, PC_get(conf, ".adv_x"), meshx.array, meshx.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_x or adv_xi in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_ve")) {
	    adv1d_init(&adv_ve, PC_get(conf, ".adv_ve"), meshve.array, meshve.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_ve in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_vi")) {
	    adv1d_init(&adv_vi, PC_get(conf, ".adv_vi"), meshvi.array, meshvi.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_vi in %s\n", argv[1]);
    }
    
    // Create the electric field, charge arrays, and the poisson solver
    poisson_solver_direct solver = new_poisson_solver_direct(meshx.array, meshx.size, delta_t, lambda, nu);
    rho = allocate_1d_array(meshx.size);
    rhoe = allocate_1d_array(meshx.size);
    rhoi = allocate_1d_array(meshx.size);
    current = allocate_1d_array(meshx.size);
    currente = allocate_1d_array(meshx.size);
    currenti = allocate_1d_array(meshx.size);
    E = allocate_1d_array(meshx.size);
    
    // Compute electric field (at initial time)
    update_rho_and_current(&meshx, &meshve, &meshvi, &pare, &pari, rhoe, rhoi, rho, currente, currenti, current);
    compute_E_from_rho_1d(solver, rho, current, E);
    
    FILE* file_diag_energy = fopen("diag_ee.txt", "w");
    fprintf(file_diag_energy, "Time | Int(Ex^2)\n");
    for (int i_time = 0; i_time < num_iteration; i_time++) {
        diag_energy(E, meshx.array, meshx.size, &ee);
        if (mpi_rank == 0) {
            fprintf(file_diag_energy, "%1.20lg %1.20lg\n", ((double)i_time)*delta_t, ee);
        }
        // Strang splitting
        // Half time-step: advection in x
        //     d_t(f_{i/e}) + v*d_x(f_{i/e}) = 0
    	advection_x(&pare, adv_xe, 0.5*delta_t, meshve.array);
    	advection_x(&pari, adv_xi, 0.5*delta_t, meshvi.array);
        // Half time step: source term for ions
        //     f_i^{n+1} = f_i^n + nu * delta_t/2 * f_e
        source_term(&pare, &pari, &meshve, &meshvi, nu, 0.5*delta_t, rho);
		// Solve Poisson
		update_rho_and_current(&meshx, &meshve, &meshvi, &pare, &pari, rhoe, rhoi, rho, currente, currenti, current);
		compute_E_from_rho_1d(solver, rho, current, E);
		// Full time-step: advection in v
		//     d_t(f_i) - 1/mu*E*d_v(f_i) = 0
		//     d_t(f_i) + E*d_v(f_i) = 0
    	advection_v(&pare, adv_ve, delta_t, E);
    	advection_v(&pari, adv_vi, delta_t, E);
        // Half time step: source term for ions
        source_term(&pare, &pari, &meshve, &meshvi, nu, 0.5*delta_t, rho);
        // Half time-step: advection in x
    	advection_x(&pare, adv_xe, 0.5*delta_t, meshve.array);
    	advection_x(&pari, adv_xi, 0.5*delta_t, meshvi.array);
		// Solve Poisson
		update_rho_and_current(&meshx, &meshve, &meshvi, &pare, &pari, rhoe, rhoi, rho, currente, currenti, current);
		compute_E_from_rho_1d(solver, rho, current, E);
    }
    
    // Output the ions / electrons density functions at the end
    diag_f(&pare, 1, meshx, meshve, 0, "fe");
    diag_f(&pari, 1, meshx, meshvi, 0, "fi");
    
    // Be clean (-:
    fclose(file_diag_energy);
    MPI_Finalize();
    return 0;
}

