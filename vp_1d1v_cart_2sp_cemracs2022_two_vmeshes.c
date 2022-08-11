//#hivlashea
//
// Compilation and run:
// ./compile_hivlashea.sh
// cd build && ./hivlashea.out ../yaml/cemracs2022_2sp.yaml

#include <mpi.h>                    // constants MPI_THREAD_FUNNELED, MPI_COMM_WORLD
                                    // functions MPI_Init_thread, MPI_Comm_size, MPI_Comm_rank
#include <math.h>                   // function  floor
#include <stdbool.h>                // type      bool
#include <stdio.h>                  // functions printf, fprintf
#include <stdlib.h>                 // type      size_t
// #include "adv1d_non_periodic_lag.h" // function no_boundary_lag_compute_index_and_alpha
#include "advect.h"                 // types     adv1d_x_t, adv1d_v_t, adv1d_[non_]periodic_lag_t
                                    // functions adv1d_x_init, adv1d_y_init, advection_x, advection_v
                                    // functions adv1d_compute_lag, adv1d_compute_i0_and_alpha
#include "diag.h"                   // functions diag_energy, diag_f
#include "init_any_expr.h"          // function  fun_1d1v, fill_array_1d
#include "mesh_1d.h"                // type      mesh_1d
                                    // function  mesh_1d_create
#include "parameter_reader.h"       // type      PC_tree_t
                                    // functions PC_get, PC_parse_path
#include "poisson_direct.h"         // type      poisson_solver_direct
                                    // function  new_poisson_solver_direct
//#include "poisson_periodic_fft.h" // type      poisson_solver_periodic_fft
                                    // function  new_poisson_solver_periodic_fft
                                    // COMMONFUN update_E_from_rho_and_current_1d
#include "remap.h"                  // type      parallel_stuff
                                    // functions init_par_variables, exchange_parallelizations, local_to_global_2d
#include "rho.h"                    // functions update_spatial_density, update_current, compute_mass
                                    // update_rho, update_rho_and_current
#include "split.h"                  // type      sll_t_splitting_coeff
#include "string_helpers.h"         // macro     ERROR_MESSAGE
#include "show.h"                   // functions show_mesh1d, show_adv

// MACROS
#define VERBOSE True // comment to hide terminal output
#define PLOTS True // comment to remove plot diagnostics
#define OUTFOLDER "../data/" // where to write .dat and .hdf5
// #define OUTFOLDER "" // where to write .dat and .hdf5

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
 * Solve the equation:
 *     f_i^{n+1} = f_i^n + nu * dt * f_e
 * It uses Lagrange interpolation with d=1 (4 points, degree 3).
 */
void source_term(parallel_stuff* electrons, parallel_stuff* ions,
        mesh_1d* meshve, mesh_1d* meshvi,
        // double nu, double dt, double* rho) {
        double nu, double dt) {
    // Lagrange interpolation to interpolate values of the electron mesh
    // on the ion mesh. The electron velocity mesh is on [-500; 500] whereas the
    // ion velocity mesh is on [-10; 10] (the spatial meh is shared). Thus
    // we need to interpolate values of fe on the ion velocity mesh.
    int d = 1;
    double* lag = malloc((2*d+2)*sizeof(double));
    int index, index_v, i_x, i_v, j; // running indexes
    double alpha, target, tmp;
    double coeff = nu * dt;
    
    if (!electrons->is_par_x) {
        exchange_parallelizations(electrons);
    }
    if (!ions->is_par_x) {
        exchange_parallelizations(ions);
    }
    local_to_global_2d(electrons, 0, 0);
    local_to_global_2d(ions, 0, 0);
    for (i_x=0; i_x < ions->size_x_par_x; i_x++) {
        for (i_v=0; i_v < ions->size_v_par_x; i_v++) {
            target = meshvi->array[i_v];
            // no_boundary_lag_compute_index_and_alpha(target, meshve, &index, &alpha);
            tmp = (double)(meshve->size) * (target - meshve->min) / (meshve->max - meshve->min);
            index = (int)floor(tmp);
            alpha = tmp - (double)index;
            adv1d_compute_lag (alpha, d, lag);
            for (j = -d; j <= d+1; j++) {
                index_v = index+j;
                if (index_v < 0 || index_v >= electrons->size_v_par_x) {
                    ERROR_MESSAGE("#The electron velocity mesh is not big enough.\n");
                }
                ions->f_parallel_in_x[i_x][i_v] += lag[j+d] * coeff * electrons->f_parallel_in_x[i_x][index_v];
            }
        }
    }
    free(lag);
}

int main(int argc, char *argv[]) {
    // MPI parallelism
    int mpi_world_size, mpi_rank;
    int mpi_thread_support;
    // Electric field and charge and current
    double Mass_e;
    double *rho, *rhoe, *rhoi;
    double *current, *currente, *currenti;
    double *E;
    // Splitting
    sll_t_splitting_coeff split;
    double delta_t = 0.;
    int num_iteration = 0; 
    // Advection
    adv1d_x_t *adv_xe,*adv_xi;
    adv1d_v_t *adv_ve,*adv_vi;
    // Physical simulation parameters
    double lambda, nu, mu;
    
    // Outputs
    double ee=0.0; // energy
    double mm=0.0; // mass conservation evaluation (mi-me-2lambda^2*E(1))
    double diffM3 = 0.0; // int_{t=0}^tn M_3(t,1) - M_3(t,-1) dt
    // int plot_frequency=1; // 1=every loop
    // int plot_frequency=10; // 1=every loop
    // int plot_frequency=50; // 1=every loop
    int plot_frequency=200; // 1=every loop

    // local variables
    int itime=0;

    #ifdef VERBOSE
        printf("Welcome in vp_1d1v_cart_2sp_cemracs2022_two_vmeshes.\n");
    #endif
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_support);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    //printf("#mpi_rank=%d mpi_world_size=%d\n", mpi_rank, mpi_world_size);

    #ifdef VERBOSE
        printf("\n[Proc %d] Reading parameters...\n", mpi_rank);
    #endif
    
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
    
    // Read the simulation parameters [WARNING: before creating Poisson solver]
    if (PC_get(conf, ".simulation_parameters")) {
	    read_simulation_parameters(PC_get(conf, ".simulation_parameters"), mpi_rank, &lambda, &nu);
    } else {
        ERROR_MESSAGE("#Missing simulation_parameters in %s\n", argv[1]);
    }
    if (PC_get(conf, ".time_parameters")) {
	    splitting(PC_get(conf, ".time_parameters"), &split);
	    delta_t = split.dt;
	    num_iteration = split.num_iteration;
    } else {
        ERROR_MESSAGE("#Missing time_parameters in %s\n", argv[1]);
    }

    #ifdef VERBOSE
        printf("[Proc %d] Done reading parameters.\n", mpi_rank);
        printf("[Proc %d] Creating Poisson solver...\n", mpi_rank);
    #endif
    
    // Create the Poisson solver
    // WARNING: comment / uncomment the following line AND the include
    // at the beginning of this file to change Poisson solver.
    bool is_periodic = false;
    poisson_solver_direct solver = new_poisson_solver_direct(meshx.array, meshx.size, lambda, nu);
//    poisson_solver_periodic_fft solver = new_poisson_solver_periodic_fft(meshx.array, meshx.size);

    #ifdef VERBOSE
        printf("[Proc %d] Done creating Poisson solver.\n", mpi_rank);
        printf("[Proc %d] Initializing simulation...\n", mpi_rank);
    #endif

    // Read the initial condition for electrons
    parallel_stuff pare;
    init_par_variables(&pare, mpi_world_size, mpi_rank, meshx.size, meshve.size, is_periodic);
    if (PC_get(conf, ".f0e")) {
	    fun_1d1v(PC_get(conf, ".f0e"), &pare, meshx.array, meshve.array);
    } else {
        ERROR_MESSAGE("#Missing f0e in %s\n", argv[1]);
    }
    // Read the initial condition for ions
    parallel_stuff pari;
    init_par_variables(&pari, mpi_world_size, mpi_rank, meshx.size, meshvi.size, is_periodic);
    if (PC_get(conf, ".f0i")) {
	    fun_1d1v(PC_get(conf, ".f0i"), &pari, meshx.array, meshvi.array);
    } else {
        ERROR_MESSAGE("#Missing f0i in %s\n", argv[1]);
    }
    
    // Read the advection parameters
    if (PC_get(conf, ".adv_xe")) {
	    adv1d_x_init(&adv_xe, PC_get(conf, ".adv_xe"), meshx.array, meshx.size, mpi_rank);
    } else if (PC_get(conf, ".adv_x")) {
	    adv1d_x_init(&adv_xe, PC_get(conf, ".adv_x"), meshx.array, meshx.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_x or adv_xe in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_xi")) {
	    adv1d_x_init(&adv_xi, PC_get(conf, ".adv_xi"), meshx.array, meshx.size, mpi_rank);
    } else if (PC_get(conf, ".adv_x")) {
	    adv1d_x_init(&adv_xi, PC_get(conf, ".adv_x"), meshx.array, meshx.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_x or adv_xi in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_ve")) {
	    adv1d_v_init(&adv_ve, PC_get(conf, ".adv_ve"), meshve.array, meshve.size, mpi_rank);
        if (adv_ve->periodic_adv) {
            mu = -1.0 / adv_ve->periodic_adv->v;
        } else {
            mu = -1.0 / adv_ve->non_periodic_adv->v;
        }
    } else {
        ERROR_MESSAGE("#Missing adv_ve in %s\n", argv[1]);
    }
    if (PC_get(conf, ".adv_vi")) {
	    adv1d_v_init(&adv_vi, PC_get(conf, ".adv_vi"), meshvi.array, meshvi.size, mpi_rank);
    } else {
        ERROR_MESSAGE("#Missing adv_vi in %s\n", argv[1]);
    }


    #ifdef VERBOSE
        if (mpi_rank==0) {
            printf("------------- Parameters -------------\n");
            show_mesh1d (&meshx, "meshx");
            show_mesh1d (&meshvi, "meshvi");
            show_mesh1d (&meshve, "meshve");

            printf("Final time %f, %d iterations (delta t=%f).\n", split.dt*split.num_iteration, split.num_iteration, split.dt);
            printf("lambda = %f, nu = %f.\n", lambda, nu);

            // show_adv (adv_xe->non_periodic_adv, "adv_xe");
            // show_adv (adv_xi->non_periodic_adv, "adv_xi");
            // show_adv (adv_ve->non_periodic_adv, "adv_ve");
            // show_adv (adv_vi->non_periodic_adv, "adv_vi");
            printf("--------------------------\n");
        }
    #endif
    
    // Create the electric field, charge arrays, and the poisson solver
    rhoe = allocate_1d_array(meshx.size); // array of int_{v} f_e(t^n,x,v) for each x
    rhoi = allocate_1d_array(meshx.size); // array of int_{v} f_i(t^n,x,v) for each x
    rho = allocate_1d_array(meshx.size);  // rhoi - rhoe 
    currente = allocate_1d_array(meshx.size); // array of int_{v} v f_e(t^n,x,v)  for each x
    currenti = allocate_1d_array(meshx.size); // array of int_{v} v f_i(t^n,x,v)  for each x
    current = allocate_1d_array(meshx.size);  // currenti - currente
    E = allocate_1d_array(meshx.size); // electrical energy at time tn
    
    if (is_periodic) {
        // Compute electric field at initial time
        update_rho_and_current(&meshx, &meshve, &meshvi, &pare, &pari,
                &Mass_e, rhoe, rhoi, rho, currente, currenti, current, is_periodic);
        // update_E_from_rho_and_current_1d(solver, delta_t, Mass_e, rho, current, E);
        update_E_from_rho_and_current_1d (meshx.size-1, (meshx.max-meshx.min)/(meshx.size-1), solver.lambda, rho, E);

    } else {
        // Read electric field (at initial time)
        if (PC_get(conf, ".E0")) {
	        fill_array_1d(PC_get(conf, ".E0"), E, meshx.array, meshx.size);
        } else {
            ERROR_MESSAGE("#Missing E0 in %s\n", argv[1]);
        }
    }

    #ifdef PLOTS
        printf("[Proc %d] Plotting initial conditions...\n", mpi_rank);
        //diag_f(&pare, 0, meshx, meshve, 0.0, "fe", OUTFOLDER, is_periodic);
        //diag_f(&pari, 0, meshx, meshvi, 0.0, "fi", OUTFOLDER, is_periodic);
        diag_1d (E, meshx.array, meshx.size, "E", OUTFOLDER, 0, 0.0);
        diag_1d (rho, meshx.array, meshx.size, "rho", OUTFOLDER, 0, 0.0);
        printf("[Proc %d] Done plotting initial conditions.\n", mpi_rank);
    #endif 

    #ifdef VERBOSE
        printf("[Proc %d] Done initializing simulation.\n", mpi_rank);
        printf("[Proc %d] Beginning of the time loop.\n", mpi_rank);
    #endif
    
    FILE* file_diag_energy = fopen("diag_ee.txt", "w");
    fprintf(file_diag_energy, "Time | Int(Ex^2)\n");
    FILE* file_diag_mass   = fopen("diag_mm.txt", "w");
    fprintf(file_diag_mass, "Time | mi - me - lambda^2 * (E(1)+E(-1))\n");
    for (itime = 0; itime < num_iteration; itime++) {
        // #ifdef VERBOSE
        //     printf("[Proc %d] Time step %d / %d.\n", mpi_rank, itime, num_iteration);
        // #endif
        
        ///////////////////////////////////////////////////////////////////////////
        // // Original Strang splitting
    	// advection_x(&pari, adv_xi, 0.5*delta_t, meshvi.array); // d_t(f_i) + v*d_x(f_i) = 0
    	// advection_x(&pare, adv_xe, 0.5*delta_t, meshve.array); // d_t(f_e) + v*d_x(f_e) = 0
        // source_term(&pare, &pari, &meshve, &meshvi, nu, 0.5*delta_t); // d_t f_i = nu * f_e
        // update_rho_and_current(&meshx, &meshve, &meshvi, &pare, &pari,
        //         &Mass_e, rhoe, rhoi, rho, currente, currenti, current, is_periodic);
        // update_E_from_rho_and_current_1d(solver, delta_t, Mass_e, rho, current, E);
    	// advection_v(&pare, adv_ve, delta_t, E); // d_t(f_e) - 1/mu*E*d_v(f_e) = 0
    	// advection_v(&pari, adv_vi, delta_t, E); // d_t(f_i) +      E*d_v(f_i) = 0
        // source_term(&pare, &pari, &meshve, &meshvi, nu, 0.5*delta_t); // source term for ions
    	// advection_x(&pare, adv_xe, 0.5*delta_t, meshve.array); // advection in x
    	// advection_x(&pari, adv_xi, 0.5*delta_t, meshvi.array); // advection in x
        ///////////////////////////////////////////////////////////////////////////


        // Strang splitting

        // Half time-step        
        // printf("adv x i 1\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
    	advection_x(&pari, adv_xi, 0.5*delta_t, meshvi.array); // d_t(f_i) + v*d_x(f_i) = 0
        // printf("adv x e 1\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
        advection_x(&pare, adv_xe, 0.5*delta_t, meshve.array); // d_t(f_e) + v*d_x(f_e) = 0

		// Full time-step
        // printf("before rho/E update\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
        update_rho_and_current(&meshx, &meshve, &meshvi, &pare, &pari, &Mass_e, rhoe, rhoi, rho, currente, currenti, current, is_periodic);
        // update_E_from_rho_and_current_1d(solver, delta_t, Mass_e, rho, current, E); // lambda^2 d_x E = rho
        update_E_from_rho_and_current_1d (meshx.size-1, (meshx.max-meshx.min)/(meshx.size-1), solver.lambda, rho, E);

        // if (itime == num_iteration-1) {
        //     double err_dxE = 0.0, dx = (meshx.max-meshx.min)/(meshx.size-1);
        //     for (int i=0; i<meshx.size-1; ++i) {
        //         printf("Ei : %e, rho_i : %e, err : %e \n", E[i], rho[i], lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
        //         err_dxE = fmax(err_dxE, lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
        //     }
        //     printf("Err_dxE : %e\n", err_dxE);

        //     printf("E : ");
        //     for (int i=0; i<meshx.size; ++i) {
        //         printf("%e\t", E[i]);
        //         // printf("Ei : %e, rho_i : %e, err : %e \n", E[i], rho[i], lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
        //         // err_dxE = fmax(err_dxE, lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
        //     }
        //     printf("\n");

        //     printf("rho : ");
        //     for (int i=0; i<meshx.size; ++i) {
        //         printf("%e\t", rho[i]);
        //         // printf("Ei : %e, rho_i : %e, err : %e \n", E[i], rho[i], lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
        //         // err_dxE = fmax(err_dxE, lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
        //     }
        //     printf("\n");
        // }

        source_term(&pare, &pari, &meshve, &meshvi, nu, delta_t); // d_t f_i = nu * f_e
        // printf("adv v i\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
    	advection_v(&pari, adv_vi, delta_t, E); // d_t(f_i) + E*d_v(f_i) = 0
        // printf("adv v e\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
    	advection_v(&pare, adv_ve, delta_t, E); // d_t(f_e) - 1/mu*E*d_v(f_e) = 0

        // Second half time step 
        // printf("adv x e 2\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
    	advection_x(&pare, adv_xe, 0.5*delta_t, meshve.array); // advection in x
        // printf("adv x i 2\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
    	advection_x(&pari, adv_xi, 0.5*delta_t, meshvi.array); // advection in x

        // diag_energy(E, meshx.array, meshx.size, &ee); 
        diag_energy(&pari, &pare, meshx, meshvi, meshve, delta_t, E, lambda, mu, &diffM3, &ee);
        diag_mass_conservation (meshx, rhoi, rhoe, lambda, E, &mm); 
        if (mpi_rank == 0) {
            fprintf(file_diag_energy, "%1.20lg %1.20lg\n", ((double)itime+1)*delta_t, ee);
            fprintf(file_diag_mass, "%1.20lg %1.20lg\n", ((double)itime+1)*delta_t, mm);
        }

        // printf("End of time step :\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");


        #ifdef PLOTS
            if ((itime < num_iteration-1) && ((itime+1) % plot_frequency == 0)) {
                printf("[Proc %d] Plotting at time t=%.2f, itime=%d/%d.\n", mpi_rank, (itime+1)*delta_t, itime+1, num_iteration);
                //diag_f(&pare, itime+1, meshx, meshve, (itime+1)*delta_t, "fe", OUTFOLDER, is_periodic);
                //diag_f(&pari, itime+1, meshx, meshvi, (itime+1)*delta_t, "fi", OUTFOLDER, is_periodic);
                diag_1d (E, meshx.array, meshx.size, "E", OUTFOLDER, itime+1, (itime+1)*delta_t);
                diag_1d (rho, meshx.array, meshx.size, "rho", OUTFOLDER, itime+1, (itime+1)*delta_t);
                printf("energy : %e, mi - me - 2 lambda^2 * E(1) : %e\n", ee, mm);
            }
        #endif // ifdef PLOTS

        // printf("End of time step after plots :\n");
        // stats_1D (rhoe, meshx.size, "rhoe"); stats_1D (rhoi, meshx.size, "rhoi");
    }

    #ifdef VERBOSE
        printf("[Proc %d] End of the time loop.\n", mpi_rank);
    #endif
    
    #ifdef PLOTS
        printf("[Proc %d] Plotting terminal values...\n", mpi_rank);
        // Output the ions / electrons density functions at the end
        //diag_f(&pare, num_iteration, meshx, meshve, num_iteration*delta_t, "fe", OUTFOLDER, is_periodic);
        //diag_f(&pari, num_iteration, meshx, meshvi, num_iteration*delta_t, "fi", OUTFOLDER, is_periodic);
        diag_1d (E, meshx.array, meshx.size, "E", OUTFOLDER, num_iteration, num_iteration*delta_t);
        diag_1d (rho, meshx.array, meshx.size, "rho", OUTFOLDER, num_iteration, num_iteration*delta_t);
        printf("[Proc %d] Done plotting terminal values.\n", mpi_rank);
    #endif 

    // printf("E : ");
    // for (int i=0; i<meshx.size; ++i) {
    //     printf("%e\t", E[i]);
    //     // printf("Ei : %e, rho_i : %e, err : %e \n", E[i], rho[i], lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
    //     // err_dxE = fmax(err_dxE, lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
    // }
    // printf("\n");

    // printf("rho : ");
    // for (int i=0; i<meshx.size; ++i) {
    //     printf("%e\t", rho[i]);
    //     // printf("Ei : %e, rho_i : %e, err : %e \n", E[i], rho[i], lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
    //     // err_dxE = fmax(err_dxE, lambda*lambda*(E[i+1] - E[i])/dx - rho[i]);
    // }
    // printf("\n");
    
    // Be clean (-:
    fclose(file_diag_energy);
    fclose(file_diag_mass);
    MPI_Finalize();

    #ifdef VERBOSE
        printf("Bye\n");
    #endif
    return 0;
}

