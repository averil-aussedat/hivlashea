#ifndef SELA_VP_1D1V_CART_REMAP
#define SELA_VP_1D1V_CART_REMAP

#include <mpi.h>           // function  MPI_Alltoallv
                           // constants MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
#include <stdbool.h>       // type      bool
#include <stdio.h>         // function  printf
#include <stdlib.h>        // functions malloc, free ((de)allocate memory)
                           //           exit (error handling)
                           // constant  EXIT_FAILURE (error handling)
                           // type      size_t
#include "array_helpers.h" // functions allocate_int_matrix, deallocate_int_matrix, allocate_2d_array
#include "math_helpers.h"  // functions imin, imax

/*****************************************************************************
 *                           Parallel layouts                                *
 *                                                                           *
 * selalib/src/parallelization/remap/sll_m_remapper.F90                      *
 *****************************************************************************/

/*
 * Index limits contained in a given processor.
 */
typedef struct box_2d box_2d;
struct box_2d {
     int i_min, i_max; // indexes for x
     int j_min, j_max; // indexes for v
};

/*
 * Prints the values of the box.
 *
 * @param box the box to print.
 */
void print_box_2d(box_2d box) {
    printf("[%d ; %d] x [%d ; %d] \n", box.i_min, box.i_max,
        box.j_min, box.j_max);
}

/*
 * Test if a box is empty.
 *
 * @param  box the box on which to check emptiness.
 * @return true iff box is empty.
 */
bool is_empty_box_2d(box_2d box) {
    return (box.i_min > box.i_max) || (box.j_min > box.j_max);
}

/*
 * Return the number of different doublets in the box.
 *
 * @param  box the box on which to count the quadruplets.
 * @return the number of different quadruplets in the box.
 */
int number_elements_box_2d(box_2d box) {
    if (is_empty_box_2d(box))
        return 0;
    else
        return (box.i_max - box.i_min + 1) * (box.j_max - box.j_min + 1);
}

/*
 * Intersects two boxes.
 *
 * @param[in]  b1, b2 : the two boxes to intersect.
 * @param[out] ans : the intersection of the two boxes.
 * @return     false iff the intersection is empty.
 */
bool intersect_boxes_2d(box_2d b1, box_2d b2, box_2d* ans) {
    ans->i_min = imax(b1.i_min, b2.i_min);
    ans->i_max = imin(b1.i_max, b2.i_max);
    ans->j_min = imax(b1.j_min, b2.j_min);
    ans->j_max = imin(b1.j_max, b2.j_max);

    return !is_empty_box_2d(*ans);
}

/*
 * Information on an array of boxes that describes the distribution
 * of data among different MPI processes.
 * 
 * FIXME : The global dimensions of the given dataset distributed
 * are for now not used and can be removed.
 * => then there is no need for a struct, a layout can directly
 * be an array of boxes.
 * => or the local sizes, now stored in parallel_stuff, can be
 * stored inside the layout.
 */
typedef struct sll_t_layout_2d sll_t_layout_2d;
struct sll_t_layout_2d {
    int global_size_x; // Global phase space size in x
    int global_size_v; // Global phase space size in v
    box_2d* boxes;     // boxes[i] describes the portion of the full
                       //   phase space that processor i has in memory
};

/*
 * MPI variables needed to change from a layout to another,
 * plus boxes that represent, for each process :
 * send_boxes[p] : the intersection between what he has and what
 *                 process p wants.
 * recv_boxes[p] : the intersection between what he wants and what
 *                 process p has.
 */
typedef struct remap_plan_2d remap_plan_2d;
struct remap_plan_2d {
    int* send_counts;
    int* send_displs;
    double* send_buf;
    int* recv_counts;
    int* recv_displs;
    double* recv_buf;
    box_2d* send_boxes;
    box_2d* recv_boxes;
};

/*
 * This is a struct conveyed from function to function throughout the
 * simulation to store everything relative to MPI parallelism.
 */
typedef struct parallel_stuff parallel_stuff;
struct parallel_stuff {
    int mpi_world_size;           // Total number of MPI processes
    int mpi_rank;                 // Rank of the local MPI process
    sll_t_layout_2d layout_par_x; // How is the data spread among processors when parallel in x?
    sll_t_layout_2d layout_par_v; // How is the data spread among processors when parallel in v?
    remap_plan_2d remap_x_to_v;   // How to remap the data from par in x to par in v?
    remap_plan_2d remap_v_to_x;   // How to remap the data from par in v to par in x?
    int size_x_par_x;             // Local ncx size when parallel in x
    int size_v_par_x;             // Local ncv size when parallel in x
    int size_x_par_v;             // Local ncx size when parallel in v
    int size_v_par_v;             // Local ncv size when parallel in v
    double** f_parallel_in_x;     // Local values of f when parallel in x
    double** f_parallel_in_v;     // Local values of f when parallel in v
    bool is_par_x;                // What is the last updated f stored?
    //double* f_1d_in;              // array[max(ncx, ncv)] for Lagrange interpolations
    //double* f_1d_out;             // array[max(ncx, ncv)] for Lagrange interpolations
    double* f_1d;              	  // array[max(ncx, ncv)] for Lagrange interpolations
    double* send_buf;             // array[loc_ncx] for MPI_Allgatherv
    double* recv_buf;             // array[ncx]         for MPI_Allgatherv
    int* recv_counts;             // array[mpi_world_size] for MPI_Allgatherv
    int* displs;                  // array[mpi_world_size] for MPI_Allgatherv
    int global_indices[2];        // Every time we need to access global indices from local  ones
    int local_indices[2];         // Every time we need to access local  indices from global ones
};

/*
 * Converts the local indices
 *   0 <= i < local_size_x
 *   0 <= j < local_size_v
 * into global indices in the full phase space.
 */
void local_to_global_2d(parallel_stuff* par_variables,
        int i, int j) {
    sll_t_layout_2d layout = par_variables->is_par_x
        ? par_variables->layout_par_x
        : par_variables->layout_par_v;
    box_2d box = layout.boxes[par_variables->mpi_rank];
    par_variables->global_indices[0] = box.i_min + i;
    par_variables->global_indices[1] = box.j_min + j;
}

/*
 * Converts the global indices
 *   0 <= i < global_size_x
 *   0 <= j < global_size_v
 * into local indices in the local phase space.
 *
 * Assumes that the global indices are handled by this process
 * (if not, will return indices that will probably cause OutOfRangeException).
 */
void global_to_local_2d(parallel_stuff* par_variables,
        int i, int j) {
    sll_t_layout_2d layout = par_variables->is_par_x
        ? par_variables->layout_par_x
        : par_variables->layout_par_v;
    box_2d box = layout.boxes[par_variables->mpi_rank];
    par_variables->local_indices[0] = i - box.i_min;
    par_variables->local_indices[1] = j - box.j_min;
}

/*
 * Creation of a plan to remap from one layout to another.
 */
remap_plan_2d new_remap_plan_2d(int mpi_world_size, 
	int mpi_rank, sll_t_layout_2d initial_layout, sll_t_layout_2d final_layout) {
    int process;
    box_2d initial_box, final_box, intersect_box;
    
    remap_plan_2d plan;
    plan.send_counts = malloc(mpi_world_size * sizeof(int));
    plan.send_displs = malloc(mpi_world_size * sizeof(int));
    plan.recv_counts = malloc(mpi_world_size * sizeof(int));
    plan.recv_displs = malloc(mpi_world_size * sizeof(int));
    plan.send_boxes = malloc(mpi_world_size * sizeof(box_2d));
    plan.recv_boxes = malloc(mpi_world_size * sizeof(box_2d));
    
    initial_box = initial_layout.boxes[mpi_rank];
    for (process = 0; process < mpi_world_size; process++) {
        final_box = final_layout.boxes[process];
        intersect_boxes_2d(initial_box, final_box, &intersect_box);
        plan.send_counts[process] = number_elements_box_2d(intersect_box);
        plan.send_boxes[process]  = intersect_box;
    }
    plan.send_displs[0] = 0;
    for (process = 1; process < mpi_world_size; process++)
        plan.send_displs[process] = plan.send_displs[process - 1] + plan.send_counts[process - 1];
    plan.send_buf = malloc((plan.send_displs[mpi_world_size - 1] + plan.send_counts[mpi_world_size - 1]) * sizeof(double));

    final_box = final_layout.boxes[mpi_rank];
    for (process = 0; process < mpi_world_size; process++) {
        initial_box = initial_layout.boxes[process];
        intersect_boxes_2d(initial_box, final_box, &intersect_box);
        plan.recv_counts[process] = number_elements_box_2d(intersect_box);
        plan.recv_boxes[process]  = intersect_box;
    }
    plan.recv_displs[0] = 0;
    for (process = 1; process < mpi_world_size; process++)
        plan.recv_displs[process] = plan.recv_displs[process - 1] + plan.recv_counts[process - 1];
    plan.recv_buf = malloc((plan.recv_displs[mpi_world_size - 1] + plan.recv_counts[mpi_world_size - 1]) * sizeof(double));

    return plan;
}

/*
 * Do the remap as stored in plan.
 */
void apply_remap_2d(parallel_stuff* par_variables, remap_plan_2d* plan) {
    int process, i, j;
    box_2d box;
    double** data_in = par_variables->is_par_x
        ? par_variables->f_parallel_in_x
        : par_variables->f_parallel_in_v;
    double** data_out = par_variables->is_par_x
        ? par_variables->f_parallel_in_v
        : par_variables->f_parallel_in_x;
    int min_values[2];
    
    int loc = 0; // first loading is at position zero
    for (process = 0; process < par_variables->mpi_world_size; process++) {
        if (plan->send_counts[process] > 0) { // send something to rank 'process'
            if (loc != plan->send_displs[process]) {
                fprintf(stderr, "ERROR DETECTED in process %d. ", par_variables->mpi_rank);
                fprintf(stderr, "discrepancy between send_displs[%d]=%d and the loading index loc=%d.\n",
                    process, plan->send_displs[process], loc);
                exit(EXIT_FAILURE);
            }
            // get the information on the box to send, get the limits,
            // convert to the local indices and find out where in the
            // buffer to start writing.
            box = plan->send_boxes[process];
            global_to_local_2d(par_variables, box.i_min, box.j_min);
            for (i = 0; i < 2; i++)
                min_values[i] = par_variables->local_indices[i];
            global_to_local_2d(par_variables, box.i_max, box.j_max);
            
            // The plan to load the send buffer is to traverse the integer
            // array with a single index (loc). When we load the buffer, each
            // data element may occupy multiple integer 'slots', hence the
            // loading index needs to be manually increased. As an advantage,
            // we can do some error checking every time we send data to a
            // different process, as we know what is the expected value of
            // the index at that point.
            for (i = min_values[0] ; i <= par_variables->local_indices[0] ; i++)
                for (j = min_values[1] ; j <= par_variables->local_indices[1] ; j++)
                    plan->send_buf[loc++] = data_in[i][j];
        }
    }
    MPI_Alltoallv(plan->send_buf, plan->send_counts, plan->send_displs, MPI_DOUBLE_PRECISION,
                  plan->recv_buf, plan->recv_counts, plan->recv_displs, MPI_DOUBLE_PRECISION,
                  MPI_COMM_WORLD);
    
    par_variables->is_par_x = !par_variables->is_par_x;
    
    loc = 0; // We load first from position 0 in the receive buffer.
    for (process = 0; process < par_variables->mpi_world_size; process++) {
        if (plan->recv_counts[process] > 0) { // we expect something from rank 'i'
            if (loc != plan->recv_displs[process]) {
                fprintf(stderr, "ERROR DETECTED in process %d. ", par_variables->mpi_rank);
                fprintf(stderr, "discrepancy between recv_displs[%d]=%d and the loading index loc=%d.\n",
                    process, plan->recv_displs[process], loc);
                exit(EXIT_FAILURE);
            }
            // get the information on the box to receive, get the limits, and
            // convert to the local indices.
            box = plan->recv_boxes[process];
            global_to_local_2d(par_variables, box.i_min, box.j_min);
            for (i = 0; i < 2; i++)
                min_values[i] = par_variables->local_indices[i];
            global_to_local_2d(par_variables, box.i_max, box.j_max);
            for (i = min_values[0] ; i <= par_variables->local_indices[0] ; i++)
                for (j = min_values[1] ; j <= par_variables->local_indices[1] ; j++)
                    data_out[i][j] = plan->recv_buf[loc++];
        }
    }
}

/*
 * Change from parallel in x to parallel in v, or the opposite.
 */
void exchange_parallelizations(parallel_stuff* par_variables) {
    //double start_time = omp_get_wtime();
    if (par_variables->is_par_x)
        apply_remap_2d(par_variables, &par_variables->remap_x_to_v);
    else
        apply_remap_2d(par_variables, &par_variables->remap_v_to_x);
    //global_timings->time_remap += omp_get_wtime() - start_time;
}

int linear_index_2d(int ncx2,
        int i, int j) {
    return i * ncx2 + j;
}

/*
 * Auxiliary function that splits a range
 * of indices described by 2 integers into a given number of intervals, the
 * function tries to partition the original interval equitably.
 */
void split_array_indices_aux(int** intervals_array,
        int start_index, int interval_segment_length,
        int min, int max) {
    if (interval_segment_length == 1) {
        // terminate recursion by filling values for this interval
        intervals_array[0][start_index] = min;
        intervals_array[1][start_index] = max;
    } else {
        // split this interval and launch new recursions
        split_array_indices_aux(intervals_array,
            start_index,                             interval_segment_length/2, min, min + (max - min) / 2);
        split_array_indices_aux(intervals_array,
            start_index + interval_segment_length/2, interval_segment_length/2, min + (max - min) / 2 + 1, max);
    }
}

int** split_array_indices(int num_elements, int num_intervals) {
    if (num_elements < num_intervals) {
        fprintf(stderr, "Requested to split %d elements into %d intervals : impossible.\n",
            num_elements, num_intervals);
        exit(EXIT_FAILURE);
    }
    int** intervals = allocate_int_matrix(2, num_intervals);
    split_array_indices_aux(intervals, 0, num_intervals, 0, num_elements - 1);
    return intervals;
}

/*
 * Creates the initial layout, depending on the number of processors available
 * and the number of cells in each direction.
 *
 * @param[in] global_ncx  number of cells in the x direction
 * @param[in] global_ncv number of cells in the v direction
 * @param[in] num_proc_x  number of MPI processes in the x direction
 * @param[in] num_proc_v number of MPI processes in the v direction
 */
sll_t_layout_2d initialize_layout_with_distributed_2d_array(
        int global_ncx, int global_ncv, 
        int num_proc_x, int num_proc_v) {
    int total_num_processors = num_proc_x  * num_proc_v;
    sll_t_layout_2d layout;
    int node;
    box_2d box;

    layout.global_size_x = global_ncx;
    layout.global_size_v = global_ncv;
    layout.boxes = malloc(total_num_processors * sizeof(box_2d));

    int** intervals_x = split_array_indices(global_ncx,  num_proc_x);
    int** intervals_v = split_array_indices(global_ncv, num_proc_v);

    for (int i = 0; i < num_proc_x; i++) {
        for (int j = 0; j < num_proc_v; j++) {
                    node = linear_index_2d(num_proc_v, i, j);
                    box.i_min = intervals_x[0][i];
                    box.i_max = intervals_x[1][i];
                    box.j_min = intervals_v[0][j];
                    box.j_max = intervals_v[1][j];
                    layout.boxes[node] = box;
        }
    }
    
    deallocate_int_matrix(intervals_x,  2, num_proc_x);
    deallocate_int_matrix(intervals_v,  2, num_proc_v);
    
    return layout;
}



void init_par_variables(parallel_stuff* par_variables, 
	    int mpi_world_size, int mpi_rank, int sizex, int sizev,
	    bool is_periodic) {
    size_t i;
    (*par_variables).mpi_world_size = mpi_world_size;
    (*par_variables).mpi_rank       = mpi_rank;
    
    // For a periodic case, we may distribute sizex-1 points among the processes,
    // and construct the last point as equal to the point in 0.
    // For a non periodic case, we need to have the last index be sizex.
    int total_nb_points_on_x = is_periodic ? sizex - 1 : sizex;
    int total_nb_points_on_v = sizev - 1;
    (*par_variables).layout_par_x = initialize_layout_with_distributed_2d_array(
        total_nb_points_on_x, total_nb_points_on_v, mpi_world_size, 1);
    (*par_variables).layout_par_v = initialize_layout_with_distributed_2d_array(
        total_nb_points_on_x, total_nb_points_on_v, 1, mpi_world_size);
    
    box_2d box = (*par_variables).layout_par_x.boxes[(*par_variables).mpi_rank];
    (*par_variables).size_x_par_x  = box.i_max - box.i_min + 1;
    (*par_variables).size_v_par_x = box.j_max - box.j_min + 1;
    box = (*par_variables).layout_par_v.boxes[(*par_variables).mpi_rank];
    (*par_variables).size_x_par_v  = box.i_max - box.i_min + 1;
    (*par_variables).size_v_par_v  = box.j_max - box.j_min + 1;
	//printf("par_variables.size_x_par_x=%d (rank=%d)\n",(*par_variables).size_x_par_x,mpi_rank);	
	//printf("par_variables.size_v_par_x=%d (rank=%d)\n",(*par_variables).size_v_par_x,mpi_rank);	
    
    (*par_variables).f_parallel_in_x = allocate_2d_array( 
    	(*par_variables).size_x_par_x, (*par_variables).size_v_par_x);
    (*par_variables).f_parallel_in_v = allocate_2d_array(
    (*par_variables).size_x_par_v, (*par_variables).size_v_par_v);
    (*par_variables).remap_x_to_v = new_remap_plan_2d(mpi_world_size, mpi_rank,
        (*par_variables).layout_par_x, (*par_variables).layout_par_v);
    (*par_variables).remap_v_to_x = new_remap_plan_2d(mpi_world_size, mpi_rank,
        (*par_variables).layout_par_v, (*par_variables).layout_par_x);
    
    int size_f_1d = imax(total_nb_points_on_x, total_nb_points_on_v);
    //(*par_variables).f_1d_in = malloc(size_f_1d * sizeof(double));
    //(*par_variables).f_1d_out = malloc(size_f_1d * sizeof(double));
    (*par_variables).f_1d = malloc(size_f_1d * sizeof(double));
    
    (*par_variables).send_buf = allocate_1d_array((*par_variables).size_x_par_x);
    (*par_variables).recv_buf = allocate_1d_array(total_nb_points_on_x);
    (*par_variables).recv_counts = malloc(mpi_world_size * sizeof(int));
    for (i = 0; i < mpi_world_size; i++)
        (*par_variables).recv_counts[i] =
            (*par_variables).layout_par_x.boxes[i].i_max - (*par_variables).layout_par_x.boxes[i].i_min + 1;
    (*par_variables).displs = malloc(mpi_world_size * sizeof(int));
    (*par_variables).displs[0] = 0;
    for (i = 1; i < mpi_world_size; i++)
        (*par_variables).displs[i] = (*par_variables).displs[i-1] + (*par_variables).recv_counts[i-1];
}

#endif // ifndef SELA_VP_1D1V_CART_REMAP
