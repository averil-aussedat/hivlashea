#ifndef SELA_VP_1D1V_CART_DIAG
#define SELA_VP_1D1V_CART_DIAG
#include "hdf5_io.h"        // function plot_f_cartesian_mesh_2d
#include "math_helpers.h"   // functions min, max
#include "remap.h"          // type parallel_stuff
                            // function exchange_parallelizations
#include "mesh_1d.h"        // type mesh_1d


/*
 * Computes int_{x,v} factor * v^momentum * f(x,v) dxdv.
 *
 * @param[in,out]   par_variables
 * @param[in]       meshx : spatial mesh
 * @param[in]       meshv : velocity mesh
 * @param[in]       factor : multiplies the integral.
 */
// double v_momentum (parallel_stuff* par_variables, mesh_1d meshx, mesh_1d meshv, double factor, int momentum) {
//     double energy = 0.0;

//     // to simplify
//     if (!par_variables->is_par_x) {
//         exchange_parallelizations(par_variables);
//     }

//     for (size_t i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
//         for (size_t i_v = 0; i_v < par_variables->size_v_par_x; ++i_v) {
//             energy += par_variables->f_parallel_in_x[i_x][i_v] * int_pow(meshv.array[i_v], momentum);
//         }
//     }
//     energy *= (meshv.max-meshv.min)/(meshv.size-1.0); // * dv
//     energy *= (meshx.max-meshx.min)/(meshx.size-1.0); // * dx
//     energy *= factor;

//     // sum over all processes
//     // MPI_AllReduce(&energy, &energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD); 
//     return energy;
// }

void diag_energy(double *E, double *x, int sizex, double *val){
    int i;
    int N;
    double dx;
    N = sizex - 1;
    dx = (x[N] - x[0])/(double)N;
    *val = 0.;
    for (i = 0; i < N; i++){
        *val += E[i]*E[i];
    }
    *val *=dx;
	*val = sqrt(*val);
    //printf("%lg %d\n", *val, sizex - 1);
}

void diag_mass_conservation (double* rhoi, double* rhoe, double lambda, double* E) {
    
}

/*
 * Save the 1d function E in a two-column file (x,E)
 *
 * @param[in] E : array of values
 * @param[in] x : array of locations (space mesh)
 * @param[in] sizex : number of points of the mesh x
 * @param[in] array_name : file name
 * @param[in] folder : saving folder 
 * @param[in] iplot : id of the plot (number)
 * @param[in] time : simulation physical time
 */
void diag_1d(double *func, double *x, int sizex, char* array_name, char* folder, int iplot, double time) {
    int i;
    FILE* file;
    char str[256];
    // int N = sizex - 1;
    // double dx = (x[N] - x[0])/(double)N;
    
    // char cplot[6]; // 4 digits + '\0' [+1 more to avoid a warning in sprintf]
    // if (iplot < 10) {
    //     sprintf(cplot, "000%d", iplot & 0xf); // 0xf = 15
    // } else if (iplot < 100) {
    //     sprintf(cplot, "00%d", iplot & 0x7f); // 0x7f = 127
    // } else if (iplot < 1000) {
    //     sprintf(cplot, "0%d", iplot & 0x3ff); // 0x3ff = 1023
    // } else {
    //     sprintf(cplot, "%d", iplot & 0x3fff); // 0x3fff = 16383
    // }
    // cplot[4] = 0;
    // sprintf(str,"%s%s%s.dat",folder,array_name,cplot);

    sprintf(str,"%s%s%06d.dat",folder,array_name,iplot);
    file=fopen(str,"w");
    fprintf(file,"%f\n", time);
    for (i=0;i<sizex;i++){
    	fprintf(file,"%1.20lg %1.20lg\n", x[i],func[i]);
    	//printf("%1.20lg %1.20lg\n", x[i],E[i]);
    }
    fclose(file);
}


void diag_f(parallel_stuff* par, int i_hdf5, mesh_1d mesh1, mesh_1d mesh2, 
		double time, char* array_name, char* folder, bool is_periodic){
	
	double (*f)[mesh2.size] = malloc((mesh1.size) * sizeof *f); // Array allocated contiguously for hdf5 outputs.
	//double *f = malloc((mesh1.size) * (mesh2.size) *sizeof(double)); // Array allocated contiguously for hdf5 outputs.
	double *f1d = malloc((par->size_x_par_x) * (mesh2.size) *sizeof(double)); // Array allocated contiguously for hdf5 outputs.
	// double *f_2 = malloc((mesh1.size-1) *sizeof(double)); // Array allocated contiguously for hdf5 outputs.
	// double *f1d_2 = malloc((par->size_x_par_x)  *sizeof(double)); // Array allocated contiguously for hdf5 outputs.
    int i,j;
    int pos;
    double minf=1e5, maxf=-1e5; // to be printed
    for (i=0;i<(par->size_x_par_x) * (mesh2.size);i++){
    	f1d[i] = 0.;
    }
    if (!par->is_par_x) {
        exchange_parallelizations(par);
    }
    for (i=0;i<par->mpi_world_size;i++) {
    	par->recv_counts[i] *=mesh2.size;
    	par->displs[i] *=mesh2.size;
    }    
    pos = 0;
    for (i = 0; i < par->size_x_par_x; i++) {
        for (j = 0; j < mesh2.size; j++) {
            f1d[pos++] = par->f_parallel_in_x[i][j];
            minf = min (minf, f1d[pos-1]);
            maxf = max (maxf, f1d[pos-1]);
        }
        if (is_periodic) {
            f1d[pos++] = par->f_parallel_in_x[i][0];
        }
        minf = min (minf, f1d[pos-1]);
        maxf = max (maxf, f1d[pos-1]);
    }
    printf("[diag_f] (min,max) %s : %e, %e.\n", array_name, minf, maxf);
//     printf("par->recv_counts=\n");
//     for(i=0;i<par->mpi_world_size;i++){
//     	printf("%d\n",par->recv_counts[i]);
//     }    
//     printf("par->displs=\n");
//     for(i=0;i<par->mpi_world_size;i++){
//     	printf("%d\n",par->displs[i]);
//     }    
//     printf("mesh1.size=%d\n",mesh1.size);
//     printf("mesh2.size=%d\n",mesh2.size);
//     printf("par->size_x_par_x=%d\n",par->size_x_par_x);
//MPI_Finalize();
//return 0;
    MPI_Allgatherv( //par->send_buf,
    	f1d,
    	//f1d_2,
    	//par->f_parallel_in_x,
    	//par->size_x_par_x, MPI_DOUBLE_PRECISION, 
     	par->size_x_par_x*mesh2.size, MPI_DOUBLE_PRECISION,
		f[0],
		//f_2,
		//par->recv_buf, 
		par->recv_counts, 
		par->displs, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD);
    for(i=0;i<par->mpi_world_size;i++){
    	par->recv_counts[i] /=mesh2.size;
    	par->displs[i] /=mesh2.size;
    }
	if (par->mpi_rank == 0){
		//plot_f_cartesian_mesh_2d(i_hdf5, f[0], mesh1, mesh2, time, array_name, folder);
        if (is_periodic) {
            for (i=0;i< (mesh2.size);i++){ 
                f[mesh1.size-1][i] = f[0][i];
            }		
        }
//     	for (i=0;i< (mesh2.size);i++){
//     	for (j=0;j< (mesh1.size);j++){
//     		printf("%d %d %1.20lg\n",j,i,f[j][i]);
//     	}}	
		plot_f_cartesian_mesh_2d(i_hdf5, f[0], mesh1, mesh2, time, array_name, folder);
	}		

    free (f);
    free (f1d);
}


#endif // ifndef SELA_VP_1D1V_CART_DIAG
