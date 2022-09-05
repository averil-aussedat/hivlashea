#ifndef SELA_VP_1D1V_CART_DIAG
#define SELA_VP_1D1V_CART_DIAG
#include "hdf5_io.h"        // function plot_f_cartesian_mesh_2d
#include "math_helpers.h"   // functions min, max
#include "remap.h"          // type parallel_stuff
                            // function exchange_parallelizations
#include "mesh_1d.h"        // type mesh_1d
#include "rho.h"            // function compute_mass

/*
 * Computes int_{x,v} factor * v^momentum * f(x,v) dxdv.
 *
 * @param[in,out]   par_variables
 * @param[in]       meshx : spatial mesh
 * @param[in]       meshv : velocity mesh
 */
double v_momentum (parallel_stuff* par_variables, mesh_1d meshx, mesh_1d meshv, int momentum) {
    double energy = 0.0;

    // to simplify
    if (!par_variables->is_par_x) {
        exchange_parallelizations(par_variables);
    }

    for (size_t i_x = 0; i_x < par_variables->size_x_par_x; i_x++) {
        for (size_t i_v = 0; i_v < par_variables->size_v_par_x; ++i_v) {
            energy += par_variables->f_parallel_in_x[i_x][i_v] * pow(meshv.array[i_v], momentum);
        }
    }
    energy *= (meshv.max-meshv.min)/(meshv.size-1.0); // * dv
    energy *= (meshx.max-meshx.min)/(meshx.size-1.0); // * dx

    // sum over all processes
    MPI_Allreduce(&energy, &energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD); 
    return energy;
}

// not sure that it works
double get_diffM3 (parallel_stuff* par_variables, mesh_1d meshx, mesh_1d meshv) {
    double M3diff = 0.0, vv=0.0; 

    if (!par_variables->is_par_x) {
        exchange_parallelizations(par_variables);
    }
    
    // si le processus local contient x=xmin : + \int_v f(t,1) v^3 dv
    if (par_variables->global_indices[0]==0) { 
        for (size_t i_v = 0; i_v < par_variables->size_v_par_x; ++i_v) {
            vv = meshv.array[i_v];
            M3diff += par_variables->f_parallel_in_x[0][i_v] * vv*vv*vv;
        }
    }
    // si le processus local contient x=xmax : - \int_v f(t,-1) v^3 dv
    // if (local_to_global_2d(par_variables, par_variables->size_x_par_x-1, 0)==meshx.size-1) { 
    if (par_variables->global_indices[1]==meshx.size-1) { 
        for (size_t i_v = 0; i_v < par_variables->size_v_par_x; ++i_v) {
            vv = meshv.array[i_v];
            M3diff -= par_variables->f_parallel_in_x[par_variables->size_x_par_x-1][i_v] * vv*vv*vv;
        }
    }
    M3diff *= (meshv.max-meshv.min)/(meshv.size-1.0);

    // set everyone's M3diff to \int_v f(t,1) v^3 - \int_v f(t,-1) v^3 dv
    MPI_Allreduce(&M3diff, &M3diff, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD); 
    return M3diff;    
}

/*
 * Computes the quantity
 *      E_kin(t) + E_pot(t) + \int_{s=0^t} M3(s,1) - M3(s,-1) ds, 
 * with
 *      E_kin(t) = \int_{x,v} v^2/2 (fi(t,x,v) + mu )
 * 
 */
void diag_energy(parallel_stuff* pari, parallel_stuff* pare, mesh_1d meshx, mesh_1d meshvi, mesh_1d meshve, 
        double dt, double *E, double lambda, double mu, double* diffM3, double *val){
    int i;
    double dx = (meshx.max-meshx.min)/(meshx.size-1.0);
    double Epot=0.0, Ekin_i=0.0, Ekin_e=0.0;
    double M3_i=0.0, M3_e=0.0;

    // potential energy : E_pot = lambda^2 / 2 * \int_{-1}^{1} E(t,x) dx
    Epot += E[0]*E[0] * 0.5;
    for (i = 1; i < meshx.size-1; i++){
        *val += E[i]*E[i];
    }
    Epot += E[meshx.size-1]*E[meshx.size-1]*0.5;
    Epot *= dx * 2 * lambda*lambda;

    // kinetic energy : E_kin = \int_{x,v} v^2/2 (fi - mu fe) dxdv
    Ekin_i = 1. / 2.0 * v_momentum (pari, meshx, meshvi, 2);
    Ekin_e = mu / 2.0 * v_momentum (pare, meshx, meshve, 2);

    // computation of M3(t,+-1) := \int_{v} v^3/2 (f_i(t,+-1,v) + mu * fe(t,+-1,v)) dv
    M3_i = get_diffM3 (pari, meshx, meshvi);
    M3_e = get_diffM3 (pare, meshx, meshve);
    *diffM3 += dt * 0.5 * (M3_i + mu * M3_e);

    *val = Ekin_e + Ekin_i + *diffM3; 
}

void diag_mass_conservation (mesh_1d meshx, double* rhoi, double* rhoe, double lambda, double* E, double* mm) {
    double mi = compute_mass (meshx.array, meshx.size, rhoi);
    double me = compute_mass (meshx.array, meshx.size, rhoe);
    *mm = mi - me - 2.0*lambda*lambda * E[meshx.size-1];
    // printf("mi : %e, me : %e, res : %e\n", mi, me, *mm);
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


void read_f(double **f,int size1,int size2,char * array_name){

	int i,j;
	FILE* file;
	double res;
	//printf("size=%d %d\n",size1,size2);
	file = fopen(array_name,"r");
	
	res = 0.;
	
	for(i=0;i<size1;i++){
		for(j=0;j<size2;j++){
			//printf("%1.20lg\n",res);
			fscanf(file,"%lf",&res);
			f[i][j] = res;
			//printf("%1.20lg\n",res);
		}
	}	
	fclose(file);
	//printf("size=%d %d\n",size1,size2);

}


void read_f_par(parallel_stuff* par, int size1, int size2,char* array_name){	
    if (!par->is_par_x) {
        exchange_parallelizations(par);
    }
	read_f(par->f_parallel_in_x,size1,size2,array_name);
    if (!par->is_par_x) {
        exchange_parallelizations(par);
    }
}


void diag_f(parallel_stuff* par, int i_hdf5, mesh_1d mesh1, mesh_1d mesh2, 
		double time, char* array_name, char* folder, bool is_periodic){
	
	double (*f)[mesh2.size] = malloc((mesh1.size) * sizeof(*f)); // Array allocated contiguously for hdf5 outputs.
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
    if(is_periodic){
    	printf("is_periodic\n");
    }else{
    	printf("is not periodic\n");
    }
    printf("par->size_x_par_x=%d is_periodic=%d\n",par->size_x_par_x,is_periodic);
    for (i = 0; i < par->size_x_par_x; i++) {
        for (j = 0; j < mesh2.size; j++) {
            f1d[pos++] = par->f_parallel_in_x[i][j];
            minf = min (minf, f1d[pos-1]);
            maxf = max (maxf, f1d[pos-1]);
        }
        if (is_periodic) {
            f1d[pos++] = par->f_parallel_in_x[i][0];
        }
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

}


#endif // ifndef SELA_VP_1D1V_CART_DIAG
