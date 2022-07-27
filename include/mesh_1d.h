#ifndef SELA_MESH_1D
#define SELA_MESH_1D

#include <stdio.h>            // function  printf
#include <stdlib.h>           // function  malloc (allocate memory)
#include "parameter_reader.h" // type      PC_tree_t
                              // functions PC_get, PC_double, PC_int
#include "string_helpers.h"   // macro     ERROR_MESSAGE

typedef struct mesh_1d mesh_1d;
struct mesh_1d {
    int size;      // number of points on the axis
                   // N.B.: this means that the mesh has (size - 1) cells
    double* array; // array containing size values
    double min;
    double max;
};

mesh_1d mesh_1d_create_unif(PC_tree_t conf, int mpi_rank, char* mesh_name) {
    double min, max;
    long N;
    
    // Read the tree created from the yaml file.
    if (PC_get(conf, ".min")) {
        PC_double(PC_get(conf, ".min"), &min);
    } else {
        ERROR_MESSAGE("#Missing the min value in the uniform mesh1d %s.\n", mesh_name);
    }
    if (PC_get(conf, ".max")) {
        PC_double(PC_get(conf, ".max"), &max);
        if (max <= min) {
            ERROR_MESSAGE("#Wrong (min, max) values in the uniform mesh1d %s: max should be > min.\n", mesh_name);
        }
    } else {
        ERROR_MESSAGE("#Missing the max value in the uniform mesh1d %s.\n", mesh_name);
    }
    if (PC_get(conf, ".N")) {
        PC_int(PC_get(conf, ".N"), &N);
        if (N <= 0) {
            ERROR_MESSAGE("#Wrong N value in the uniform mesh1d %s: N should be > 0.\n", mesh_name);
        }
    } else {
        ERROR_MESSAGE("#Missing the N value in the uniform mesh1d %s.\n", mesh_name);
    }
    
    // Creation of the mesh.
    mesh_1d* mesh_p = malloc(sizeof(mesh_1d));
    mesh_p->size = (int)(N + 1);
    mesh_p->array = malloc(mesh_p->size * sizeof(double));
    double delta = (max - min) / (double)N;
    for (int i = 0; i < mesh_p->size; i++){
        mesh_p->array[i] = min + ((double)i) * delta;
    }
    mesh_p->min  = min;
    mesh_p->max  = max;
    // Print the values.
    if (mpi_rank == 0) {
        printf("#Reading the mesh1d %s.\n", mesh_name);
        printf("#min = %1.20g, max = %1.20g, N = %ld.\n", min, max, N);
        printf("#array = [%1.20lg, %1.20lg, ...]\n", mesh_p->array[0], mesh_p->array[1]); // Possible since N > 0.
    }
    
    return *mesh_p;
}


int gausspt(int d, double *res){
//d:=0:L:=[fsolve(evalf(expand(LegendreP(d+1,x))))];map(x->0.5*(1+x),L);
    //double* *res;
    //*res = malloc(sizeof(double)*(d+1));
    switch (d)
    {
    	case 0:
    		res[0] = 0.5;
    		break;
    	case 1:
    		res[0] = 0.21132486540518711774;
    		res[1] = 0.78867513459481288225;
    		break;
    	case 2:
    		res[0] = 0.11270166537925831148;
    		res[1] = 0.5;
    		res[2] = 0.88729833462074168850;
    		break;
		default:
            ERROR_MESSAGE("#Wrong d %d: d should be >= 0. and <=2\n", d);
			
    }
	return 0;
}


mesh_1d mesh_1d_create_gaussunif(PC_tree_t conf, int mpi_rank, char* mesh_name) {
    double min, max;
    long N;
    long d;
    double *gauss;
    // Read the tree created from the yaml file.
    if (PC_get(conf, ".min")) {
        PC_double(PC_get(conf, ".min"), &min);
    } else {
        ERROR_MESSAGE("#Missing the min value in the uniform mesh1d %s.\n", mesh_name);
    }
    if (PC_get(conf, ".max")) {
        PC_double(PC_get(conf, ".max"), &max);
        if (max <= min) {
            ERROR_MESSAGE("#Wrong (min, max) values in the uniform mesh1d %s: max should be > min.\n", mesh_name);
        }
    } else {
        ERROR_MESSAGE("#Missing the max value in the uniform mesh1d %s.\n", mesh_name);
    }
    if (PC_get(conf, ".N")) {
        PC_int(PC_get(conf, ".N"), &N);
        if (N <= 0) {
            ERROR_MESSAGE("#Wrong N value in the uniform mesh1d %s: N should be > 0.\n", mesh_name);
        }
    } else {
        ERROR_MESSAGE("#Missing the N value in the uniform mesh1d %s.\n", mesh_name);
    }

    if (PC_get(conf, ".d")) {
        PC_int(PC_get(conf, ".d"), &d);
        if (d <= -1) {
            ERROR_MESSAGE("#Wrong d value in the uniform mesh1d %s: d should be > 0.\n", mesh_name);
        }
    } else {
        ERROR_MESSAGE("#Missing the d value in the uniform mesh1d %s.\n", mesh_name);
    }
	gauss = malloc(sizeof(double)*(d+1));
    gausspt(d,gauss);

    printf("gauss=%1.20lg d=%ld\n",gauss[0],d);
    // Creation of the mesh.
    mesh_1d* mesh_p = malloc(sizeof(mesh_1d));
    mesh_p->size = (int)(N*(d+1));
    mesh_p->array = malloc(mesh_p->size * sizeof(double));
    double delta = (max - min) / (double)N;
    int k =0;
    
   for (int i = 0; i < N; i++){
    	for (int j = 0; j < d+1; j++){
        	mesh_p->array[k] = min + ((double)i+gauss[j]) * delta;
        	k++;
        	//printf("k=%d mesh_p->size=%d\n",k,mesh_p->size);
        }	
    }
    mesh_p->min = min;
    mesh_p->max = max;
    //exit(-1);
    // Print the values.
    if (mpi_rank == 0) {
        printf("#Reading the mesh1d %s.\n", mesh_name);
        printf("#min = %1.20g, max = %1.20g, N = %ld.\n", min, max, N);
        printf("#array = [%1.20lg, %1.20lg, ...]\n", mesh_p->array[0], mesh_p->array[1]); // Possible since N > 0.
    }
    return *mesh_p;
}




mesh_1d mesh_1d_create(PC_tree_t conf, int mpi_rank) {
    if (PC_get(conf, ".unif")) {
        return mesh_1d_create_unif(PC_get(conf, ".unif"), mpi_rank, conf->key);
    }else if (PC_get(conf, ".gaussunif")) {
        return mesh_1d_create_gaussunif(PC_get(conf, ".gaussunif"), mpi_rank, conf->key);
    } else {
        ERROR_MESSAGE("#Bad mesh1d type for %s. The only supported type is unif.\n", conf->key);
    }
}

#endif // ifndef SELA_MESH_1D
