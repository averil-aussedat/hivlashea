#ifndef SELA_ARRAY_HELPERS
#define SELA_ARRAY_HELPERS

#include <stdlib.h>         // functions malloc, free ((de)allocate memory)
                            // type      size_t
#include "string_helpers.h" // macro     ERROR_MESSAGE

/*****************************************************************************
 *                           Matrix functions                                *
 *                                INT**                                      *
 *****************************************************************************/

/*
 * Matrix allocation - int** format.
 * @param[in] nbRow is the number of rows of the matrix to be allocated.
 * @param[in] nbCol is the number of cols of the matrix to be allocated.
 * @return    a newly allocated (nbRow x nbCol) matrix.
 */
int** allocate_int_matrix(int nbRow, int nbCol) {
    size_t i;
    int** a = malloc(nbRow * sizeof(int*));
    if (a == NULL) {
        ERROR_MESSAGE("allocate_int_matrix(%d, %d) : malloc error.\n", nbRow, nbCol);
    }
    for (i = 0; i < nbRow; i++) {
        a[i] = malloc(nbCol * sizeof(int));
        if (a[i] == NULL) {
            ERROR_MESSAGE("allocate_int_matrix(%d, %d) : malloc error.\n", nbRow, nbCol);
        }
    }
    return a;
}

/*
 * Matrix deallocation - int** format.
 * @param[in, out] a is the (nbRow x nbCol) matrix to deallocate.
 * @param[in]      nbRow is the number of rows of the matrix to be deallocated.
 * @param[in]      nbCol is the number of cols of the matrix to be deallocated.
 */
void deallocate_int_matrix(int** a, int nbRow, int nbCol) {
    size_t i;
    for (i = 0; i < nbRow; i++)
        free(a[i]);
    free(a);
}

/*****************************************************************************
 *                          2d-Array functions                               *
 *                             DOUBLE**                                      *
 *****************************************************************************/

/*
 * 2d-Array allocation - double** format.
 * @param[in] n1 is the number of cells in the dimension 1 of the 2d-Array to be allocated.
 * @param[in] n2 is the number of cells in the dimension 2 of the 2d-Array to be allocated.
 * @return    a newly allocated (n1 x n2) 2d-Array.
 */
double** allocate_2d_array(int n1, int n2) {
    size_t i;
    double** a = malloc(n1 * sizeof(double*));
    if (a == NULL) {
        ERROR_MESSAGE("allocate_2d_array(%d, %d) : malloc error.\n", n1, n2);
    }
    for (i = 0; i < n1; i++) {
        a[i] = malloc(n2 * sizeof(double));
        if (a[i] == NULL) {
            ERROR_MESSAGE("allocate_2d_array(%d, %d) : malloc error.\n", n1, n2);
        }
    }
    return a;
}

/*
 * 2d-Array deallocation - double** format.
 * @param[in, out] a is the (n1 x n2) 2d-Array to deallocate.
 * @param[in]      n1 is the number of cells in the dimension 1 of the 2d-Array to be deallocated.
 * @param[in]      n2 is the number of cells in the dimension 2 of the 2d-Array to be deallocated.
 */
void deallocate_2d_array(double** a, int n1, int n2) {
    size_t i;
    for (i = 0; i < n1; i++) {
        free(a[i]);
    }
    free(a);
}


/*****************************************************************************
 *                          3d-Array functions                               *
 *                             DOUBLE**                                      *
 *****************************************************************************/

/*
 * 2d-Array allocation - double** format.
 * @param[in] n1 is the number of cells in the dimension 1 of the 2d-Array to be allocated.
 * @param[in] n2 is the number of cells in the dimension 2 of the 2d-Array to be allocated.
 * @return    a newly allocated (n1 x n2) 2d-Array.
 */
double*** allocate_3d_array(int n1, int n2, int n3) {
    size_t i,j;
    double*** a = malloc(n1 * sizeof(double**));
    if (a == NULL) {
        ERROR_MESSAGE("allocate_3d_array(%d, %d, %d) : malloc error.\n", n1, n2, n3);
    }
    for (i = 0; i < n1; i++) {
        a[i] = malloc(n2 * sizeof(double*));
        if (a[i] == NULL) {
            ERROR_MESSAGE("allocate_3d_array(%d, %d, %d) : malloc error.\n", n1, n2, n3);
        }
        for (j = 0; j < n2; j++) {
        	a[i][j] = malloc(n3 * sizeof(double));
        	if (a[i][j] == NULL) {
            	ERROR_MESSAGE("allocate_3d_array(%d, %d, %d) : malloc error.\n", n1, n2, n3);
        	}
        }	
    }
    return a;
}

/*
 * 2d-Array deallocation - double** format.
 * @param[in, out] a is the (n1 x n2) 2d-Array to deallocate.
 * @param[in]      n1 is the number of cells in the dimension 1 of the 2d-Array to be deallocated.
 * @param[in]      n2 is the number of cells in the dimension 2 of the 2d-Array to be deallocated.
 */
void deallocate_3d_array(double*** a, int n1, int n2, int n3) {
    size_t i,j;
    for (i = 0; i < n1; i++) {
     	for (j = 0; j < n2; j++) {
        	free(a[i][j]);
    	}
       	free(a[i]);
    }
    free(a);
}



/*****************************************************************************
 *                          1d-Array functions                               *
 *                             DOUBLE*                                       *
 *****************************************************************************/

/*
 * 1d-Array allocation - double* format.
 * @param[in] n1 is the number of cells in the dimension 1 of the 2d-Array to be allocated.
 * @return    a newly allocated (n1) 1d-Array.
 */
double* allocate_1d_array(int n1) {
    double* a = malloc(n1 * sizeof(double));
    if (a == NULL) {
        ERROR_MESSAGE("allocate_1d_array(%d) : malloc error.\n", n1);
    }
    return a;
}

/*
 * 1d-Array deallocation - double* format.
 * @param[in, out] a is the (n1 ) 1d-Array to deallocate.
 * @param[in]      n1 is the number of cells in the dimension 1 of the 1d-Array to be deallocated.
 */
void deallocate_1d_array(double* a, int n1) {
    free(a);
}

/*****************************************************************************
 *                           Matrix functions                                *
 *                              DOUBLE**                                     *
 *****************************************************************************/

/*
 * Matrix allocation - double** format.
 * @param[in] nbRow is the number of rows of the matrix to be allocated.
 * @param[in] nbCol is the number of cols of the matrix to be allocated.
 * @return    a newly allocated (nbRow x nbCol) matrix.
 */
double** allocate_matrix(int nbRow, int nbCol) {
    size_t i;
    double** a = malloc(nbRow * sizeof(double*));
    if (a == NULL) {
        ERROR_MESSAGE("allocate_matrix(%d, %d) : malloc error.\n", nbRow, nbCol);
    }
    for (i = 0; i < nbRow; i++) {
        a[i] = malloc(nbCol * sizeof(double));
        if (a[i] == NULL) {
            ERROR_MESSAGE("allocate_matrix(%d, %d) : malloc error.\n", nbRow, nbCol);
        }
    }
    return a;
}




/*
 * Matrix deallocation - double** format.
 * @param[in, out] a is the (nbRow x nbCol) matrix to deallocate.
 * @param[in]      nbRow is the number of rows of the matrix to be deallocated.
 * @param[in]      nbCol is the number of cols of the matrix to be deallocated.
 */
void deallocate_matrix(double** a, int nbRow, int nbCol) {
    size_t i;
    for (i = 0; i < nbRow; i++)
        free(a[i]);
    free(a);
}

#endif // ifndef SELA_ARRAY_HELPERS
