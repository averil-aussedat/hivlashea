#ifndef SELA_HDF5_IO
#define SELA_HDF5_IO

#include <stdio.h>          // functions fprintf, fopen, fclose
#include "compiler_test.h"  // constant  SELA_C11
#include "hdf5.h"           // types     hid_t, herr_t, H5D_layout_t, hsize_t
                            // constants H5F_ACC_TRUNC, H5P_DEFAULT, H5D_CONTIGUOUS, H5T_IEEE_F64BE, H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE
                            //           H5S_ALL, H5P_DATASET_CREATE
                            // functions H5Fcreate, H5Screate_simple, H5Pset_layout, H5Pclose
#include "mesh_1d.h"        // type      mesh_1d
#include "string_helpers.h" // constant  DBL_DECIMAL_DIG
#include "variadic.h"       // macros    VARIADIC, NUMARG32

// Big Endian or Little Endian
// https://stackoverflow.com/questions/2100331/c-macro-definition-to-determine-big-endian-or-little-endian-machine
#define IS_BIG_ENDIAN (!*(unsigned char *)&(uint16_t){1})

/*
 * Create a new HDF5 file using the default properties.
 */
hid_t hdf5_ser_file_create(char* filename) {
    return H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

/*
 * Terminate access to an HDF5 file.
 */
void hdf5_ser_file_close(hid_t file_id) {
    H5Fclose(file_id);
}

/*
 * Write a 2d array of double in an HDF5 file.
 */
void hdf5_ser_write_dble_array_2d(hid_t file_id, hsize_t array_dims[2], const void* array, char* dsetname) {
    // Create dataspace. Setting maximum size to (void*)0 sets the maximum size to be the current size.
    hid_t dataspace_id        = H5Screate_simple(2, array_dims, (void*)0);
    // Create the dataset creation property list, set the layout to contiguous (compact doesn't work).
    hid_t dc_property_list_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(dc_property_list_id, H5D_CONTIGUOUS);
    // Create the dataset using the default properties.
    hid_t double_type_id      = IS_BIG_ENDIAN ? H5T_IEEE_F64BE : H5T_IEEE_F64LE;
    hid_t dataset_id          = H5Dcreate(file_id, dsetname, double_type_id, dataspace_id,
        H5P_DEFAULT, dc_property_list_id, H5P_DEFAULT);
    // Write the data to the dataset.
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
    // Close and release resources.
    H5Sclose(dataspace_id);
    H5Pclose(dc_property_list_id);
    H5Dclose(dataset_id);
}

/*
 * Create the XML file and begin to write first lines.
 *
 * @return a pointer to the created file.
 */
FILE* xml_file_create(char* file_name) {
    FILE* output_file = fopen(file_name, "w");
    fprintf(output_file, "<?xml version='1.0' ?>\n");
    fprintf(output_file, "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n");
    fprintf(output_file, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n");
    fprintf(output_file, "<Domain>\n");
    return output_file;
}

/*
 * Write the description of a scalar field on a 2D mesh.
 *
 * The file named filename must exist.
 *
 * @param[in] output_file pointer to the xml file
 * @param[in] filename    file name where the heavy data are 
 * @param[in] nnodes_x1   number of nodes number along direction 1
 * @param[in] nnodes_x2   number of nodes number along direction 2
 */
void xml_dataitem_2d(FILE* output_file, char* filename, int nnodes_x1, int nnodes_x2) {
    fprintf(output_file, "<DataItem Dimensions='%d %d' NumberType='Float' Precision='8' Format='HDF'>\n", nnodes_x1, nnodes_x2);
    fprintf(output_file, "%s\n", filename);
    fprintf(output_file, "</DataItem>\n");
}

/*
 * Write the description of a 2D strutured grid
 * mesh with its nodes coordinates contains in filename-x1 and filename-x2.
 *
 * The file named x*filename-x*dsetname.h5 with dataset x*dsetname must exist.
 * Low level version where you have to set dataset names in hdf5 files
 *
 * @param[in] output_file pointer to the xml file
 * @param[in] x1filename  file name where the coordinates x1 are
 * @param[in] x2filename  file name where the coordinates x2 are
 * @param[in] x1dsetname  dataset name of coordinates x1 are
 * @param[in] x2dsetname  dataset name of coordinates x2 are
 * @param[in] nnodes_x1   number of nodes number along direction 1
 * @param[in] nnodes_x2   number of nodes number along direction 2
 * @param[in] gridtype    "Uniform" or "Collection"
 */
void xml_grid_geometry_2d_low_level(FILE* output_file, char* x1filename, int nnodes_x1,
        char* x2filename, int nnodes_x2, char* x1dsetname, char* x2dsetname, char* gridtype) {
    fprintf(output_file, "<Grid Name='mesh' GridType='%s'>\n", gridtype);
    fprintf(output_file, "<Topology TopologyType='2DSMesh' NumberOfElements='%d %d'/>\n", nnodes_x1, nnodes_x2);
    fprintf(output_file, "<Geometry GeometryType='X_Y'>\n");
    char dataset_name[111];
    sprintf(dataset_name, "%s:/%s", x1filename, x1dsetname);
    xml_dataitem_2d(output_file, dataset_name, nnodes_x1, nnodes_x2);
    sprintf(dataset_name, "%s:/%s", x2filename, x2dsetname);
    xml_dataitem_2d(output_file, dataset_name, nnodes_x1, nnodes_x2);
    fprintf(output_file, "</Geometry>\n");
}

/*
 * Write the description of a 2D strutured grid
 * mesh with its nodes coordinates contains in mesh_name-x1 and mesh_name-x2.
 *
 * The file named filename.h5 with dataset x1 and x2 must exist.
 *
 * @param[in] output_file pointer to the xml file
 * @param[in] mesh_name   file name where the coordinates data are
 * @param[in] nnodes_x1   number of nodes number along direction 1
 * @param[in] nnodes_x2   number of nodes number along direction 2
 */
void xml_grid_geometry_2d_high_level(FILE* output_file, char* mesh_name, int nnodes_x1, int nnodes_x2) {
    char mesh_name_x1[117];
    char mesh_name_x2[117];
    sprintf(mesh_name_x1, "%s-x1.h5", mesh_name);
    sprintf(mesh_name_x2, "%s-x2.h5", mesh_name);
    xml_grid_geometry_2d_low_level(output_file,
                                   mesh_name_x1, nnodes_x1,
                                   mesh_name_x2, nnodes_x2,
                                   "x1", "x2", "Uniform");
}

/*
 * Open an XDMF format file for a 2d plot.
 *
 * @param[in] mesh_name file name where the coordinates data are
 * @param[in] nnodes_x1 number of nodes number along direction 1
 * @param[in] nnodes_x2 number of nodes number along direction 2
 */
FILE* xdmf_open_2d(char* file_name, char* mesh_name, int nnodes_x1, int nnodes_x2) {
    FILE* output_file = xml_file_create(file_name);
    xml_grid_geometry_2d_high_level(output_file, mesh_name, nnodes_x1, nnodes_x2);
    return output_file;
}

/*
 * Write the description of a scalar field on a 2D mesh.
 *
 * The file named filename.h5 with dataset fieldname must exist.
 *
 * @param[in] output_file pointer to the xml file
 * @param[in] fieldname   dataset name where the heavy data are 
 * @param[in] filename    file name where the heavy data are 
 * @param[in] npoints_1   number of nodes or cells number along direction 1
 * @param[in] npoints_2   number of nodes or cells number along direction 2
 * @param[in] center      values are centered on ("Node" or "Cell")
 */
void xml_field_2d(FILE* output_file, char* fieldname, char* filename,
        int npoints_1, int npoints_2, char* center) {
    fprintf(output_file, "<Attribute Name='%s' AttributeType='Scalar' Center='%s'>\n", fieldname, center);
    xml_dataitem_2d(output_file, filename, npoints_1, npoints_2);
    fprintf(output_file, "</Attribute>\n");
}

/*
 * Write 2d array in hdf5 file and the matching line in an XDMF file.
 *
 * @param[in] array_dims dimensions of the array to write
 * @param[in] array      array to write
 * @param[in] array_name name of the array
 * @param[in] xdmf_file  pointer to the xml file
 * @param[in] center     values are centered on ("Node" or "Cell")
 */
void xdmf_array_2d(char* mesh_name, hsize_t array_dims[2], const void* array,
        char* array_name, FILE* xdmf_file, char* center) {
    char file_name[111];
    
    sprintf(file_name, "%s-%s.h5", mesh_name, array_name);
    hid_t hfile_id = hdf5_ser_file_create(file_name);
    sprintf(file_name, "/%s", array_name);
    //printf("array_dims=%d %d\n",(int)array_dims[0],(int)array_dims[1]);
    hdf5_ser_write_dble_array_2d(hfile_id, array_dims, array, file_name);
    hdf5_ser_file_close(hfile_id);
    sprintf(file_name, "%s-%s.h5:/%s", mesh_name, array_name, array_name);
    xml_field_2d(xdmf_file, array_name, file_name, (int)array_dims[0], (int)array_dims[1], center);
}

/*
 * Close the XML file and finish to write last lines.
 *
 * @param[in, out] xdmf_file, the FILE* or your xml file
 */
void xdmf_close(FILE** xdmf_file) {
    fprintf(*xdmf_file, "</Grid>\n");
    fprintf(*xdmf_file, "</Domain>\n");
    fprintf(*xdmf_file, "</Xdmf>\n");
    fclose(*xdmf_file);
}

/*
 * Add the the good value of time in VisIt plot.
 *
 * @param[in] xdmf_file, the FILE* or your xml file
 * @param[in] time, the input time
 */
void xdmf_set_time(FILE* xdmf_file, double time) {
    fprintf(xdmf_file, "<Time Value='%.*g'/>\n", DBL_DECIMAL_DIG, time);
}

// When calling with 4 arguments, the macro will set time to -1 and array_name to "f".
// There is a ambiguity with 5 arguments. Either time or array_name is missing. This
// is checked by the type of the 5th argument (that should be either double either char*).
// WARNING: This check needs C11 support from your compiler (_Generic). If your compiler
// does not support C11, you can use instead the P99 macros from Jens Gustedt:
// https://gustedt.wordpress.com/2012/01/02/emulating-c11-compiler-features-with-gcc-_generic/
// But the easiest choice is maybe to fill all the arguments of the function when you call it.
#define plot_f_cartesian_mesh_2d_4(a, b, c, d      ) a, b, c, d, -1, "f"
#if defined(SELA_C11)
#define plot_f_cartesian_mesh_2d_5(a, b, c, d, e   ) a, b, c, d, \
  _Generic((e), char*: -1, default: e), \
  _Generic((e), char*:  e, default: "f")
#endif
#define plot_f_cartesian_mesh_2d_6(a, b, c, d, e, f) a, b, c, d, e, f
#define plot_f_cartesian_mesh_2d_7(a, b, c, d, e, f, g) a, b, c, d, e, f, g
#define plot_f_cartesian_mesh_2d(...) VARIADIC(plot_f_cartesian_mesh_2d, NUMARG32(__VA_ARGS__), __VA_ARGS__)

/*
 * Plot 2d distribution function for VisIt.
 * This routine will create a file named [array_name][iplot].xmf
 * Note: also calls the VARIADIC macro, but you wouldn't notice unless this portion of code doesn't compile :)
 *
 * @param[in] iplot   plot counter.
 * @param[in] f       function values.
 * @param[in] mesh1, mesh2 the two 1d meshes.
 * @param[in] time the plot time. Optional, defaults to -1.
 * @param[in] array_name a name for the array. Optional, defaults to "f".
 * @param[in] folder : name for the folder
 */
void plot_f_cartesian_mesh_2d(int iplot, const void* f, mesh_1d mesh1, mesh_1d mesh2,
        double time, char* array_name, char* folder) {
    hid_t hfile_id;
    char cplot[6]; // 4 digits + '\0' [+1 more to avoid a warning in sprintf]
    char mesh_name[111];
    char file_name[117];
    int nnodes_x1 = mesh1.size;
    int nnodes_x2 = mesh2.size;
    hsize_t array_dims[2] = { nnodes_x1, nnodes_x2 };
    
    if (iplot < 10)
        sprintf(cplot, "000%d", iplot & 0xf); // 0xf = 15
    else if (iplot < 100)
        sprintf(cplot, "00%d", iplot & 0x7f); // 0x7f = 127
    else if (iplot < 1000)
        sprintf(cplot, "0%d", iplot & 0x3ff); // 0x3ff = 1023
    else
        sprintf(cplot, "%d", iplot & 0x3fff); // 0x3fff = 16383
    cplot[4] = 0;
    sprintf(mesh_name, "cartesian_mesh-%s", array_name);
    
    if (iplot == 1) {
		double (*x1)[nnodes_x2] = malloc((nnodes_x1) * sizeof *x1); 
		double (*x2)[nnodes_x2] = malloc((nnodes_x1) * sizeof *x2); 
        //double x1[nnodes_x1][nnodes_x2];
        //double x2[nnodes_x1][nnodes_x2];
        for (int i = 0; i < nnodes_x1; i++) {
            for (int j = 0; j < nnodes_x2; j++) {
                x1[i][j] = mesh1.array[i];
                x2[i][j] = mesh2.array[j];
            }
        }
        // Creation of the hdf5 file for the x1 mesh.
        sprintf(file_name, "%s%s-x1.h5", folder, mesh_name);
        hfile_id = hdf5_ser_file_create(file_name);
        hdf5_ser_write_dble_array_2d(hfile_id, array_dims, x1, "/x1");
        hdf5_ser_file_close(hfile_id);
        // Creation of the hdf5 file for the x2 mesh.
        sprintf(file_name, "%s%s-x2.h5", folder, mesh_name);
        hfile_id = hdf5_ser_file_create(file_name);
        hdf5_ser_write_dble_array_2d(hfile_id, array_dims, x2, "/x2");
        hdf5_ser_file_close(hfile_id);
    }
    // Creation of the hdf5 file for the f array.
    sprintf(file_name, "%s%s%s.xmf", folder, array_name, cplot);
    FILE* xdmf_file = xdmf_open_2d(file_name, mesh_name, nnodes_x1, nnodes_x2);
    sprintf(file_name, "%s%s%s", folder, array_name, cplot);
    xdmf_set_time(xdmf_file, time);
    xdmf_array_2d(file_name, array_dims, f, "values", xdmf_file, "Node");
    xdmf_close(&xdmf_file);
}

#endif // ifndef SELA_HDF5_IO
