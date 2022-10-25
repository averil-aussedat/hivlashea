#ifndef SELA_VP_1D1V_CART_SPLIT
#define SELA_VP_1D1V_CART_SPLIT

#include <stdbool.h>          // type      bool
                              // constants true, false
#include <string.h>           // function  strcmp
#include "parameter_reader.h" // type      PC_tree_t
                              // functions PC_get, PC_double, PC_int, PC_string
#include "string_helpers.h"   // macro     ERROR_MESSAGE

/*****************************************************************************
 *                              Time splitting                               *
 *                                                                           *
 * selalib/src/time_integration/splitting_methods/                           *
 *                                            sll_m_time_splitting_coeff.F90 *
 *****************************************************************************/

#define sll_p_lie_tv             -2
#define sll_p_lie_vt             -1
#define sll_p_strang_tvt         0
#define sll_p_strang_vtv         1
#define sll_p_triple_jump_tvt    2
#define sll_p_triple_jump_vtv    3
#define sll_p_order6_vtv         4
#define sll_p_order6_tvt         5

typedef struct sll_t_splitting_coeff sll_t_splitting_coeff;
struct sll_t_splitting_coeff {
    int     split_case;
    double* split_step;
    int     nb_split_step;
    bool    split_begin_T;
    double dt;
    int num_iteration;
    int plot_frequency;
};

sll_t_splitting_coeff sll_f_new_time_splitting_coeff(  
        int split_case, double dt, int num_iteration, int plot_frequency) {
    sll_t_splitting_coeff split;
    split.dt = dt;
    split.num_iteration = num_iteration;
    split.plot_frequency = plot_frequency;
    split.split_case = split_case;
    switch (split_case) {
        case sll_p_lie_tv:
            split.nb_split_step = 2;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = true;
            split.split_step[0] = 1.;
            split.split_step[1] = 1.;
            break;
        case (sll_p_lie_vt):
            split.nb_split_step = 2;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = false;
            split.split_step[0] = 1.;
            split.split_step[1] = 1.;
            break;
        case (sll_p_strang_tvt): // Strang splitting TVT
            split.nb_split_step = 3;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = true;
            split.split_step[0] = 0.5;
            split.split_step[1] = 1.;
            split.split_step[2] = split.split_step[0];
            break;
        case (sll_p_strang_vtv): // Strang splitting VTV
            split.nb_split_step = 3;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = false;
            split.split_step[0] = 0.5;
            split.split_step[1] = 1.;
            split.split_step[2] = split.split_step[0];
            break;
        case (sll_p_triple_jump_tvt): // triple jump TVT
            split.nb_split_step = 7;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = true;
            split.split_step[0] =  0.675603595979829;
            split.split_step[1] =  1.351207191959658;
            split.split_step[2] = -0.17560359597982855;
            split.split_step[3] = -1.702414383919315;
            split.split_step[4] = split.split_step[2];
            split.split_step[5] = split.split_step[1];
            split.split_step[6] = split.split_step[0];
            break;
        case (sll_p_triple_jump_vtv): // triple jump VTV
            split.nb_split_step = 7;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = false;
            split.split_step[0] =  0.675603595979829;
            split.split_step[1] =  1.351207191959658;
            split.split_step[2] = -0.17560359597982855;
            split.split_step[3] = -1.702414383919315;
            split.split_step[4] = split.split_step[2];
            split.split_step[5] = split.split_step[1];
            split.split_step[6] = split.split_step[0];
            break;
        case (sll_p_order6_vtv): // Order 6 VTV (O6-11 of Blanes)
            split.nb_split_step = 23;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = false;
            split.split_step[ 0] =  0.0414649985182624;
            split.split_step[ 1] =  0.123229775946271;
            split.split_step[ 2] =  0.198128671918067;
            split.split_step[ 3] =  0.290553797799558;
            split.split_step[ 4] = -0.0400061921041533;
            split.split_step[ 5] = -0.127049212625417;
            split.split_step[ 6] =  0.0752539843015807;
            split.split_step[ 7] = -0.246331761062075;
            split.split_step[ 8] = -0.0115113874206879;
            split.split_step[ 9] =  0.357208872795928;
            split.split_step[10] =  0.23666992478693111;
            split.split_step[11] =  0.20477705429147008;
            split.split_step[12] = split.split_step[10];
            split.split_step[13] = split.split_step[ 9];
            split.split_step[14] = split.split_step[ 8];
            split.split_step[15] = split.split_step[ 7];
            split.split_step[16] = split.split_step[ 6];
            split.split_step[17] = split.split_step[ 5];
            split.split_step[18] = split.split_step[ 4];
            split.split_step[19] = split.split_step[ 3];
            split.split_step[20] = split.split_step[ 2];
            split.split_step[21] = split.split_step[ 1];
            split.split_step[22] = split.split_step[ 0];
            break;
        case (sll_p_order6_tvt): // Order 6 TVT (O6-14 of Blanes)
            split.nb_split_step = 29;
            split.split_step = malloc(split.nb_split_step * sizeof(double));
            split.split_begin_T = true;
            split.split_step[ 0] =  0.0378593198406116;
            split.split_step[ 1] =  0.09171915262446165;
            split.split_step[ 2] =  0.102635633102435;
            split.split_step[ 3] =  0.183983170005006;
            split.split_step[ 4] = -0.0258678882665587;
            split.split_step[ 5] = -0.05653436583288827;
            split.split_step[ 6] =  0.314241403071447;
            split.split_step[ 7] =  0.004914688774712854;
            split.split_step[ 8] = -0.130144459517415;
            split.split_step[ 9] =  0.143761127168358;
            split.split_step[10] =  0.106417700369543;
            split.split_step[11] =  0.328567693746804;
            split.split_step[12] = -0.00879424312851058;
            split.split_step[13] = -0.196411466486454234;
            split.split_step[14] =  0.20730506905689536;
            split.split_step[15] = split.split_step[13];
            split.split_step[16] = split.split_step[12];
            split.split_step[17] = split.split_step[11];
            split.split_step[18] = split.split_step[10];
            split.split_step[19] = split.split_step[ 9];
            split.split_step[20] = split.split_step[ 8];
            split.split_step[21] = split.split_step[ 7];
            split.split_step[22] = split.split_step[ 6];
            split.split_step[23] = split.split_step[ 5];
            split.split_step[24] = split.split_step[ 4];
            split.split_step[25] = split.split_step[ 3];
            split.split_step[26] = split.split_step[ 2];
            split.split_step[27] = split.split_step[ 1];
            split.split_step[28] = split.split_step[ 0];
            break;
        default:
            ERROR_MESSAGE("Bad time-splitting case.\n");
    }
    return split;
}

void free_splitting_coeff(sll_t_splitting_coeff* split) {
    free(split->split_step);
}

void splitting(PC_tree_t conf, sll_t_splitting_coeff* split) {
    int splitcase;
    double dt;
    long num_iteration;
    int plot_frequency;
    char* a_string = "STRANG_XVX"; // Default case.
    if (PC_get(conf, ".case")) {
        PC_string(PC_get(conf, ".case"), &a_string);
    }
    if (PC_get(conf, ".dt")) {
        PC_double(PC_get(conf, ".dt"), &dt);
    } else {
        ERROR_MESSAGE("#Error in splitting %s: dt parameter missing.\n", conf->key);
    }
    if (PC_get(conf, ".num_iteration")) {
        PC_int(PC_get(conf, ".num_iteration"), &num_iteration);
    } else {
        ERROR_MESSAGE("#Error in splitting %s: num_iteration parameter missing.\n", conf->key);
    }
    if (PC_get(conf, ".plot_frequency")) {
        PC_int(PC_get(conf, ".plot_frequency"), &plot_frequency);
    } else {
        ERROR_MESSAGE("#Error in splitting %s: plot_frequency parameter missing.\n", conf->key);
    }
    if (strcmp(a_string, "LIE_XV") == 0) {
        splitcase = sll_p_lie_tv;
    } else if (strcmp(a_string, "LIE_VX") == 0) {
        splitcase = sll_p_lie_vt;
    } else if (strcmp(a_string, "STRANG_XVX") == 0) {
        splitcase = sll_p_strang_tvt;
    } else if (strcmp(a_string, "STRANG_VXV") == 0) { // Strang splitting VTV
        splitcase = sll_p_strang_vtv;
    } else if (strcmp(a_string, "TRIPLE_JUMP_XVX") == 0) { // triple jump TVT
        splitcase = sll_p_triple_jump_tvt;
    } else if (strcmp(a_string, "TRIPLE_JUMP_VXV") == 0) { // triple jump VTV
        splitcase = sll_p_triple_jump_vtv;
    } else if (strcmp(a_string, "ORDER6_VXV") == 0) { // Order 6 VTV (O6-11 of Blanes)
        splitcase = sll_p_order6_vtv;
    } else if (strcmp(a_string, "ORDER6_XVX") == 0) { // Order 6 TVT (O6-14 of Blanes)
        splitcase = sll_p_order6_tvt;
    } else {
        ERROR_MESSAGE("#Error in splitting %s: unknown case.\n#Possible cases are 'LIE_XV', 'LIE_VX', 'STRANG_XVX', 'STRANG_XVX', 'TRIPLE_JUMP_XVX', 'TRIPLE_JUMP_VXV', 'ORDER6_VXV' or 'ORDER6_XVX'.\n", conf->key);
    }
    
    *split = sll_f_new_time_splitting_coeff(splitcase, dt, num_iteration, plot_frequency);
}

#endif // ifndef SELA_VP_1D1V_CART_SPLIT
