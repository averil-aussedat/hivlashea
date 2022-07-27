#ifndef SELA_PARAMETER_READER
#define SELA_PARAMETER_READER

/*
 * Emulate a subset of the yaml language format, as used with CEA's library
 * paraconf (https://gitlab.maisondelasimulation.fr/jbigot/libparaconf)
 *    > The indentation must be 2 white spaces.
 *    > We only handle mappings (a set of key: value) possibilities, where:
 *          - a key must be a string
 *          - a value must be a string or a mapping
 */

#include <ctype.h>          // function  isspace
#include <errno.h>          // variable  errno (e.g. set by fgets)
                            // constant  ERANGE
#include <limits.h>         // constant  INT_MAX
#include <stdbool.h>        // type      bool
                            // constants true, false
#include <stdio.h>          // functions printf, fgets, fopen
#include <stdlib.h>         // functions strtol (converts string to long), malloc, free ((de)allocate memory)
#include <string.h>         // functions memcpy, strcmp, strlen, strstr
#include "string_helpers.h" // macro     ERROR_MESSAGE
                            // functions trim_spaces, split_in_two

// Helpers to read parameters from file
#ifndef SELA_MAX_DEPTH
#    define SELA_MAX_DEPTH 256
#endif
#ifndef SELA_MAX_LINE_LENGTH
#    define SELA_MAX_LINE_LENGTH 2048
#endif
#define YAML_DELIM ":"

/*
 * Remove yaml comments from a line.
 *
 * See https://yaml.org/spec/1.2/spec.html#id2780069:
 * "Comments begin with the number sign (#), can start anywhere on a line and
 *  continue until the end of the line. Comments must be separated from other
 *  tokens by whitespace characters."
 */
void remove_comments(char* str) {
    if (*str == '#') {
        *str = '\0';
        return;
    }
    while (*str != '\0') {
        if (isspace((unsigned char)*str) && *(str+1) == '#')
            break;
        str++;
    }
    *str = '\0';
}

/*
 * A node in a yaml document. Here, we only parse mappings and scalars.
 *
 * See https://yaml.org/spec/1.2/spec.html#id2764044
 * "A YAML node represents a single native data structure. Such nodes have
 *  content of one of three kinds: scalar, sequence, or mapping. In addition,
 *  each node has a tag which serves to restrict the set of possible values the
 *  content can have.
 *  > Scalar: The content of a scalar node is an opaque datum that can be
 *    presented as a series of zero or more Unicode characters. 
 *  > Sequence: The content of a sequence node is an ordered series of zero or
 *    more nodes. In particular, a sequence may contain the same node more than
 *    once. It could even contain itself (directly or indirectly). 
 *  > Mapping: The content of a mapping node is an unordered set of key: value
 *    node pairs, with the restriction that each of the keys is unique. YAML
 *    places no further restrictions on the nodes. In particular, keys may be
 *    arbitrary nodes, the same node may be used as the value of several key:
 *    value pairs, and a mapping could even contain itself as a key or a value
 *    (directly or indirectly)."
 */
typedef struct sela_node {
    struct sela_node* first_child;
    struct sela_node* next_sibling;
    char* key;
    char* value;
    long long_value;
    double double_value;
} sela_node;
typedef struct sela_node* PC_tree_t;

/*
 * Read a long from a string.
 *
 * @param[in]  string the string that should contain an int to read.
 * @param[out] long_value the value of the long read (0 if nothing could be read).
 */
bool long_from_string(char* string, long int* long_value) {
    errno = 0;
    char* strtol_after_read;
    *long_value = strtol(string, &strtol_after_read, 10);
    return errno != ERANGE && string != strtol_after_read && *long_value >= 0;
}

/*
 * Read an int from a string.
 *
 * @param[in]  string the string that should contain an int to read.
 * @param[out] int_value the value of the int read (0 if nothing could be read).
 */
bool int_from_string(char* string, int* int_value) {
    long int long_value;
    if (long_from_string(string, &long_value) && long_value <= INT_MAX) {
        *int_value = (int)long_value;
        return true;
    }
    *int_value = 0;
    return false;
}

/*
 * Read a double from a string.
 *
 * @param[in]  string the string that should contain a double to read.
 * @param[out] double_value the value of the double read (0. if nothing could be read).
 */
bool double_from_string(char* string, double* double_value) {
    errno = 0;
    char* strtod_after_read;
    *double_value = strtod(string, &strtod_after_read);
    return errno != ERANGE && string != strtod_after_read;
}

/*
 * Auxiliary function to print a sela_node on the standard output. It prints
 * the indentation of the current node, the key, value, and the tree nodes
 * that follow it.
 *
 * @param[in] parameters a pointer to the sela_node.
 * @param[in] tree_depth the depth of the current sela_node in the full tree to be printed.
 */
void print_yaml_from_parameters_aux(PC_tree_t parameters, int tree_depth) {
    if (parameters) {
        if (tree_depth > 0) {
            // The root node is only part of the representation, it is not present in the yaml file.
            // Print the indentation.
            for (int i = 1; i < tree_depth; i++) {
                printf("  ");
            }
            // Print the key, value (as a string).
            printf("%s: %s\n", parameters->key, parameters->value);
        }
        // Print the children.
        print_yaml_from_parameters_aux(parameters->first_child,  tree_depth + 1);
        if (tree_depth > 0) {
            // Print the siblings.
            // Sometimes we can invoke print_yaml_from_parameters on a node which
            // is not the root node. In that case, we do not want to print the
            // siblings of this node. This test does not change anything for the
            // root node which does not have siblings anyway.
            print_yaml_from_parameters_aux(parameters->next_sibling, tree_depth);
        }
    }
}

/*
 * Prints the current sela_node on the standard output. It prints the same
 * content contained in the original yaml file --- excluding comments ---
 * if we are printing the root node.
 *
 * @param[in] parameters a pointer to the sela_node.
 */
void print_yaml_from_parameters(PC_tree_t parameters) {
    print_yaml_from_parameters_aux(parameters, 0);
}

/*
 * Reads a parameter file, written in a subset of the yaml language format
 * (see beginning of this file).
 *
 * @param filename the name of the file containing the yaml description of the parameters.
 */
PC_tree_t PC_parse_path(const char* filename) {
    char line[SELA_MAX_LINE_LENGTH + 1]; // Each line read
    int count = 0;                       // Number of the last line read
    char* parameter_key;                 // Key of the current yaml entry
    char* parameter_value;               // Value of the current yaml entry
    int tree_depth = 1;
    PC_tree_t last_node_per_level[SELA_MAX_DEPTH];
    // Root node
    last_node_per_level[0] = malloc(sizeof(sela_node));
    last_node_per_level[0]->first_child = (void*)0;
    last_node_per_level[0]->next_sibling = (void*)0;
    last_node_per_level[0]->key = (char*)0;
    last_node_per_level[0]->value = (char*)0;
    last_node_per_level[0]->long_value = 0;
    last_node_per_level[0]->double_value = 0.;
    // Last node seen are set to (void*)0 for all other levels.
    for (int i = 1; i < SELA_MAX_DEPTH; i++) {
        last_node_per_level[i] = (void*)0;
    }
    
    // Read the parameter file.
    FILE* fp = fopen(filename, "r");
    if (!fp) { // Error in file opening
        ERROR_MESSAGE("I could not open %s.\n", filename);
    }
    while (fgets(line, SELA_MAX_LINE_LENGTH, fp)) { // Successful line reading
        count++;
        line[strlen(line) - 1] = '\0'; // Useful when we need to print the line
        remove_comments(line);         // Allows the user to put comments in the parameter file
        char* trimmed_line = malloc((SELA_MAX_LINE_LENGTH + 1) * sizeof(char));
        memcpy(trimmed_line, line, SELA_MAX_LINE_LENGTH + 1);
        bool is_empty_line = !trim_spaces(trimmed_line)[0];
        free(trimmed_line);
        if (is_empty_line) // Empty line
            continue;
        int previous_tree_depth = tree_depth;
        tree_depth = 1;
        for (int i = 0; i < previous_tree_depth + 1; i++) {
            if (line[2 * i] == ' ') {
                if (line[2 * i + 1] == ' ') {
                    tree_depth++;
                } else {
                    ERROR_MESSAGE("Indentation problem (odd number of leading spaces) on line %d (%s).\n", count, line);
                }
            } else {
                break;
            }
        }
        split_in_two(line, YAML_DELIM, &parameter_key, &parameter_value);
        if (!parameter_value) {// YAML_DELIM was not found in the line
            ERROR_MESSAGE("Missing %s on line %d (%s).\n", YAML_DELIM, count, line);
        }
        
        // Creation of the tree node.
        PC_tree_t current_node = malloc(sizeof(sela_node));
        current_node->first_child = (void*)0;
        current_node->next_sibling = (void*)0;
        int key_length = strlen(parameter_key) + 1;
        current_node->key = malloc(key_length * sizeof(char));
        memcpy(current_node->key, parameter_key, key_length);
        int value_length = strlen(parameter_value) + 1;
        current_node->value = malloc(value_length * sizeof(char));
        memcpy(current_node->value, parameter_value, value_length);
        long_from_string(parameter_value, &(current_node->long_value));
        double_from_string(parameter_value, &(current_node->double_value));
        
        // If we go back in depth, it is a new branch, reset the last nodes seen.
        if (tree_depth < previous_tree_depth) {
            for (int i = tree_depth + 1; i < SELA_MAX_DEPTH; i++) {
                last_node_per_level[i] = (void*)0;
            }
        }
        // Find our previous sibling (if it exists)
        PC_tree_t previous_sibling = last_node_per_level[tree_depth];
        if (previous_sibling) {
            // If we have a sibling, bind ourselves to this sibling
            previous_sibling->next_sibling = current_node;
        } else {
            // else, find our direct ancestor (it should exist if the file is correctly written)
            PC_tree_t ancestor = last_node_per_level[tree_depth - 1];
            if (ancestor) {
                // if we have a direct ancestor, bind ourselves to this direct ancestor
                ancestor->first_child = current_node;
            } else {
                // if we do not have a direct ancestor, there is a problem in the yaml file
                ERROR_MESSAGE("Indentation problem (too much leading spaces) on line %d (%s).\n", count, line);
            }
        }
        last_node_per_level[tree_depth] = current_node;
    }
/* DEBUG */
//    print_yaml_from_parameters(last_node_per_level[0]);
/* END DEBUG */
    return last_node_per_level[0];
}

/*
 * Outputs the sela_node whose key is str among the siblings
 * of the current sela_node, or (void*)0 if not found.
 *
 * @param[in] parameters a pointer to the sela_node.
 * @param[in] str the string to search for.
 * @return    the first sibling to have str as key (if found) or
 *            (void*)0 if not found.
 */
PC_tree_t sela_find(PC_tree_t parameters, char* str) {
    if (!parameters)
        return (void*)0;
    if (strcmp(parameters->key, str) == 0) {
        return parameters;
    }
    return sela_find(parameters->next_sibling, str);
}

/*
 * Outputs the sela_node whose key is '.' + str among the children
 * of the current sela_node, or (void*)0 if not found.
 *
 * @param[in] parameters a pointer to the sela_node.
 * @param[in] str the string to search for.
 * @return    the first child to have str as key (if found) or
 *            (void*)0 if not found.
 */
PC_tree_t PC_get(PC_tree_t parameters, char* str) {
    if (!str || strlen(str) == 0 || str[0] != '.' || !parameters)
        return (void*)0;
    return sela_find(parameters->first_child, str+1);
}

/*
 * Outputs the string value contained in the current sela_node.
 *
 * @param[in]  parameters a pointer to the sela_node.
 * @param[out] out a pointer to the output string.
 */
void PC_string(PC_tree_t parameters, char** out) {
    *out = parameters->value;
}

/*
 * Outputs the double value contained in the current sela_node.
 *
 * @param[in]  parameters a pointer to the sela_node.
 * @param[out] out a pointer to the output double.
 */
void PC_double(PC_tree_t parameters, double* out) {
    *out = parameters->double_value;
}

/*
 * Outputs the integer value contained in the current sela_node.
 *
 * @param[in]  parameters a pointer to the sela_node.
 * @param[out] out a pointer to the output integer.
 */
void PC_int(PC_tree_t parameters, long* out) {
    *out = parameters->long_value;
}

/*
 * Free the memory allocated for the current sela_node.
 *
 * @param[in] parameters_p a pointer to the pointer on the sela_node.
 */
void PC_tree_destroy(PC_tree_t* parameters_p) {
    if (*parameters_p) {
        PC_tree_destroy(&(*parameters_p)->first_child);
        PC_tree_destroy(&(*parameters_p)->next_sibling);
        free((*parameters_p)->value);
        free((*parameters_p)->key);
        free(parameters_p);
    }
}

#endif // ifndef SELA_PARAMETER_READER
