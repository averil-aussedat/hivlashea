#ifndef SELA_VP_1D1V_CART_INIT_ANY_EXPR
#define SELA_VP_1D1V_CART_INIT_ANY_EXPR

#include <math.h>             // function  cos, exp, sqrt
                              // constant  M_PI
#include <stdbool.h>          // type      bool
                              // constants true, false
#include <stdlib.h>           // functions malloc, free ((de)allocate memory)
                              // type      size_t
#include <string.h>           // functions strcmp, strlen
#include "parameter_reader.h" // type      PC_tree_t
                              // functions PC_get, PC_double, PC_string
#include "string_helpers.h"   // functions split_in_two, count_chars
                              // macro     ERROR_MESSAGE
#include "remap.h"            // type      parallel_stuff
#include "tinyexpr.h"         // types     te_expr, te_variable
                              // functions te_compile, te_eval, te_free
                              // constants TE_VARIABLE, TE_CLOSURE0, functions
#include "variadic.h"         // macros    VARIADIC, NUMARG32

/*
 * A user may specify any expression which follows these conditions:
 *     - the accepted mathematical operators are
 *           "+": addition
 *           "-": substraction
 *           "*": multiplication
 *           "/": division
 *           "^": exponentiation
 *           "%": modulo
 *     - the accepted mathematical constants are
 *           "e":  2.71828182845904523536
 *           "pi": 3.14159265358979323846
 *     - the accepted mathematical functions are
 *           "abs":   x -> absolute value of x
 *           "acos":  x -> arc cosine of x
 *           "asin":  x -> arc sine of x
 *           "atan":  x -> arc tangent of x
 *           "atan2": (x, y) -> arc tangent of y/x (uses x and y signs to determine the quadrant)
 *           "ceil":  x -> smallest integer value >= x
 *           "cos":   x -> cosine of x
 *           "cosh":  x -> hyperbolic cosine of x
 *           "exp":   x -> e^x
 *           "fac":   x -> x!
 *           "floor": x -> largest integer value <= x
 *           "ln":    x -> natural logarithm of x (base e)
 *           "log":   x -> natural logarithm of x (base e) --- because in math.h, log = ln
 *           "log10": x -> logarithm of x (base 10)
 *           "ncr":   (n, r) -> n choose r --- number of combinations of r objects among n
 *           "npr":   (n, r) -> n permutations r --- number of permutations of r objects among n
 *           "pow":   (a, b) -> a^b
 *           "sin":   x -> sine of x
 *           "sinh":  x -> hyperbolic sine of x
 *           "sqrt":  x -> square root of x
 *           "tan":   x -> tangent of x
 *           "tanh":  x -> hyperbolic tangent of x
 *     - the accepted logical constants are
 *           "true":  1.
 *           "false": 0.
 *     - the accepted logical operators are
 *           ">":  greater
 *           ">=": greater or equal
 *           "<":  lower
 *           "<=": lower or equal
 *           "==": equal
 *           "&&": and
 *           "||": or
 *           "!":  not
 *     - the accepted logical function is
 *           "if": (test, a, b) -> if (test) then a else b
 *     - the accepted variables and/or parameters are any string of the form [a-z][a-z0-9_]*
 *       not already listed before as a mathematical constant or a mathematical function.
 */

// When calling with 2 or 3 arguments, the macro will set nb_reserved_names to 0 and reserved_names to (char**)0.
#define check_valid_variable_name_2(a, b      ) a, b, 0, (char**)0
#define check_valid_variable_name_3(a, b, c   ) a, b, 0, (char**)0
#define check_valid_variable_name_4(a, b, c, d) a, b, c, d
#define check_valid_variable_name(...) VARIADIC(check_valid_variable_name, NUMARG32(__VA_ARGS__), __VA_ARGS__)

/*
 * Checks that a parameter name given by the user is valid with respect to the current function.
 * Note: also calls the VARIADIC macro, but you wouldn't notice unless this portion of code doesn't compile :)
 *
 * @param[in] name the parameter name to be checked.
 * @param[in] fun_name the name of the function in which it appears.
 * @param[in] nb_reserved_names number of other reserved names. Optional, defaults to 0.
 * @param[in] reserved_names other reserved names. Optional, defaults to (char**)0.
 */
void check_valid_variable_name(char* name, char* fun_name, int nb_reserved_names, char** reserved_names) {
    // non-empty parameter name
    if (!name || name[0] == '\0') {
        ERROR_MESSAGE("#Error in function %s: a parameter name cannot be empty.\n", fun_name);
    }
    // tinyexpr only accepts names of the form [a-z][a-z0-9_]*
    if (name[0] < 'a' || name[0] > 'z') {
        ERROR_MESSAGE("#Error in function %s: a parameter name must begin by 'a-z' and contain only 'a-z0-9_'.\n", fun_name);
    }
    for (int i = 1; i < strlen(name); i++) {
        if (!((name[i] >= 'a' && name[i] <= 'z') || (name[i] >= '0' && name[i] <= '9') || (name[i] == '_'))) {
            ERROR_MESSAGE("#Error in function %s: a parameter name must begin by 'a-z' and contain only 'a-z0-9_'.\n", fun_name);
        }
    }
    // tinyexpr has reserved function names
    size_t nb_functions = sizeof(functions) / sizeof(functions[0]) - 1; // Last member of functions is { 0, 0, 0, 0 } and is not counted.
    for (int i = 0; i < nb_functions; i++) {
        if (!strcmp(name, functions[i].name)) {
            ERROR_MESSAGE("#Error in function %s: '%s' cannot be a parameter name, it is reserved.\n", fun_name, name);
        }
    }
    // the user has already used some parameter names
    for (int i = 0; i < nb_reserved_names; i++) {
        if (!strcmp(name, reserved_names[i])) {
            ERROR_MESSAGE("#Error in function %s: '%s' cannot be used multiple times as a parameter name.\n", fun_name, name);
        }
    }
}

#define DEFINITION_DELIM "->"
#define VARIABLES_DELIM ","

typedef struct tiny_expr {
    // Arity of the function, followed by the arrays of variable names and values.
    int arity;
    char** variable_names;
    double* variable_values;
    // The te_expr uses those values for evaluation, and also the values of previous parameters.
    te_variable* vars;
    te_expr* expr;
} tiny_expr;

double cloture1(void *context, double a) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 1) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}
double cloture2(void *context, double a, double b) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 2) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        tiny_expression->variable_values[1] = b;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}
double cloture3(void *context, double a, double b, double c) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 3) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        tiny_expression->variable_values[1] = b;
        tiny_expression->variable_values[2] = c;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}
double cloture4(void *context, double a, double b, double c, double d) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 4) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        tiny_expression->variable_values[1] = b;
        tiny_expression->variable_values[2] = c;
        tiny_expression->variable_values[3] = d;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}
double cloture5(void *context, double a, double b, double c, double d, double e) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 5) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        tiny_expression->variable_values[1] = b;
        tiny_expression->variable_values[2] = c;
        tiny_expression->variable_values[3] = d;
        tiny_expression->variable_values[4] = e;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}
double cloture6(void *context, double a, double b, double c, double d, double e, double f) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 6) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        tiny_expression->variable_values[1] = b;
        tiny_expression->variable_values[2] = c;
        tiny_expression->variable_values[3] = d;
        tiny_expression->variable_values[4] = e;
        tiny_expression->variable_values[5] = f;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}
double cloture7(void *context, double a, double b, double c, double d, double e, double f, double g) {
    if (context) {
        tiny_expr* tiny_expression = (tiny_expr*)context;
        if (tiny_expression->arity != 7) {
            ERROR_MESSAGE("#Evaluating a cloture with the wrong arity.\n");
        }
        tiny_expression->variable_values[0] = a;
        tiny_expression->variable_values[1] = b;
        tiny_expression->variable_values[2] = c;
        tiny_expression->variable_values[3] = d;
        tiny_expression->variable_values[4] = e;
        tiny_expression->variable_values[5] = f;
        tiny_expression->variable_values[6] = g;
        return te_eval(tiny_expression->expr);
    }
    ERROR_MESSAGE("#Evaluating a cloture with no context.\n");
}

void any_fun_from_string(char* definition, char* fun_name,
        int nb_previous_params, char** previous_params_names,
        double* previous_params_values, int* previous_params_arities, tiny_expr* previous_params_functions,
        tiny_expr* tiny_expression) {
    // Getting the variable string.
    char* variables_string = (char*)0;
    char* expression       = (char*)0;
    split_in_two(definition, DEFINITION_DELIM, &variables_string, &expression);
    if (!expression) {
        // There is no function arrow, this is a constant --- maybe depending on previous parameters.
        tiny_expression->arity = 0;
        expression = variables_string;
        variables_string = (char*)0;
    } else {
        // There is a function arrow, parse the variables.
        tiny_expression->arity = count_chars(variables_string, ',') + 1;
        if (tiny_expression->arity > 7) {
            ERROR_MESSAGE("#Error in function %s: functions of arity > 7 are not handled.\n", fun_name);
        }
        int length = strlen(variables_string);
        int nb_opening_parentheses = count_chars(variables_string, '(');
        int nb_closing_parentheses = count_chars(variables_string, ')');
        if (nb_opening_parentheses > 1 || nb_closing_parentheses > 1) {
            ERROR_MESSAGE("#Error in function %s: one pair of parentheses at most is allowed to define the variables.\n", fun_name);
        } else if (nb_opening_parentheses != nb_closing_parentheses) {
            ERROR_MESSAGE("#Error in function %s: parentheses mismatch in the definition of the variables.\n", fun_name);
        } else if (nb_opening_parentheses == 1) {
            if (variables_string[0] != '(' || variables_string[length-1] != ')') {
                ERROR_MESSAGE("#Error in function %s: the parentheses can only be used to surround the variable names.\n", fun_name);
            }
            variables_string[length-1] = '\0';
            variables_string++;
        }
    }
    // Converting the string to names.
    tiny_expression->variable_names  = malloc((tiny_expression->arity) * sizeof(char*));
    tiny_expression->variable_values = malloc((tiny_expression->arity) * sizeof(double));
    for (int i = 1; i < tiny_expression->arity; i++) {
        split_in_two(variables_string, VARIABLES_DELIM, &tiny_expression->variable_names[i - 1], &variables_string);
        check_valid_variable_name(tiny_expression->variable_names[i - 1], fun_name, nb_previous_params, previous_params_names);
    }
    if (tiny_expression->arity > 0) {
        tiny_expression->variable_names[tiny_expression->arity - 1] = variables_string;
        check_valid_variable_name(tiny_expression->variable_names[tiny_expression->arity - 1], fun_name, nb_previous_params, previous_params_names);
    }
    // The parameters names and values are given as input, we build the total number of bound
    // names from the variable names and previous parameter names.
    int nb_parameters = tiny_expression->arity + nb_previous_params;
    
    // Building the tinyexpression
    tiny_expression->vars = malloc((nb_parameters) * sizeof(te_variable));
    for (int i = 0; i < nb_previous_params; i++) {
        tiny_expression->vars[i].name    = previous_params_names[i];
        if (previous_params_arities[i] == 0)
            tiny_expression->vars[i].address = &previous_params_values[i];
        else if (previous_params_arities[i] == 1)
            tiny_expression->vars[i].address = cloture1;
        else if (previous_params_arities[i] == 2)
            tiny_expression->vars[i].address = cloture2;
        else if (previous_params_arities[i] == 3)
            tiny_expression->vars[i].address = cloture3;
        else if (previous_params_arities[i] == 4)
            tiny_expression->vars[i].address = cloture4;
        else if (previous_params_arities[i] == 5)
            tiny_expression->vars[i].address = cloture5;
        else if (previous_params_arities[i] == 6)
            tiny_expression->vars[i].address = cloture6;
        else if (previous_params_arities[i] == 7)
            tiny_expression->vars[i].address = cloture7;
        tiny_expression->vars[i].type    = previous_params_arities[i] == 0
            ? TE_VARIABLE
            : TE_CLOSURE0 + previous_params_arities[i];
        tiny_expression->vars[i].context = previous_params_arities[i] == 0
            ? (void*)0
            : &previous_params_functions[i];
    }
    for (int i = 0; i < tiny_expression->arity; i++) {
        tiny_expression->vars[nb_previous_params + i].name    = tiny_expression->variable_names[i];
        tiny_expression->vars[nb_previous_params + i].address = &tiny_expression->variable_values[i];
        tiny_expression->vars[nb_previous_params + i].type    = TE_VARIABLE;
        tiny_expression->vars[nb_previous_params + i].context = (void*)0;
    }
    
    /* This will compile the expression and check for errors. */
    int err;
    tiny_expression->expr = te_compile(expression, tiny_expression->vars, nb_parameters, &err);
    
    if (!tiny_expression->expr) {
        /* Show the user where the error is at. */
        ERROR_MESSAGE("#When evaluating:\n\t%s\n\t%*s^\n#There was an error near here.\n", expression, err-1, "");
    }
}

/*
 * Builds the parallel representation of an initial function provided by the user.
 */
void fun_1d1v_user_defined(PC_tree_t conf, parallel_stuff* par_variables, 
        double* array1, double *array2, char* fun_name) {
    // Getting the definition.
    char* definition;
    if (PC_get(conf, ".definition")) {
        PC_string(PC_get(conf, ".definition"), &definition);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'definition' field.\n", fun_name);
    }
    // Getting the variable string.
    char* variables_string = (char*)0;
    char* expression       = (char*)0;
    split_in_two(definition, DEFINITION_DELIM, &variables_string, &expression);
    int length = strlen(variables_string);
    if (variables_string[0] != '(' || variables_string[length-1] != ')') {
        ERROR_MESSAGE("#Error in function %s: the variable names must be enclosed in parentheses.\n", fun_name);
    }
    variables_string[length-1] = '\0';
    variables_string++;
    // Converting the string to names.
    char* variables[2];
    split_in_two(variables_string, VARIABLES_DELIM, &variables[0], &variables[1]);
    if (!variables[1] || strlen(variables[1]) == 0 || strlen(variables[0]) == 0) {
        ERROR_MESSAGE("#Error in function %s: too few variables (there should be 2).\n", fun_name);
        if (strstr(variables[1], VARIABLES_DELIM)) {
            ERROR_MESSAGE("#Error in function %s: too much variables (there should be 2).\n", fun_name);
        }
    }
    // Getting the parameters names, values, and expressions.
    // When the parameter has arity 0, its value is used. Otherwise, its expression is used.
    int nb_parameters = 0;
    char** parameter_names;
    double* parameter_values;
    int* parameter_arities;
    tiny_expr* parameter_functions;
    if (PC_get(conf, ".parameters")) {
        PC_tree_t parameters_tree = PC_get(conf, ".parameters")->first_child;
        // Getting the number of parameters.
        while (parameters_tree) {
            parameters_tree = parameters_tree->next_sibling;
            nb_parameters++;
        }
        if (nb_parameters == 0) {
            ERROR_MESSAGE("#Error in function %s: you must provide parameter names and values or remove the 'parameters' field.\n", fun_name);
        }
        // Getting the names and values of parameters.
        parameters_tree = PC_get(conf, ".parameters")->first_child;
        parameter_names     = malloc(nb_parameters * sizeof(char*));
        parameter_values    = malloc(nb_parameters * sizeof(double));
        parameter_arities   = malloc(nb_parameters * sizeof(int));
        parameter_functions = malloc(nb_parameters * sizeof(tiny_expr));
        int id_parameter = 0;
        while (parameters_tree) {
            check_valid_variable_name(parameters_tree->key, fun_name);
            if (!strcmp(parameters_tree->key, fun_name)) {
                ERROR_MESSAGE("#Error in function %s: no parameter can have the name of the function.\n", fun_name);
            }
            parameter_names[id_parameter] = parameters_tree->key;
            any_fun_from_string(parameters_tree->value, fun_name,
                id_parameter, parameter_names,
                parameter_values, parameter_arities, parameter_functions,
                &parameter_functions[id_parameter]);
            parameter_arities[id_parameter] = parameter_functions[id_parameter].arity;
            parameters_tree = parameters_tree->next_sibling;
            id_parameter++;
        }
    } else {
        // It is perfectly valid to have no parameter.
    }
    
    // Building the tinyexpression
    double x, v;
    te_variable* vars = malloc((2 + nb_parameters) * sizeof(te_variable));
    for (int id_parameter = 0; id_parameter < nb_parameters; id_parameter++) {
        vars[id_parameter].name    = parameter_names[id_parameter];
        if (parameter_arities[id_parameter] == 0)
            vars[id_parameter].address = &parameter_values[id_parameter];
        else if (parameter_arities[id_parameter] == 1)
            vars[id_parameter].address = cloture1;
        else if (parameter_arities[id_parameter] == 2)
            vars[id_parameter].address = cloture2;
        else if (parameter_arities[id_parameter] == 3)
            vars[id_parameter].address = cloture3;
        else if (parameter_arities[id_parameter] == 4)
            vars[id_parameter].address = cloture4;
        else if (parameter_arities[id_parameter] == 5)
            vars[id_parameter].address = cloture5;
        else if (parameter_arities[id_parameter] == 6)
            vars[id_parameter].address = cloture6;
        else if (parameter_arities[id_parameter] == 7)
            vars[id_parameter].address = cloture7;
        vars[id_parameter].type    = parameter_arities[id_parameter] == 0
            ? TE_VARIABLE
            : TE_CLOSURE0 + parameter_arities[id_parameter];
        vars[id_parameter].context = parameter_arities[id_parameter] == 0
            ? (void*)0
            : &parameter_functions[id_parameter];
    }
    vars[nb_parameters    ].name    = variables[0];
    vars[nb_parameters    ].address = &x;
    vars[nb_parameters    ].type    = TE_VARIABLE;
    vars[nb_parameters    ].context = (void*)0;
    vars[nb_parameters + 1].name    = variables[1];
    vars[nb_parameters + 1].address = &v;
    vars[nb_parameters + 1].type    = TE_VARIABLE;
    vars[nb_parameters + 1].context = (void*)0;
    
    /* This will compile the expression and check for errors. */
    int err;
    te_expr* n = te_compile(expression, vars, 2 + nb_parameters, &err);
    
    if (n) {
        par_variables->is_par_x = true;
        local_to_global_2d(par_variables, 0, 0);
        for (int i = 0; i < par_variables->size_x_par_x; i++) {
            for (int j = 0; j < par_variables->size_v_par_x; j++) {
                x = array1[i + par_variables->global_indices[0]];
                v = array2[j + par_variables->global_indices[1]];
                // Evaluation of all parameters.
                for (int id_parameter = 0; id_parameter < nb_parameters; id_parameter++) {
                    double param_value = te_eval(parameter_functions[id_parameter].expr);
                    parameter_values[id_parameter] = param_value;
                }
                par_variables->f_parallel_in_x[i][j] = te_eval(n);
            }
        }
        
        te_free(n);
    } else {
        /* Show the user where the error is at. */
        ERROR_MESSAGE("#When evaluating:\n\t%s\n\t%*s^\n#There was an error near here.\n", expression, err-1, "");
    }
    
    // Cleaning
    if (nb_parameters > 0) {
        free(parameter_names);
        free(parameter_values);
        free(parameter_arities);
        for (int id_parameter = 0; id_parameter < nb_parameters; id_parameter++) {
            free(parameter_functions[id_parameter].variable_names);
            free(parameter_functions[id_parameter].variable_values);
        }
        free(parameter_functions);
    }
    free(vars);
}

/*
 * Builds the parallel representation of the Landau 1d1v initial function.
 */
void fun_1d1v_landau(PC_tree_t conf, parallel_stuff* par_variables,
        double* array1, double* array2, char* fun_name) {
    double res;
    double mode, eps; // parameters
    double x, v;      // variables
    
    // Getting the parameters values from the yaml file.
    if (PC_get(conf, ".mode")) {
        PC_double(PC_get(conf, ".mode"), &mode);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'mode' field.\n", fun_name);
    }
    if (PC_get(conf, ".eps")) {
        PC_double(PC_get(conf, ".eps"), &eps);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'eps' field.\n", fun_name);
    }
    
    par_variables->is_par_x = true;
    local_to_global_2d(par_variables, 0, 0);
    for (int i = 0; i < par_variables->size_x_par_x; i++) {
        for (int j = 0; j < par_variables->size_v_par_x; j++) {
            x = array1[i + par_variables->global_indices[0]];
            v = array2[j + par_variables->global_indices[1]];
            res = (1. + eps * cos(mode * x)) * exp(-v * v / 2.) / sqrt(2. * M_PI);
            par_variables->f_parallel_in_x[i][j] = res;
        }
    }
}

void fun_1d1v_beam(PC_tree_t conf, parallel_stuff* par_variables,
        double* array1, double* array2, char* fun_name) {
    double res;
    double alpha; // parameters
    double x_plus,x_minus;
    double x, v;      // variables
    
    // Getting the parameters values from the yaml file.
    if (PC_get(conf, ".alpha")) {
        PC_double(PC_get(conf, ".alpha"), &alpha);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'alpha' field.\n", fun_name);
    }
    
    par_variables->is_par_x = true;
    local_to_global_2d(par_variables, 0, 0);
    for (int i = 0; i < par_variables->size_x_par_x; i++) {
        for (int j = 0; j < par_variables->size_v_par_x; j++) {
            x = array1[i + par_variables->global_indices[0]];
            v = array2[j + par_variables->global_indices[1]];
 		    x_plus = (x+1.2)/1e-10;//0.3;
            x_minus = (x-1.2)/1e-10;//0.3;
            res = 4./sqrt(2.*M_PI*alpha);
		    res = res*(0.5*erf(x_plus)-0.5*erf(x_minus));
            res = res*exp(-v*v/(2.*alpha));   
            par_variables->f_parallel_in_x[i][j] = res;
            if(fabs(res)>4){
            printf("%1.20lg\n",res);
            }
        }
    }
}



//   function sll_f_beam_initializer_2d( x, vx, params ) result(res)
//     sll_real64 :: res
//     sll_real64, intent(in) :: x
//     sll_real64, intent(in) :: vx
//  
//     sll_real64, dimension(:), intent(in) :: params
//     sll_real64 :: alpha
//     sll_real64 :: x_plus
//     sll_real64 :: x_minus
//     
// 
// !!$    if( .not. present(params) ) then
// !!$       print *, '#sll_f_beam_initializer_2d, error: the params array must ', &
// !!$            'be passed. params(1) = alpha'
// !!$       stop
// !!$    end if
//     SLL_ASSERT(size(params)>=1)
//     alpha = params(1) 
//     
//     x_plus = (x+1.2_f64)/0.3_f64
//     x_minus = (x-1.2_f64)/0.3_f64
//     res = 4._f64/sqrt(2._f64*sll_p_pi*alpha)
//     res = res*(0.5_f64*erf(x_plus)-0.5_f64*erf(x_minus))
//     res = res*exp(-vx**2/(2*alpha))    
// 
//   end function sll_f_beam_initializer_2d



void fun_1d1v_malkov(PC_tree_t conf, parallel_stuff* par_variables,
        double* array1, double* array2, char* fun_name) {
    double res;
    double rho0, z0, alpha2, a, ap,res0; // parameters
    double z, u;      // variables
    
    // Getting the parameters values from the yaml file.
    if (PC_get(conf, ".rho0")) {
        PC_double(PC_get(conf, ".rho0"), &rho0);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'rho0' field.\n", fun_name);
    }
    if (PC_get(conf, ".z0")) {
        PC_double(PC_get(conf, ".z0"), &z0);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'z0' field.\n", fun_name);
    }
    if (PC_get(conf, ".alpha2")) {
        PC_double(PC_get(conf, ".alpha2"), &alpha2);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'alpha2' field.\n", fun_name);
    }
    if (PC_get(conf, ".a")) {
        PC_double(PC_get(conf, ".a"), &a);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'a' field.\n", fun_name);
    }
    if (PC_get(conf, ".ap")) {
        PC_double(PC_get(conf, ".ap"), &ap);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'ap' field.\n", fun_name);
    }
    if (PC_get(conf, ".res0")) {
        PC_double(PC_get(conf, ".res0"), &res0);
    } else {
        ERROR_MESSAGE("#Error in function %s: missing the 'res0' field.\n", fun_name);
    }
    
    par_variables->is_par_x = true;
    local_to_global_2d(par_variables, 0, 0);
    for (int i = 0; i < par_variables->size_x_par_x; i++) {
        for (int j = 0; j < par_variables->size_v_par_x; j++) {
            z = array1[i + par_variables->global_indices[0]];
            z = 0.;
            u = array2[j + par_variables->global_indices[1]];
            res = alpha2*(z0*z0-z*z/(a*a))-a*a*(u-(ap/a)*z)*(u-(ap/a)*z);
            if(res>res0){
            	res = pow(res,-0.5);
            }else{
            	res = 0.;
            }
            res = (rho0/M_PI)*res;
            //res = (1. + eps * cos(mode * x)) * exp(-v * v / 2.) / sqrt(2. * M_PI);
            par_variables->f_parallel_in_x[i][j] = res;
        }
    }
}

/*
 * Builds the parallel representation of any initial function.
 */
void fun_1d1v(PC_tree_t conf, parallel_stuff* par_variables, 
        double* array1, double *array2) {
    if (PC_get(conf, ".Landau")) {
        fun_1d1v_landau(PC_get(conf, ".Landau"), par_variables, array1, array2, conf->key);
    } else if (PC_get(conf, ".Beam")) {
        fun_1d1v_beam(PC_get(conf, ".Beam"), par_variables, array1, array2, conf->key);
    } else if (PC_get(conf, ".Malkov")) {
        fun_1d1v_malkov(PC_get(conf, ".Malkov"), par_variables, array1, array2, conf->key);
    } else if (PC_get(conf, ".userDefined")) {
        fun_1d1v_user_defined(PC_get(conf, ".userDefined"), par_variables, array1, array2, conf->key);
    } else {
        ERROR_MESSAGE("#Error in function %s: unkwown function type.\n#Possible function types are 'Landau', 'Beam', 'Malkov' or 'userDefined'.\n", conf->key);
    }
}

#endif // ifndef SELA_VP_1D1V_CART_INIT_ANY_EXPR
