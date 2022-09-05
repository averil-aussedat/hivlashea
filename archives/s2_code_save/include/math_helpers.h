#ifndef SELA_MATH_HELPERS
#define SELA_MATH_HELPERS

#define DMIN(a,b) ((a)<=(b)?(a):(b))
#define DMAX(a,b) ((a)<=(b)?(b):(a))

/*
 * Computes the maximum of two doubles.
 *
 * @param  a, b.
 * @return max(a, b).
 */
static inline double max(double a, double b) {
    return a > b ? a : b;
}

/*
 * Computes the minimum of two doubles.
 *
 * @param  a, b.
 * @return min(a, b).
 */
static inline double min(double a, double b) {
    return a < b ? a : b;
}

/*
 * Computes the maximum of two ints.
 *
 * @param  a, b.
 * @return max(a, b).
 */
static inline int imax(int a, int b) {
    return a > b ? a : b;
}

/*
 * Computes the minimum of two ints.
 *
 * @param  a, b.
 * @return min(a, b).
 */
static inline int imin(int a, int b) {
    return a < b ? a : b;
}

/*
 * Integer exponentiation.
 *
 * @param[in] base is the number to get exponentiated.
 * @param[in] exp is the exponent
 * @return    base**exp = base^exp
 */
int int_pow(int base, int exp) {
    int result = 1;
    while (exp) {
        if (exp & 1)
            result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

#endif // ifndef SELA_MATH_HELPERS
