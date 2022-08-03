#ifndef SELA_STRING_HELPERS
#define SELA_STRING_HELPERS

/*************************
 * Error messages helper *
 *************************/
#include <stdio.h>  // function fprintf (output strings on a stream)
                    // constant stderr (standard error output stream)
#include <stdlib.h> // function exit (error handling)
                    // constant EXIT_FAILURE (error handling)
                    // type     size_t
#include <string.h> // functions strlen, strstr
#include <ctype.h>  // function  isspace
#define ERROR_MESSAGE(...)            \
    do {                              \
        fprintf(stderr, __VA_ARGS__); \
        exit(EXIT_FAILURE);           \
    } while(0)


/**************************************************
 * Real number helpers                            *
 * (when using printf("%.*g", PRECISION, number). *
 **************************************************/
#include <float.h> // constant DBL_DECIMAL_DIG (17 : number of decimals a double has, only in C11)
                   // constant FLT_DECIMAL_DIG ( 9 : number of decimals a float  has, only in C11)
// According to https://en.wikipedia.org/wiki/Floating_point#Internal_representation
#if !defined(FLT_DECIMAL_DIG)
#    define FLT_DECIMAL_DIG 9
#endif
#if !defined(DBL_DECIMAL_DIG)
#    define DBL_DECIMAL_DIG 17
#endif


/************************
 * Other string helpers *
 ************************/
/*
 * Removes leading and trailing spaces. Changes are made in the memory
 * allocated for str.
 * https://stackoverflow.com/questions/122616/how-do-i-trim-leading-trailing-whitespace-in-a-standard-way
 *
 * @param[in] str the initial string
 * @return    str without leading and trailing spaces.
 */
char* trim_spaces(char* str) {
    size_t len = 0;
    char* frontp = str;
    char* endp = (void*)0;
    
    if (str == (void*)0) { return (void*)0; }
    if (str[0] == '\0')  { return str; }
    
    len = strlen(str);
    endp = str + len;
    
    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while (isspace((unsigned char) *frontp)) { ++frontp; }
    if (endp != frontp) {
        while (isspace((unsigned char) *(--endp)) && endp != frontp) {}
    }
    
    if (frontp != str && endp == frontp) {
        *str = '\0';
    } else if (str + len - 1 != endp) {
        *(endp + 1) = '\0';
    }
    
    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer. Note the reuse
     * of endp to mean the front of the string buffer now.
     */
    endp = str;
    if (frontp != str) {
        while (*frontp) { *endp++ = *frontp++; }
        *endp = '\0';
    }
    
    return str;
}

/*
 * Split a string in two with respect to a delimiter.
 * If the delimiter is "@" and the string is "firstname.surname@site.com", then
 * *before will be "firstname.surname", *after will be "site.com".
 * If the delimiter is not found, then *before will be str and *after will be (char*)0.
 *
 * @param[in]  str    the string to split
 * @param[in]  delim  the delimiter
 * @param[out] before the portion of the string before the delimiter
 * @param[out] after  the portion of the string after the delimiter
 */
void split_in_two(char* str, const char* delim, char** before, char** after) {
    *before = str;
    *after = strstr(str, delim); // Points where delim is, if it belongs to this string
    if (!*after) {               // delim was not found in the string
        return;
    }
    **after = '\0';                  // Ending the string here for *before
    *after = *after + strlen(delim);
    *after = trim_spaces(*after);    // Remove leading and trailing spaces (after delim)
    *before = trim_spaces(*before);  // Remove leading and trailing spaces (before delim)
}

/*
 * Count the number of times a given character appears in a given string.
 * https://stackoverflow.com/questions/7349053/counting-the-number-of-times-a-character-occurs-in-a-string-in-c
 *
 * @param[in] str the string in which to search for c.
 * @param[in] c the character to search for.
 * @return    the number of times c is appears in str.
 */
int count_chars(const char* str, char c) {
    int count = 0;
    for(; *str; count += (*str++ == c)) ;
    return count;
}

#endif // ifndef SELA_STRING_HELPERS
