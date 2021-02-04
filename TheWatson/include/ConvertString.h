
#ifndef __CONVERT_TO
#define __CONVERT_TO

/*  Convert string to bool, long integer or double
 *
 *    All three functions get two char pointers and an integer pointer
 *    as arguments (optarg, varname and control)
 *
 *    optarg:   the string to be converted to a number
 *
 *    varname:  optional variable name, if not NULL print it in front
 *              of possible errors for easier identification.
 *
 *    control:  if control is NULL let the function handle the error,
 *              else set its value to errvar (either errno or):
 *                  1     No conversion
 *                  2     Only part of string was converted to number
 *
 *    Depending on control the function either returns or exits (NULL)
 *    the program on first encountered error therefore, if control is
 *    set varname has no effect.
 */

long int convertstring_to_bool(char* optarg, char* varname, int *control);
long int convertstring_to_long(char* optarg, char* varname, int *control);
double convertstring_to_double(char* optarg, char* varname, int *control);

#endif
