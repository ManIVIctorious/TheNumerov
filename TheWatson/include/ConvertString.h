
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

/*  splitstring_array() receives a string, the address to a char ** pointer
 *  which should be initialised to NULL, a list of delimiters as well as the
 *  preserve_array flag.
 *
 *  the function splits the string "array" at its delimiters "delimit" and
 *  stores the results in a (growing) array of strings "stringlist".
 *
 *  If the preserve_array option is non zero, the first <preserve_array>
 *  bytes of "array" will be copied to an auxiliary string, else the original
 *  string will be altered.
 */
int splitstring_to_array(char* array, char** *stringlist, const char* delimit, size_t preserve_array);

#endif
