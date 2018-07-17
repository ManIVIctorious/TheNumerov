
/* "controlfile.h" */
/*{{{
    This header provides:
    * the definition of a value for maximum line length
    * the prototype of the ParseControlFile() function
    * the definition of the struct "keywords":
        explanation of struct entries:
          keyword:
            The actual keyword string
          identifier:
            Can be used to individually identify the given struct
            by using integer values or chars (e.g. 'a')
            or to group a set of them, e.g. all integer values
          set:
            An integer which should be initialized to 0. The control
            file parsing function will increment this value every time
            the given keyword was set
          value:
            A char array of MaxLineLength. The values corresponding
            to the given keyword will be stored in this array for
            later conversion.

    Invocation:
      A pointer "keywordlist" of the first element of an array of struct
      keywords is passed as an argument of the parsing function.
      The last entry of this element has to be filled with zeros.

        struct keywords keywordlist[] = {
            { "first_keyword",   0,  'a',  0 },
            { "second_keyword",  0,  'b',  0 },
            { "third_keyword",   0,  'c',  0 },
            { 0 , 0 , 0 , 0 }
        };

//}}}*/

#ifndef _CONTROL_FILE_H
#define _CONTROL_FILE_H

#define _MaxLineLength_ 2048

struct keywords {
    const char * keyword;
    int          set;
    const int    identifier;
    char         value[_MaxLineLength_];
};

int ControlFileParser(char* filename, struct keywords* keywordlist, int verbosity_flag);

#endif
