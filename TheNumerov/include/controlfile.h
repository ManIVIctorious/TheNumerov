
/* "controlfile.h" */
/*{{{
    This header provides:
    * the definition of a value for maximum line length and maximum value length
    * the definition of the struct "keywords":
        struct entries:
          char * keyword:
            The actual keyword string
          int set:
            Should be initialised to zero.
            The tokeniser increments this value every time the given keyword was set
          const int identifier:
            Can be used to individually identify the given struct by using integer
            values or chars (e.g. 'a') or to group a set of them, e.g. all integer values
          char ** value:
            List of strings, should be initialised to NULL.
            Each time a given keyword-value pair is found the buffer gets extended by
            an additional entry in which the value string is stored for later use.
            The type conversion and setting of values is performed in the invoking function
            which is also responsible for the freeing of the allocated memory.
            For this purpose the helper function free_keywordlistvalues() is provided.

    * the prototype of the TokeniseControlFile() function
        This function reads the control file. After removing comments and empty lines it
        splits the buffer at every occurrence of semi colon ';'.
        This is then split again at the occurrence of an equality sign '=', resulting in
        a keyword and its respective value.
        Both are stripped from leading and trailing whitespaces and then the keyword is
        (case insensitively) compared to all entries of the keyword list.
        If the keyword was recognised its set counter is incremented and an additional
        buffer of size _MaxValLength_ is allocated for the value list.
        Last the value is copied to the new buffer.

        The function checks for various types of invalid entries and can be invoked with
        all warnings being treated as errors. The following things should be checked on
        invocation side:
            - Multiple use of keywords (all values are stored in the char ** value list)
            - The tokeniser allows for whitespaces within keywords and values.
              In case of the keywords this is handled by the removal of whitespaces in
              the keywords of the provided list.
              If whitespaces in the values are undesired this must be checked in the
              calling function.
            - Semi colon ';' is an invalid character for both, keyword and value
              (it is used as separator for multiple keywords per line)
            - The keywordlist should be freed after all conversions are preformed.
              For this purpose the helper function free_keywordlistvalues() is provided.

    Invocation:
      A pointer "keywordlist" of the first element of an array of struct
      keywords is passed as an argument of the parsing function.
      The last entry of this element has to be filled with zeros or NULL, respectively.
      (strictly speaking for the tokeniser to work only the first struct member, i.e.
      .keyword of the last entry has to be NULL)

        struct keywords keywordlist[] = {
            { "first_keyword",   0,  'a', NULL },
            { "second_keyword",  0,  'b', NULL },
            { "third_keyword",   0,  'c', NULL },
            { NULL , 0 , 0 , NULL }
        };

//}}}*/

#ifndef _CONTROL_FILE_H
#define _CONTROL_FILE_H

// maximum line and key-value length of control file
#define _MaxLineLength_ 2048
#define _MaxValLength_   512

struct keywords {
    char       * keyword;
    int          set;
    const int    identifier;
    char      ** value;
};

int  TokeniseControlFile(char* filename, struct keywords* kwlist, int ExitOnError);
void free_keywordlistvalues(struct keywords* kwlist);

#endif
