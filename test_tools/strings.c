
#include "strings.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// this is essentially strdup, but strdup is not standard C, so we use this...
char* string_copy(const char* s)
{
    size_t n = strlen(s) + 1;
    char* copy = (char*)malloc(n);
    strcpy(copy, s);
    return copy;
}

char* trim_whitespace(char* str)
{
    // find first non-whitespace
    char* beg = str;
    while (isspace((unsigned char)*beg))
        beg++;

    if (*beg == 0)
    { // All spaces?
        // make the string empty.
        *str = '\0';
        return str;
    }

    // find end of non-whitespace
    char* end = beg;
    while (*end && !isspace((unsigned char)*end))
        end++;

    // Write new null terminator, in case we haven't reached the end of str
    *(end + 1) = 0;

    // move the non-whitespace token to the front of str.
    char* dst;
    for (dst = str; *beg; beg++, dst++)
        *dst = *beg;
    *dst = '\0';

    return str;
}
