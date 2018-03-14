
#ifndef STRINGS_H
#define STRINGS_H

/* Here are just some utility routines for manipulating strings. */
// remove leading whitespace and update str to the first non-whitespace token.
// modifies the str inplace and returns it as well (so you can use it when
// chaining functions).
char *trim_whitespace(char *str);



// this is essentially strdup, but strdup is not standard C, so we use this...
char *string_copy(const char *s);

#endif
