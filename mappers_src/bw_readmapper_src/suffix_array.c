
#include "suffix_array.h"
#include "strings.h"
#include "pair_stack.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

struct suffix_array *empty_suffix_array()
{
    struct suffix_array *sa =
        (struct suffix_array*)malloc(sizeof(struct suffix_array));
    sa->string = 0;
    sa->length = 0;
    sa->array = 0;
    
    sa->c_table = 0;
    
    return sa;
}

static struct suffix_array *allocate_sa(char *string)
{
    struct suffix_array *sa = empty_suffix_array();
    sa->string = strdup(string);
    sa->length = strlen(string);
    sa->array = (size_t*)malloc(sa->length * sizeof(size_t));
    
    sa->c_table = 0;
    
    return sa;
}

static // Wrapper of strcmp needed for qsort
int construction_cmpfunc(const void *a, const void *b)
{
    return strcmp(*(char **)a, *(char **)b);
}

struct suffix_array *qsort_sa_construction(char *string)
{
    struct suffix_array *sa = allocate_sa(string);
    
    char **suffixes = malloc(sa->length * sizeof(char *));
    for (int i = 0; i < sa->length; ++i)
        suffixes[i] = (char *)string + i;
    
    qsort(suffixes, sa->length, sizeof(char *), construction_cmpfunc);
    
    for (int i = 0; i < sa->length; i++)
        sa->array[i] = suffixes[i] - string;
    
    return sa;
}

void compute_c_table(struct suffix_array *sa)
{
    // I know we do not use all the characters, but this is easier
    // than preprocessing the string and matching it to a smaller set
    // and for the sizes of data I can handle, in any case, this won't
    // be the main problem.
    sa->c_table = (size_t*)calloc(C_TABLE_SIZE, sizeof(size_t));
    for (int i = 0; i < sa->length; i++) {
        size_t index = (size_t)sa->string[i];
        sa->c_table[index]++;
    }
}


void delete_suffix_array(struct suffix_array *sa)
{
    if (sa->string)  free(sa->string);
    if (sa->array)   free(sa->array);
    if (sa->c_table) free(sa->c_table);
    free(sa);
}

