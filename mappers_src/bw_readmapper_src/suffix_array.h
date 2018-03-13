
#ifndef SUFFIX_ARRAY_H
#define SUFFIX_ARRAY_H

#include <stddef.h>
#include <assert.h>

#define C_TABLE_SIZE 256

struct suffix_array {
    // length of the array
    size_t length;
    // the actual suffix array
    size_t *array;
    
    // used in bw search
    size_t  *c_table;
    size_t   c_table_no_symbols;
    char    *c_table_symbols;
    size_t  *c_table_symbols_inverse; // reverse map +1 (to recognize misses)
    size_t *o_table;
};

struct suffix_array *empty_suffix_array(void);
struct suffix_array *qsort_sa_construction(const char *string);

void compute_c_table(struct suffix_array *sa, const char *string);
void compute_o_table(struct suffix_array *sa, const char *string);

size_t o_table_index(struct suffix_array *sa, char symbol, size_t idx);

void delete_suffix_array(struct suffix_array *sa);

#endif
