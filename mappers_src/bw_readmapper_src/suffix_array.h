
#ifndef SUFFIX_ARRAY_H
#define SUFFIX_ARRAY_H

#include <stddef.h>
#include <assert.h>

#define C_TABLE_SIZE 256

struct suffix_array {
    // the suffix array owns this, so copy it if you want to keep it.
    // it will be freed when the suffix array is freed.
    char *string;
    // length of the array
    size_t length;
    // the actual suffix array
    size_t *array;
    
    // used in bw search
    size_t *c_table;
};

struct suffix_array *empty_suffix_array(void);
struct suffix_array *qsort_sa_construction(char *string);
void compute_c_table(struct suffix_array *sa);

void delete_suffix_array(struct suffix_array *sa);

#endif
