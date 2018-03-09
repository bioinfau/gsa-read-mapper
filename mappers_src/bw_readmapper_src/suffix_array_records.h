
#ifndef SUFFIX_ARRAY_RECORDS_H
#define SUFFIX_ARRAY_RECORDS_H

#include "string_vector.h"
#include "size_vector.h"
#include "fasta.h"
#include <stdio.h>


struct suffix_array_records {
    struct string_vector *names;
    struct size_vector *suffix_arrays;
};

struct suffix_array_records *empty_suffix_array_records();
void delete_suffix_array_records(struct suffix_array_records *records);

int read_suffix_array_records(struct suffix_array_records *records,
                              struct fasta_records *fasta_records,
                              FILE *file);


#endif
