
#ifndef SUFFIX_ARRAY_RECORDS_H
#define SUFFIX_ARRAY_RECORDS_H

#include <stdio.h>
#include "fasta.h"
#include "suffix_array.h"

struct suffix_array_records {
    struct string_vector *names;
    struct suffix_array **suffix_arrays;
};

struct suffix_array_records *empty_suffix_array_records(void);
struct suffix_array_records *build_suffix_array_records(struct fasta_records *fasta_records);

void delete_suffix_array_records(struct suffix_array_records *records);

int write_suffix_array_records(struct suffix_array_records *records,
                               struct fasta_records *fasta_records,
                               const char *filename_prefix);

int read_suffix_array_records(struct suffix_array_records *records,
                              struct fasta_records *fasta_records,
                              const char *filename_prefix);


#endif
