
#include "suffix_array_records.h"
#include <stdlib.h>
#include <string.h>

struct suffix_array_records *empty_suffix_array_records()
{
    struct suffix_array_records *records =
        (struct suffix_array_records*)malloc(sizeof(struct suffix_array_records));
    records->names = empty_string_vector(10); // arbitrary size...
    records->suffix_arrays = 0;
    return records;
}

void delete_suffix_array_records(struct suffix_array_records *records)
{
    delete_string_vector(records->names);
    if (records->suffix_arrays) {
        for (int i = 0; i < records->names->used; i++) {
            free(records->suffix_arrays[i]);
        }
        free(records->suffix_arrays);
    }
    free(records);
}

#define NAME_BUFFER_SIZE 1024
int read_suffix_array_records(struct suffix_array_records *records,
                              struct fasta_records *fasta_records,
                              FILE *file)
{
    int no_records = fasta_records->names->used;
    records->suffix_arrays = malloc(sizeof(size_t) * no_records);
    char seq_name[NAME_BUFFER_SIZE];
    for (int i = 0; i < fasta_records->names->used; i++) {
        fscanf(file, "%1024s", (char*)&seq_name);
        if (strcmp(seq_name, fasta_records->names->strings[i]) != 0) {
            fprintf(stderr, "The preprocessed sequence read is %s while the FASTA record is %s. This is an error!\n",
                    seq_name, fasta_records->names->strings[i]);
            return 1;
        }
        add_string_copy(records->names, seq_name);
        size_t *sa = records->suffix_arrays[i] = malloc(sizeof(size_t) * fasta_records->seq_sizes->sizes[i]);
        for (int j = 0; j < fasta_records->seq_sizes->sizes[i]; j++) {
            size_t index;
            fscanf(file, "%lu", &index);
            sa[j] = index;
        }
    }
    return 0;
}
