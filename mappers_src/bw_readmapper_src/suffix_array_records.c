
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

struct suffix_array_records *build_suffix_array_records(struct fasta_records *fasta_records)
{
    size_t no_records = fasta_records->names->used;
    struct suffix_array_records *records = empty_suffix_array_records();
    records->suffix_arrays = (struct suffix_array **)malloc(sizeof(struct suffix_array*)*no_records);
    
    fprintf(stderr, "Building suffix arrays.\n");
    for (int i = 0; i < no_records; i++) {
        const char *seq_name = fasta_records->names->strings[i];
        add_string_copy(records->names, seq_name);
        fprintf(stderr, "building suffix array for %s.\n", seq_name);
        records->suffix_arrays[i] = qsort_sa_construction(
            fasta_records->sequences->strings[i]
        );
        fprintf(stderr, "building c-table for %s.\n", seq_name);
        compute_c_table(records->suffix_arrays[i]);
    }
    fprintf(stderr, "Done.\n");
    
    return records;
}

void delete_suffix_array_records(struct suffix_array_records *records)
{
    delete_string_vector(records->names);
    if (records->suffix_arrays) {
        for (int i = 0; i < records->names->used; i++) {
            delete_suffix_array(records->suffix_arrays[i]);
        }
        free(records->suffix_arrays);
    }
    free(records);
}

static char *make_file_name(const char *prefix, const char *suffix) {
    size_t prefix_length = strlen(prefix);
    size_t suffix_length = strlen(suffix);
    size_t string_length = prefix_length + 1 + suffix_length + 1;
    
    char *buffer = (char*)malloc(string_length);
    char *c = buffer;
    for (int i = 0; i < prefix_length; i++, c++) {
        *c = prefix[i];
    }
    *c = '.'; c++;
    for (int i = 0; i < suffix_length; i++, c++) {
        *c = suffix[i];
    }
    *c = 0;
    
    return buffer;
}

int write_suffix_array_records(struct suffix_array_records *records,
                               struct fasta_records *fasta_records,
                               const char *filename_prefix)
{
    fprintf(stderr, "Writig preprocessed data to files.\n");
    char *filename = make_file_name(filename_prefix, "suffix_arrays");
    fprintf(stderr, "writing suffix array to %s.\n", filename);
    FILE *sa_file = fopen(filename, "w");
    for (int i = 0; i < records->names->used; i++) {
        struct suffix_array *sa = records->suffix_arrays[i];
        fprintf(sa_file, "%s", records->names->strings[i]);
        for (int j = 0; j < sa->length; j++) {
            fprintf(sa_file, " %lu", sa->array[j]);
        }
        fprintf(sa_file, "\n");
    }
    fclose(sa_file);
    free(filename);
    
    filename = make_file_name(filename_prefix, "c_tables");
    fprintf(stderr, "writing c-table to %s.\n", filename);
    sa_file = fopen(filename, "w");
    for (int i = 0; i < records->names->used; i++) {
        size_t *c_table = records->suffix_arrays[i]->c_table;
        assert(c_table);
        size_t no_non_empty = 0;
        for (int j = 0; j < C_TABLE_SIZE; j++) {
            if (c_table[j] != 0) no_non_empty++;
        }
        fprintf(sa_file, "%s", records->names->strings[i]);
        fprintf(sa_file, " %lu", no_non_empty);
        for (int j = 0; j < C_TABLE_SIZE; j++) {
            if (c_table[j] != 0)
                fprintf(sa_file, " %c %lu", (char)j, c_table[j]);
        }
        fprintf(sa_file, "\n");
    }
    fclose(sa_file);
    free(filename);
    
    fprintf(stderr, "Done.\n");
    
    return 0;
}

#define NAME_BUFFER_SIZE 1024
static int read_c_table_records(struct suffix_array_records *records,
                                struct fasta_records *fasta_records,
                                const char *filename_prefix)
{
    char *filename = make_file_name(filename_prefix, "c_tables");
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open file %s.\n", filename);
        exit(1);
    }
    fprintf(stderr, "Reading c-table from %s.\n", filename);
    
    int no_records = fasta_records->names->used;
    for (int i = 0; i < no_records; i++) {
        char seq_name[NAME_BUFFER_SIZE];
        fscanf(file, "%1024s", (char*)&seq_name);
        if (strcmp(seq_name, fasta_records->names->strings[i]) != 0) {
            fprintf(stderr, "The preprocessed c-table sequence read is %s while the FASTA record is %s. This is an error!\n",
                    seq_name, fasta_records->names->strings[i]);
            return 1;
        }
        fprintf(stderr, "reading c-table for sequence %s.\n", seq_name);
        size_t *c_table = calloc(256, sizeof(size_t));
        size_t c_table_size;
        fscanf(file, "%lu", &c_table_size);
        fprintf(stderr, "... contains %lu non-zero records.\n",
                c_table_size);
        
        for (int j = 0; j < c_table_size; j++) {
            char symbol[NAME_BUFFER_SIZE]; size_t count;
            fscanf(file, "%1024s %lu", (char*)&symbol, &count);
            assert(strlen(symbol) == 1);
            char symbol_c = symbol[0];
            c_table[(size_t)symbol_c] = count;
            fprintf(stderr, "... %c -> %lu\n", symbol_c, count);
        }
        
        if (records->suffix_arrays[i]->c_table)
            free(records->suffix_arrays[i]->c_table);
        records->suffix_arrays[i]->c_table = c_table;
    }
    
    fprintf(stderr, "Done.\n");
    
    return 0;
}


int read_suffix_array_records(struct suffix_array_records *records,
                              struct fasta_records *fasta_records,
                              const char *filename_prefix)
{
    char *filename = make_file_name(filename_prefix, "suffix_arrays");
    FILE *file = fopen(filename, "r");
    
    if (!file) {
        fprintf(stderr, "Could not open file %s.\n", filename);
        exit(1);
    }
    fprintf(stderr, "Reading suffix arrays from %s.\n", filename);
    
    int no_records = fasta_records->names->used;
    assert(records->suffix_arrays == 0);
    
    records->suffix_arrays =
        (struct suffix_array**)malloc(sizeof(struct suffix_array*) * no_records);
    
    char seq_name[NAME_BUFFER_SIZE];
    for (int i = 0; i < fasta_records->names->used; i++) {
        fscanf(file, "%1024s", (char*)&seq_name);
        if (strcmp(seq_name, fasta_records->names->strings[i]) != 0) {
            fprintf(stderr, "The preprocessed sequence read is %s while the FASTA record is %s. This is an error!\n",
                    seq_name, fasta_records->names->strings[i]);
            return 1;
        }
        fprintf(stderr, "reading suffix array for sequence %s.\n",
                seq_name);
        add_string_copy(records->names, seq_name);
        struct suffix_array *sa = empty_suffix_array();
        size_t length = fasta_records->seq_sizes->sizes[i];
        sa->length = length;
        sa->array = malloc(sizeof(size_t) * length);
        for (int j = 0; j < length; j++) {
            size_t index;
            fscanf(file, "%lu", &index);
            sa->array[j] = index;
        }
        records->suffix_arrays[i] = sa;
    }
    
    fclose(file);
    free(filename);
    
    fprintf(stderr, "Done.\n");
    
    read_c_table_records(records, fasta_records, filename_prefix);
    
    return 0;
}


