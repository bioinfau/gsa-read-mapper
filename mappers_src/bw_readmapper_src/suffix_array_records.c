
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
        const char *string = fasta_records->sequences->strings[i];
        add_string_copy(records->names, seq_name);
        fprintf(stderr, "building suffix array for %s.\n", seq_name);
        records->suffix_arrays[i] = qsort_sa_construction(
            fasta_records->sequences->strings[i]
        );
        fprintf(stderr, "building c-table for %s.\n", seq_name);
        compute_c_table(records->suffix_arrays[i], string);
        
        fprintf(stderr, "building o-table for %s.\n", seq_name);
        compute_o_table(records->suffix_arrays[i], string);
    }
    fprintf(stderr, "Done.\n");
    
    return records;
}

void delete_suffix_array_records(struct suffix_array_records *records)
{
    
    if (records->suffix_arrays) {
        for (int i = 0; i < records->names->used; i++) {
            fprintf(stderr, "deleting suffix array for sequence %s.\n",
                    records->names->strings[i]);
            delete_suffix_array(records->suffix_arrays[i]);
        }
        free(records->suffix_arrays);
    }
    delete_string_vector(records->names);
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
        fprintf(sa_file, " %lu", sa->length);
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
        size_t  c_table_no_symbols = records->suffix_arrays[i]->c_table_no_symbols;
        char    *c_table_symbols = records->suffix_arrays[i]->c_table_symbols;
        
        assert(c_table);
        assert(c_table_no_symbols);
        assert(c_table_symbols);
        
        fprintf(sa_file, "%s", records->names->strings[i]);
        fprintf(sa_file, " %lu", c_table_no_symbols);
        for (int j = 0; j < c_table_no_symbols; j++) {
            char symbol = c_table_symbols[j];
            fprintf(sa_file, " %c %lu", symbol, c_table[(size_t)symbol]);
        }
        fprintf(sa_file, "\n");
    }
    fclose(sa_file);
    free(filename);
    
    filename = make_file_name(filename_prefix, "o_tables");
    fprintf(stderr, "writing o-table to %s.\n", filename);
    sa_file = fopen(filename, "w");
    for (int i = 0; i < records->names->used; i++) {
        size_t *o_table = records->suffix_arrays[i]->o_table;
        size_t  o_table_size = records->suffix_arrays[i]->c_table_no_symbols * records->suffix_arrays[i]->length;
        
        assert(o_table);
        assert(o_table_size);
        
        fprintf(sa_file, "%s", records->names->strings[i]);
        fprintf(sa_file, " %lu", o_table_size);
        for (size_t i = 0; i < o_table_size; i++) {
            fprintf(sa_file, " %lu", o_table[i]);
        }
        fprintf(sa_file, "\n");
    }
    fclose(sa_file);
    free(filename);
    
    fprintf(stderr, "Done.\n");
    
    return 0;
}

#define NAME_BUFFER_SIZE 1024
static int read_o_table_records(struct suffix_array_records *records,
                                struct fasta_records *fasta_records,
                                const char *filename_prefix)
{
    char *filename = make_file_name(filename_prefix, "o_tables");
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open file %s.\n", filename);
        exit(1);
    }
    fprintf(stderr, "Reading o-table from %s.\n", filename);
    
    int no_records = fasta_records->names->used;
    for (int i = 0; i < no_records; i++) {
        char seq_name[NAME_BUFFER_SIZE];
        fscanf(file, "%1024s", (char*)&seq_name);
        if (strcmp(seq_name, fasta_records->names->strings[i]) != 0) {
            fprintf(stderr, "The preprocessed o-table sequence read is %s while the FASTA record is %s. This is an error!\n",
                    seq_name, fasta_records->names->strings[i]);
            return 1;
        }
        fprintf(stderr, "reading o-table for sequence %s.\n", seq_name);

        size_t no_symbols = records->suffix_arrays[i]->c_table_no_symbols;
        size_t string_length = records->suffix_arrays[i]->length;
        size_t o_table_size;
        fscanf(file, "%lu", &o_table_size);
        fprintf(stderr, "... size %lu [%lu x %lu].\n",
                o_table_size, no_symbols, string_length);
        if (o_table_size != no_symbols * string_length) {
            fprintf(stderr, "Unexpected size of o-table!\n");
            exit(1);
        }
        
        size_t *o_table = malloc(sizeof(size_t)*o_table_size);
        for (size_t i = 0; i < o_table_size; i++) {
            size_t val;
            fscanf(file, "%lu", &val);
            o_table[i] = val;
        }
        
        if (records->suffix_arrays[i]->o_table)
            free(records->suffix_arrays[i]->o_table);
        records->suffix_arrays[i]->o_table = o_table;
        
#if 0
        struct suffix_array *sa = records->suffix_arrays[i];
        for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
            char symbol = sa->c_table_symbols[i];
            printf("O(%c,) =", (symbol == 0) ? '$' : symbol);
            for (size_t j = 0; j < sa->length; ++j) {
                size_t idx = o_table_index(sa, symbol, j);
                printf(" %lu", sa->o_table[idx]);
            }
            printf("\n");
        }
#endif

    }
    
    fprintf(stderr, "Done.\n");
    
    
    
    return 0;
}

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
        size_t c_table_size;
        fscanf(file, "%lu", &c_table_size);
        fprintf(stderr, "... contains %lu non-zero records.\n",
                c_table_size);
        
        size_t *c_table = calloc(256, sizeof(size_t));
        size_t  c_table_no_symbols = c_table_size;
        size_t  current_symbol_index = 0;
        char   *c_table_symbols = malloc(c_table_size);
        
        for (int j = 0; j < c_table_size; j++) {
            char symbol[NAME_BUFFER_SIZE]; size_t count;
            fscanf(file, "%1024s %lu", (char*)&symbol, &count);
            // the symbol should be '\0' or a single character
            assert(strlen(symbol) <= 1);
            char symbol_c = symbol[0];
            c_table[(size_t)symbol_c] = count;
            c_table_symbols[current_symbol_index++] = symbol_c;
            fprintf(stderr, "... %c -> %lu\n",
                    (symbol_c == 0) ? '$' : symbol_c, count);
        }
        
        int *c_table_symbols_inverse = calloc(C_TABLE_SIZE, sizeof(int));
        for (int j = 0; j < c_table_no_symbols; j++) {
            char symbol = c_table_symbols[j];
            c_table_symbols_inverse[symbol] = j + 1;
        }
        for (int j = 0; j < c_table_no_symbols; j++) {
            char symbol = c_table_symbols[j];
            int index = c_table_symbols_inverse[symbol];
            fprintf(stderr, "symbol %c has index %d.\n",
                    (symbol == '\0') ? '$' : symbol, index - 1);
        }
        
        if (records->suffix_arrays[i]->c_table)
            free(records->suffix_arrays[i]->c_table);
        records->suffix_arrays[i]->c_table = c_table;
        if (records->suffix_arrays[i]->c_table_symbols)
            free(records->suffix_arrays[i]->c_table_symbols);
        records->suffix_arrays[i]->c_table_no_symbols = c_table_no_symbols;
        records->suffix_arrays[i]->c_table_symbols = c_table_symbols;
        if (records->suffix_arrays[i]->c_table_symbols_inverse)
            free(records->suffix_arrays[i]->c_table_symbols_inverse);
        records->suffix_arrays[i]->c_table_symbols_inverse = c_table_symbols_inverse;
        
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
        size_t length = fasta_records->seq_sizes->sizes[i] + 1;
        size_t file_length;
        fscanf(file, "%lu", &file_length);
        if (length != file_length) {
            fprintf(stderr, "We expected to see a suffix array of length %lu but find one of length %lu.\n", length, file_length);
            exit(1);
        }
        sa->length = length;
        sa->array = malloc(sizeof(size_t) * length);
        for (int j = 0; j < length; j++) {
            size_t index;
            fscanf(file, "%lu", &index);
            sa->array[j] = index;
        }
        records->suffix_arrays[i] = sa;
        
#if 0
        const char *string = fasta_records->sequences->strings[i];
        fprintf(stderr, "suffix array:\n");
        for (int j = 0; j < sa->length; j++) {
            fprintf(stderr, "sa[%3d] == %4lu %s\n",
                    j, sa->array[j], string + sa->array[j]);
        }
#endif

    }
    
    fclose(file);
    free(filename);
    
    fprintf(stderr, "Done.\n");
    
    

    read_c_table_records(records, fasta_records, filename_prefix);
    read_o_table_records(records, fasta_records, filename_prefix);
    
    return 0;
}


