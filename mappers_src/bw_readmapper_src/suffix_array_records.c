
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
    for (size_t i = 0; i < no_records; i++) {
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
        for (size_t i = 0; i < records->names->used; i++) {
            fprintf(stderr, "deleting suffix array for sequence %s.\n",
                    records->names->strings[i]);
            delete_suffix_array(records->suffix_arrays[i]);
        }
        free(records->suffix_arrays);
    }
    delete_string_vector(records->names);
    free(records);
}

static char *make_file_name(const char *prefix,
                            const char *suffix,
                            const char *seq_suffix) {
    size_t prefix_length = strlen(prefix);
    size_t suffix_length = strlen(suffix);
    size_t string_length = prefix_length + 1 + suffix_length + 1;

    size_t seq_suffix_length = 0;
    if (seq_suffix) {
        seq_suffix_length = strlen(seq_suffix);
        string_length += seq_suffix_length + 1;
    }
    
    char *buffer = (char*)malloc(string_length);
    char *c = buffer;
    for (size_t i = 0; i < prefix_length; i++, c++) {
        *c = prefix[i];
    }
    *c = '.'; c++;
    for (size_t i = 0; i < suffix_length; i++, c++) {
        *c = suffix[i];
    }
    if (seq_suffix) {
        *c = '.'; c++;
        for (size_t i = 0; i < seq_suffix_length; i++, c++) {
            *c = seq_suffix[i];
        }
    }
    *c = 0;
    
    assert(strlen(buffer) + 1 == string_length);
    
    return buffer;
}

int write_suffix_array_records(struct suffix_array_records *records,
                               struct fasta_records *fasta_records,
                               const char *filename_prefix)
{
    fprintf(stderr, "Writing preprocessed data to files.\n");
    
    for (size_t i = 0; i < records->names->used; i++) {
        const char *seq_name = records->names->strings[i];
        struct suffix_array *sa = records->suffix_arrays[i];
        
        char *filename = make_file_name(filename_prefix, "suffix_arrays", seq_name);
        
        fprintf(stderr, "writing suffix array to %s.\n", filename);
        FILE *file = fopen(filename, "wb");
        free(filename);
        
        fwrite(sa->array, sizeof(size_t), sa->length, file);
        fclose(file);
    }
    
    char *filename = make_file_name(filename_prefix, "c_tables", 0);
    fprintf(stderr, "writing c-table to %s.\n", filename);
    FILE *sa_file = fopen(filename, "w");
    for (size_t i = 0; i < records->names->used; i++) {
        size_t *c_table = records->suffix_arrays[i]->c_table;
        size_t  c_table_no_symbols = records->suffix_arrays[i]->c_table_no_symbols;
        char    *c_table_symbols = records->suffix_arrays[i]->c_table_symbols;
        
        assert(c_table);
        assert(c_table_no_symbols);
        assert(c_table_symbols);
        
        fprintf(sa_file, "%s", records->names->strings[i]);
        fprintf(sa_file, " %lu", c_table_no_symbols);
        for (size_t j = 0; j < c_table_no_symbols; j++) {
            char symbol = c_table_symbols[j];
            fprintf(sa_file, " %c %lu", symbol, c_table[(size_t)symbol]);
        }
        fprintf(sa_file, "\n");
    }
    fclose(sa_file);
    free(filename);

    for (size_t i = 0; i < records->names->used; i++) {
        const char *seq_name = records->names->strings[i];
        
        size_t *o_table = records->suffix_arrays[i]->o_table;
        size_t  o_table_size = records->suffix_arrays[i]->c_table_no_symbols * records->suffix_arrays[i]->length;
        
        char *filename = make_file_name(filename_prefix, "o_tables", seq_name);
        fprintf(stderr, "writing o-table to %s.\n", filename);
        
        FILE *file = fopen(filename, "wb");
        free(filename);
        fwrite(o_table, sizeof(size_t), o_table_size, file);
        fclose(file);
        
    }
    
    fprintf(stderr, "Done.\n");
    
    return 0;
}

static int read_o_table_records(struct suffix_array_records *records,
                                struct fasta_records *fasta_records,
                                const char *filename_prefix)
{
    assert(records->suffix_arrays != 0);
    
    for (size_t i = 0; i < fasta_records->names->used; i++) {
        const char *seq_name = fasta_records->names->strings[i];
        struct suffix_array *sa = records->suffix_arrays[i];
        char *filename = make_file_name(filename_prefix, "o_tables", seq_name);
        
        FILE *file = fopen(filename, "rb");
        if (!file) {
            fprintf(stderr, "Could not open file %s.\n", filename);
            exit(1);
        }
        fprintf(stderr, "reading o-table from %s.\n", filename);
        
        size_t  o_table_size = records->suffix_arrays[i]->c_table_no_symbols * records->suffix_arrays[i]->length;
        fprintf(stderr, "...allocating o-table size: [%lu x %lu] (%lu)\n",
                sa->c_table_no_symbols, sa->length,
                o_table_size);
        assert(o_table_size > 0);
        sa->o_table = malloc(o_table_size * sizeof(size_t));
        if (!sa->o_table) {
            fprintf(stderr, "...could not allocate memory for o-table.\n");
            exit(1);
        }
        fread(sa->o_table, sizeof(size_t), o_table_size, file);
        
        fclose(file);
        free(filename);

#if 0
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

#define NAME_BUFFER_SIZE 1024
static int read_c_table_records(struct suffix_array_records *records,
                                struct fasta_records *fasta_records,
                                const char *filename_prefix)
{
    char *filename = make_file_name(filename_prefix, "c_tables", 0);
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Could not open file %s.\n", filename);
        exit(1);
    }
    fprintf(stderr, "Reading c-table from %s.\n", filename);
    
    size_t no_records = fasta_records->names->used;
    for (size_t i = 0; i < no_records; i++) {
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
        
        for (size_t j = 0; j < c_table_size; j++) {
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
        
        size_t *c_table_symbols_inverse = calloc(C_TABLE_SIZE, sizeof(size_t));
        for (size_t j = 0; j < c_table_no_symbols; j++) {
            char symbol = c_table_symbols[j];
            c_table_symbols_inverse[(int)symbol] = j + 1;
        }
        for (size_t j = 0; j < c_table_no_symbols; j++) {
            char symbol = c_table_symbols[j];
            size_t index = c_table_symbols_inverse[(int)symbol];
            fprintf(stderr, "symbol %c has index %lu.\n",
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
    size_t no_records = fasta_records->names->used;
    assert(records->suffix_arrays == 0);
    
    records->suffix_arrays =
        (struct suffix_array**)malloc(sizeof(struct suffix_array*) * no_records);
    
    for (size_t i = 0; i < fasta_records->names->used; i++) {
        const char *seq_name = fasta_records->names->strings[i];
        size_t seq_length = fasta_records->seq_sizes->sizes[i];
        struct suffix_array *sa =
            records->suffix_arrays[i] =
            empty_suffix_array();
        
        char *filename = make_file_name(filename_prefix, "suffix_arrays", seq_name);
        FILE *file = fopen(filename, "rb");
        if (!file) {
            fprintf(stderr, "Could not open file %s.\n", filename);
            exit(1);
        }
        
        sa->length = seq_length + 1;
        sa->array = malloc(sizeof(size_t) * sa->length);
    
        fprintf(stderr, "Reading suffix arrays from %s [length %lu].\n",
                filename, sa->length);
        fread(sa->array, sizeof(size_t), sa->length, file);
        
        fclose(file);
        free(filename);

    }
    

    fprintf(stderr, "Done.\n");
    
    read_c_table_records(records, fasta_records, filename_prefix);
    read_o_table_records(records, fasta_records, filename_prefix);

    return 0;
}


