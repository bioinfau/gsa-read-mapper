
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
    sa->length = 0;
    sa->array = 0;
    
    sa->c_table = 0;
    sa->c_table_no_symbols = 0;
    sa->c_table_symbols = 0;
    sa->c_table_symbols_inverse = 0;
    
    sa->o_table = 0;
    
    return sa;
}

static struct suffix_array *allocate_sa(const char *string)
{
    struct suffix_array *sa = empty_suffix_array();
    sa->length = strlen(string) + 1; // + 1 for empty string
    sa->array = (size_t*)malloc(sa->length * sizeof(size_t));
    
    return sa;
}

static // Wrapper of strcmp needed for qsort
int construction_cmpfunc(const void *a, const void *b)
{
    return strcmp(*(char **)a, *(char **)b);
}

struct suffix_array *qsort_sa_construction(const char *string)
{
    struct suffix_array *sa = allocate_sa(string);
    
    char **suffixes = malloc(sa->length * sizeof(char *));
    for (size_t i = 0; i < sa->length; ++i)
        suffixes[i] = (char *)string + i;
    
    qsort(suffixes, sa->length, sizeof(char *), construction_cmpfunc);
    
    for (size_t i = 0; i < sa->length; i++)
        sa->array[i] = (size_t)(suffixes[i] - string);
    
#if 0
    fprintf(stderr, "suffix array:\n");
    for (size_t i = 0; i < sa->length; i++) {
        fprintf(stderr, "sa[%d] == %4lu %s\n",
                i, sa->array[i], string + sa->array[i]);
    }
#endif
    
    return sa;
}

void compute_c_table(struct suffix_array *sa, const char *string)
{
    // I know we do not use all the characters, but this is easier
    // than preprocessing the string and matching it to a smaller set
    // and for the sizes of data I can handle, in any case, this won't
    // be the main problem.
    sa->c_table = (size_t*)calloc(C_TABLE_SIZE, sizeof(size_t));
    sa->c_table_symbols = (char*)calloc(C_TABLE_SIZE, 1);
    sa->c_table_symbols_inverse = (size_t*)calloc(C_TABLE_SIZE, sizeof(size_t));
    sa->c_table_no_symbols = 0;
    
    // first, count the occurrances of each symbol
    for (size_t i = 0; i < sa->length; i++) {
        size_t index = (size_t)string[i];
        sa->c_table[index]++;
    }
    
    // then, collect the non-empty bins
    for (int c = 0; c < C_TABLE_SIZE; c++) {
        if (sa->c_table[c] != 0) {
            sa->c_table_symbols[sa->c_table_no_symbols++] = (char)c;
        }
    }
    
    // get the reverse of that table...
    for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
        char symbol = sa->c_table_symbols[i];
        sa->c_table_symbols_inverse[(int)symbol] = i + 1;
    }
    
    // finally, adjust the table to the actual c-table
    size_t count = 0;
    // start loop at 1 to not count '$'
    sa->c_table['\0'] = 0;
    for (size_t i = 1; i < sa->c_table_no_symbols; i++) {
        char symbol = sa->c_table_symbols[i];
        size_t tmp = sa->c_table[(int)symbol];
        sa->c_table[(int)symbol] = count;
        count += tmp;
    }
}

void compute_o_table(struct suffix_array *sa, const char *string)
{
    // These must be computed first
    assert(sa->array);
    assert(sa->c_table);
    assert(sa->c_table_no_symbols);
    assert(sa->c_table_symbols);
    assert(sa->c_table_symbols_inverse);
    
    fprintf(stderr, "...building b table.\n");
    char *b = malloc(sa->length);
    for (size_t i = 0; i < sa->length; i++) {
        size_t sa_index = sa->array[i];
        if (sa_index == 0) {
            b[i] = '\0';
        } else {
            b[i] = string[sa_index - 1];
        }
    }
#if 0
    for (size_t i = 0; i < sa->length; i++) {
        fprintf(stderr, "b[%lu] == %c\n", i, b[i]);
    }
#endif
    
    size_t o_table_size = sa->c_table_no_symbols * sa->length;
    fprintf(stderr, "...allocating o-table size: [%lu x %lu] (%lu)\n",
            sa->c_table_no_symbols, sa->length,
            o_table_size);
    assert(o_table_size > 0);
    sa->o_table = malloc(o_table_size * sizeof(size_t));
    
    fprintf(stderr, "...building o-table.\n");
    // first column of the o-table is a special case because
    // we cannot look at b[i-1].
    size_t idx = 0;
    char b_symbol = b[0];
    for (size_t j = 0; j < sa->c_table_no_symbols; j++) {
        char row_symbol = sa->c_table_symbols[j];
        sa->o_table[idx] = (b_symbol == row_symbol) ? 1 : 0;
        idx += sa->length;
    }
    
    for (size_t i = 1; i < sa->length; i++) {
        b_symbol = b[i];
        for (size_t j = 0; j < sa->c_table_no_symbols; j++) {
            char row_symbol = sa->c_table_symbols[j];
            idx = o_table_index(sa, row_symbol, i);
            sa->o_table[idx] = (b_symbol == row_symbol) ?
                sa->o_table[idx-1] + 1 :
                sa->o_table[idx-1];
        }
    }
    free(b);
    fprintf(stderr, "...Done\n");

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

size_t o_table_index(struct suffix_array *sa, char symbol, size_t idx)
{
    // we are off by one so we can use 0
    // to indicate no index (easier with calloc)
    size_t symbol_idx = sa->c_table_symbols_inverse[(int)symbol];
    assert(symbol_idx > 0);
    size_t real_symbol_idx = symbol_idx - 1;
    return real_symbol_idx * sa->length + idx;
}


void delete_suffix_array(struct suffix_array *sa)
{
    if (sa->array)                   free(sa->array);
    
    if (sa->c_table)                 free(sa->c_table);
    if (sa->c_table_symbols)         free(sa->c_table_symbols);
    if (sa->c_table_symbols_inverse) free(sa->c_table_symbols_inverse);
    
    if (sa->o_table)                 free(sa->o_table);
    
    free(sa);
}

