
#include "suffix_array.h"
#include "strings.h"
#include "pair_stack.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>


static struct suffix_array *allocate_sa(char *string)
{
    struct suffix_array *sa =
        (struct suffix_array*)malloc(sizeof(struct suffix_array));
    sa->string = string;
    sa->length = strlen(string);
    sa->array = (size_t*)malloc(sa->length * sizeof(size_t));
    
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


void delete_suffix_array(struct suffix_array *sa)
{
    if(sa->string) free(sa->string);
    free(sa->array);
    free(sa);
}

// when searching, we cannot simply use bsearch because we want
// to get a lower bound if the key isn't in the array -- bsearch
// would give us NULL in that case.
size_t lower_bound_search(struct suffix_array *sa, const char *key)
{
    int low = 0;
    int high = (int)sa->length;
    int mid = 0;
    int cmp = 0;
    size_t key_len = strlen(key);
    
    while (low < high) {
        mid = low + (high-low) / 2;
        cmp = strncmp(key, sa->string + sa->array[mid], key_len);
        if (cmp < 0) {
            high = mid - 1;
        } else if (cmp > 0) {
            low = mid + 1;
        } else {
            // a hit, search down until we get the smallest hit...
            for (int i = mid - 1; i >= 0; --i) {
                if (strncmp(sa->string + sa->array[i], key, key_len) < 0)
                    return i + 1;
            }
            return 0; // if we get here, we didn't find a smaller string, so we
                      // have to return 0
        }
    }
    
    // we didn't find the key -- we are either at the smallest upper bound
    // or highest lower bound. The relative order of mid and high tells us which
    assert(cmp != 0);
    //printf("fell through ... %d %d %d\n", low, mid, high);
    if (high < 0) return 0;
    if (low >= sa->length) return sa->length - 1;
    
    if (high < mid) {
        // we moved down so mid points to a larger string. This means
        // that the largest that is smaller must be one below mid (which is high)
        return mid - 1;
    } else {
        // we moved up, so mid is smaller. we have two possible cases, then
        // either low points to a match or mid is the largest that is smaller
        assert(low > mid);
        if (strncmp(sa->string + sa->array[low], key, key_len) == 0)
            return low;
        else
            return mid;
    }
    assert(false); // we should never get here.
}

void suffix_array_bsearch_match(const char *text, size_t n,
                                const char *pattern, size_t m,
                                match_callback_func callback,
                                void *callback_data)
{
    struct suffix_array *sa = qsort_sa_construction(string_copy(text));
    size_t lb = lower_bound_search(sa, pattern);
    for (size_t i = lb; i < sa->length; ++i) {
        if (strncmp(pattern, sa->string + sa->array[i], m) != 0)
             break;
        callback(sa->array[i], callback_data);
    }
    delete_suffix_array(sa);
}
