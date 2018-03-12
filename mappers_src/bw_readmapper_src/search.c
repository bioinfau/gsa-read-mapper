
#include "search.h"
#include "sam.h"

#include <strings.h>

// FIXME: only exact matching for now
void search(const char *read_name, const char *read, size_t read_idx,
            const char *quality,
            const char *ref_name, size_t L, size_t R, size_t d,
            struct suffix_array *sa,
            FILE *samfile)
{
#if 0
    printf("\"%s\" : \"%s\", L == %lu, R == %lu, %lu\n",
           read, read + read_idx, L, R, read_idx);
#endif
    
    // FIXME
    char cigar[1000];
    sprintf((char*)&cigar, "%luM", strlen(read));
    
    if (read_idx > 0) {
        // we have not reached the beginning of the read, so
        // update L and R and recurse
        char a = read[read_idx - 1];
        
        // quick test to abort if we see a symbol not in the reference
        if (sa->c_table_symbols_inverse[a] == 0)
            return; // no match (fixme when approximative)
        if (L == 0)
            L = sa->c_table[a] + 1;
        else
            L = sa->c_table[a] + 1 + sa->o_table[o_table_index(sa, a, L-1)];
        R = sa->c_table[a] + sa->o_table[o_table_index(sa, a, R)];
        return search(read_name, read, read_idx - 1, quality,
                      ref_name, L, R, d, sa, samfile);
        
    } else {
        assert(read_idx == 0);
        
        // we have matched to the end and can output
        // all sequences between L and R
        for (size_t i = L; i <= R; i++) {
            size_t index = sa->array[i];
            sam_line(samfile,
                     read_name,
                     ref_name,
                     index + 1, // + 1 for 1-indexing in SAM format.
                     cigar,
                     read,
                     quality);
        }
    }
    
}
