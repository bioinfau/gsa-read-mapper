
#include "search.h"
#include "sam.h"
#include "cigar.h"

#include <strings.h>

// FIXME: only exact matching for now
void search(const char *read_name, const char *read, size_t read_idx,
            const char *quality,
            const char *ref_name, size_t L, size_t R, int d,
            char *cigar, char *cigar_buffer,
            struct suffix_array *sa,
            FILE *samfile)
{
#if 1
    printf("\"%s\" : \"%s\", d == %d, L == %lu, R == %lu, %lu\n",
           read, read + read_idx, d, L, R, read_idx);
    printf("cigar buffer = \"%s\".\n", cigar_buffer + 1);
#endif
    
    if (read_idx > 0) {
        // we have not reached the beginning of the read, so
        // update L and R and recurse
        
        char a = read[read_idx - 1];
        size_t new_L, new_R;
        
        for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
            char b = sa->c_table_symbols[i];
            if (sa->c_table_symbols_inverse[b] == 0)
                return; // no match with this character
                
            if (L == 0)
                new_L = sa->c_table[b] + 1;
            else
                new_L = sa->c_table[b] + 1 + sa->o_table[o_table_index(sa, b, L-1)];
            new_R = sa->c_table[b] + sa->o_table[o_table_index(sa, b, R)];
            
            // matches/mismatches
            if (a == b) {
                *cigar_buffer = 'M';
                search(read_name, read, read_idx - 1, quality,
                       ref_name, new_L, new_R, d,
                       cigar, cigar_buffer - 1,
                       sa, samfile);
            }
            if (d > 0 && a != b) {
                *cigar_buffer = 'M';
                search(read_name, read, read_idx - 1, quality,
                       ref_name, new_L, new_R, d - 1,
                       cigar, cigar_buffer - 1,
                       sa, samfile);
            }
            
            // insertion -- use L and R adjusted with 'b' but do not decrease index
            if (d > 0) {
                *cigar_buffer = 'I';
                search(read_name, read, read_idx, quality,
                       ref_name, new_L, new_R, d - 1,
                       cigar, cigar_buffer - 1,
                       sa, samfile);
            }
        }

        if (d > 0) {
            // deletion -- just skip a character but continue with same L and R
            *cigar_buffer = 'D';
            search(read_name, read, read_idx - 1, quality,
                   ref_name, L, R, d - 1,
                   cigar, cigar_buffer - 1,
                   sa, samfile);
        }
        
    } else {
        // we have reached the beginning of the read -- report what we have found
        assert(read_idx == 0);
        
        // we have matched to the end and can output
        // all sequences between L and R
        simplify_cigar(cigar_buffer + 1, cigar);
        
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
