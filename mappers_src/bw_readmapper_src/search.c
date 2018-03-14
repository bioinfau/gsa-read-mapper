
#include "search.h"
#include "sam.h"
#include "cigar.h"

#include <strings.h>


void search(const char *read_name, const char *read, size_t read_idx,
            const char *quality,
            const char *ref_name, size_t L, size_t R, int d,
            char *cigar, char *cigar_buffer,
            char *match_buffer, // fixme: debug
            struct suffix_array *sa,
            FILE *samfile)
{
#if 1
    fprintf(stderr, "\"%s\" : \"%s\" [%s;%s](len=%lu), read_idx == %lu, d == %d, L == %lu, R == %lu, %lu\n",
            read, read + read_idx,
            cigar_buffer + 1, match_buffer + 1,
            strlen(cigar_buffer + 1),
            read_idx, d, L, R, read_idx);
#endif
    
    assert(d >= 0); // if it get's negative we've called too deeply
    
    if (read_idx == 0) {
        // We have reached the beginning of the read.
        // Report what we have found
        
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
        return; // all done.
    }
    
    
    // We have not reached the beginning of the read, so
    // update L and R and recurse
    
    // ---MATCHING----------------------------------------------
    // Get `a` as an exact match...
    char a = read[read_idx - 1];
    size_t new_L, new_R;
    if (L == 0)
        new_L = sa->c_table[(int)a] + 1;
    else
        new_L = sa->c_table[(int)a] + 1 + sa->o_table[o_table_index(sa, a, L-1)];
    new_R = sa->c_table[(int)a] + sa->o_table[o_table_index(sa, a, R)];
    
    fprintf(stderr, "Attempting to match %c, %lu [%lu], %lu [%lu]\n",
            a,
            new_L, sa->array[new_L],
            new_R, sa->array[new_R]);

#ifdef EXTENDED_CIGAR
    *cigar_buffer = '=';
#else
    *cigar_buffer = 'M';
#endif
    *match_buffer = a; // FIXME
    
    search(read_name, read, read_idx - 1, quality,
           ref_name, new_L, new_R, d,
           cigar, cigar_buffer - 1,
           match_buffer - 1,
           sa, samfile);
    
    if (d > 0) {
        // ---SUBSTITUTION------------------------------------------
        for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
            char b = sa->c_table_symbols[i];
            if (a == b || b == '\0') continue;
            
            if (sa->c_table_symbols_inverse[(int)b] == 0)
                continue; // no match with this character
            
            if (L == 0)
                new_L = sa->c_table[(int)b] + 1;
            else
                new_L = sa->c_table[(int)b] + 1 + sa->o_table[o_table_index(sa, b, L-1)];
            new_R = sa->c_table[(int)b] + sa->o_table[o_table_index(sa, b, R)];
            
            fprintf(stderr, "Attempting substitution (%c -> %c), %lu [%lu], %lu [%lu]\n",
                    a, b,
                    new_L, sa->array[new_L],
                    new_R, sa->array[new_R]);
            if (new_L > new_R) continue;
            
#ifdef EXTENDED_CIGAR
            *cigar_buffer = 'X';
#else
            *cigar_buffer = 'M';
#endif
            *match_buffer = b; // FIXME
            search(read_name, read, read_idx - 1, quality,
                   ref_name, new_L, new_R, d - 1,
                   cigar, cigar_buffer - 1,
                   match_buffer - 1,
                   sa, samfile);
        }
        
        // ---DELETION----------------------------------------------
        for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
            char b = sa->c_table_symbols[i];
            if (sa->c_table_symbols_inverse[(int)b] == 0)
                return; // no match with this character
            
            if (L == 0)
                new_L = sa->c_table[(int)b] + 1;
            else
                new_L = sa->c_table[(int)b] + 1 + sa->o_table[o_table_index(sa, b, L-1)];
            new_R = sa->c_table[(int)b] + sa->o_table[o_table_index(sa, b, R)];
            
            if (new_L > new_R) continue;
            
            *cigar_buffer = 'D';
            *match_buffer = b; // FIXME
            search(read_name, read, read_idx, quality,
                   ref_name, new_L, new_R, d - 1,
                   cigar, cigar_buffer - 1,
                   match_buffer - 1,
                   sa, samfile);
        }
        
        // ---INSERTION---------------------------------------------
        *cigar_buffer = 'I';
        search(read_name, read, read_idx - 1, quality,
               ref_name, L, R, d - 1,
               cigar, cigar_buffer - 1,
               match_buffer, // FIXME
               sa, samfile);

    }
}

