
#include "cigar.h"
#include "sam.h"
#include "search.h"

#include <strings.h>

void search(const char *read_name, const char *read, size_t read_idx,
            const char *quality, const char *ref_name, size_t L, size_t R,
            int d, char *cigar, char *cigar_buffer, struct suffix_array *sa,
            FILE *samfile, struct options *options)
{
    assert(d >= 0); // if it get's negative we've called too deeply

    if (read_idx == 0) {
        // We have reached the beginning of the read.
        // Report what we have found

        // we have matched to the end and can output
        // all sequences between L and R
        simplify_cigar(cigar_buffer + 1, cigar);

        for (size_t i = L; i <= R; i++) {
            size_t index = sa->array[i];
            sam_line(samfile, read_name, ref_name,
                     index + 1, // + 1 for 1-indexing in SAM format.
                     cigar, read, quality);
        }

        // For completeness of the d-edit-cloud, we still need to
        // explore deletions...
        // ---DELETION----------------------------------------------
        if (d > 0) {
            size_t new_L, new_R;
            for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
                char b = sa->c_table_symbols[i];
                if (sa->c_table_symbols_inverse[(int)b] == 0)
                    return; // no match with this character

                if (L == 0)
                    new_L = sa->c_table[(int)b] + 1;
                else
                    new_L = sa->c_table[(int)b] + 1 +
                            sa->o_table[o_table_index(sa, b, L - 1)];
                new_R =
                    sa->c_table[(int)b] + sa->o_table[o_table_index(sa, b, R)];

                if (new_L > new_R)
                    continue;

                *cigar_buffer = 'D';
                search(read_name, read, read_idx, quality, ref_name, new_L,
                       new_R, d - 1, cigar, cigar_buffer - 1, sa, samfile, options);
            }
        }

        return; // all done.
    }

    // else: read_idx > 0
    // We have not reached the beginning of the read, so
    // update L and R and recurse

    // ---MATCHING----------------------------------------------
    // Get `a` as an exact match...
    char a = read[read_idx - 1];
    size_t new_L, new_R;
    if (L == 0)
        new_L = sa->c_table[(int)a] + 1;
    else
        new_L =
            sa->c_table[(int)a] + 1 + sa->o_table[o_table_index(sa, a, L - 1)];
    new_R = sa->c_table[(int)a] + sa->o_table[o_table_index(sa, a, R)];

    if (options->extended_cigars)
        *cigar_buffer = '=';
    else
        *cigar_buffer = 'M';

    search(read_name, read, read_idx - 1, quality, ref_name, new_L, new_R, d,
           cigar, cigar_buffer - 1, sa, samfile, options);

    if (d > 0) {
        // ---SUBSTITUTION------------------------------------------
        for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
            char b = sa->c_table_symbols[i];
            if (a == b || b == '\0')
                continue;

            if (sa->c_table_symbols_inverse[(int)b] == 0)
                continue; // no match with this character

            if (L == 0)
                new_L = sa->c_table[(int)b] + 1;
            else
                new_L = sa->c_table[(int)b] + 1 +
                        sa->o_table[o_table_index(sa, b, L - 1)];
            new_R = sa->c_table[(int)b] + sa->o_table[o_table_index(sa, b, R)];
            if (new_L > new_R)
                continue;

            if (options->extended_cigars)
                *cigar_buffer = 'X';
            else
                *cigar_buffer = 'M';

            search(read_name, read, read_idx - 1, quality, ref_name, new_L,
                   new_R, d - 1, cigar, cigar_buffer - 1, sa, samfile, options);
        } // end for

        // ---DELETION----------------------------------------------
        for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
            char b = sa->c_table_symbols[i];
            if (sa->c_table_symbols_inverse[(int)b] == 0)
                return; // no match with this character

            if (L == 0)
                new_L = sa->c_table[(int)b] + 1;
            else
                new_L = sa->c_table[(int)b] + 1 +
                        sa->o_table[o_table_index(sa, b, L - 1)];
            new_R = sa->c_table[(int)b] + sa->o_table[o_table_index(sa, b, R)];

            if (new_L > new_R)
                continue;

            *cigar_buffer = 'D';
            search(read_name, read, read_idx, quality, ref_name, new_L, new_R,
                   d - 1, cigar, cigar_buffer - 1, sa, samfile, options);
        } // end for

        // ---INSERTION---------------------------------------------
        *cigar_buffer = 'I';
        search(read_name, read, read_idx - 1, quality, ref_name, L, R, d - 1,
               cigar, cigar_buffer - 1, sa, samfile, options);
        
    } // end if (d > 0)
}
