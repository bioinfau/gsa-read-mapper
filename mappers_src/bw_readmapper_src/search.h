
#ifndef SEARCH_H
#define SEARCH_H

#include "suffix_array_records.h"

#include <stdlib.h>
#include <stdio.h>

void search(const char *read_name, const char *read, size_t read_idx,
            const char *quality,
            const char *ref_name, size_t L, size_t R, size_t d,
            struct suffix_array *sa,
            FILE *samfile);

#endif
