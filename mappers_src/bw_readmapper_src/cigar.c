
#include "cigar.h"

#include <string.h>
#include <stdio.h>

static const char *scan(const char *cigar)
{
    const char *p = cigar;
    while (*p == *cigar)
        ++p;
    return p;
}

void simplify_cigar(const char *from, char *to)
{
    while (*from) {
        const char *next = scan(from);
        to = to + sprintf(to, "%lu%c", next - from, *from);
        from = next;
    }
    *to = '\0';
}
