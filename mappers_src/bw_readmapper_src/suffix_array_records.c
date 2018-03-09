
#include "suffix_array_records.h"
#include <stdlib.h>

struct suffix_array_records *empty_suffix_array_records()
{
    struct suffix_array_records *records =
        (struct suffix_array_records*)malloc(sizeof(struct suffix_array_records));
    records->names = empty_string_vector(10); // arbitrary size...
    records->suffix_arrays = empty_size_vector(10); // arbitrary size...
    return records;
}

void delete_suffix_array_records(struct suffix_array_records *records)
{
    delete_string_vector(records->names);
    delete_size_vector(records->suffix_arrays);
    free(records);
}

int read_suffix_array_records(struct suffix_array_records *records, FILE *file)
{
    // FIXME
    return 0;
}
