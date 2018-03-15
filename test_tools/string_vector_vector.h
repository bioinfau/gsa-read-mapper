
#ifndef STRING_VECTOR_VECTOR_H
#define STRING_VECTOR_VECTOR_H

#include <stdlib.h>

struct string_vector_vector {
    struct string_vector **string_vectors;
    size_t size;
    size_t used;
};

struct string_vector_vector *empty_string_vector_vector(size_t initial_size);
void delete_string_vector_vector(struct string_vector_vector *v);

// add a new string vector to the end of the vector vector and return its index
size_t append_vector(struct string_vector_vector *v);
// add a copy of a string to the vector at index
void add_string_copy_to_vector(struct string_vector_vector *v, size_t index, const char *s);

#endif
