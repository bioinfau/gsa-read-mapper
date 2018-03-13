

#ifndef SIZE_VECTOR_H
#define SIZE_VECTOR_H

#include <stddef.h>

struct size_vector {
    size_t *sizes;
    size_t size;
    size_t used;
};

struct size_vector *empty_size_vector(size_t initial_size);
void delete_size_vector(struct size_vector *v);
struct size_vector *add_size(struct size_vector *v, size_t size);

#endif
