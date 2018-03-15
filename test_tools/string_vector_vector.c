
#include "string_vector.h"
#include "string_vector_vector.h"
#include <assert.h>
#include <stdlib.h>

struct string_vector_vector* empty_string_vector_vector(size_t initial_size)
{
    struct string_vector_vector* v = malloc(sizeof(struct string_vector_vector));
    v->string_vectors = malloc(sizeof(struct string_vector*) * initial_size);
    v->size = initial_size;
    v->used = 0;
    return v;
}

void delete_string_vector_vector(struct string_vector_vector* v)
{
    for (size_t i = 0; i < v->used; ++i)
    {
        delete_string_vector(v->string_vectors[i]);
    }
    free(v->string_vectors);
    free(v);
}

// add a new string vector to the end of the vector vector and return its index
size_t append_vector(struct string_vector_vector* v)
{
    if (v->size == v->used)
    {
        v->string_vectors =
            realloc(v->string_vectors, 2 * v->size * sizeof(struct string_vector*));
        v->size = 2 * v->size;
    }
    v->string_vectors[v->used++] = empty_string_vector(1);
    return v->used - 1;
}

// add a copy of a string to the vector at index
void add_string_copy_to_vector(struct string_vector_vector* v, size_t index, const char* s)
{
    assert(index <= v->used);
    add_string_copy(v->string_vectors[index], s);
}
