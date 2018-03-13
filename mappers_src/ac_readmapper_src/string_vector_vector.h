
#ifndef STRING_VECTOR_VECTOR_H
#define STRING_VECTOR_VECTOR_H

struct string_vector_vector {
    struct string_vector **string_vectors;
    int size;
    int used;
};

struct string_vector_vector *empty_string_vector_vector(int initial_size);
void delete_string_vector_vector(struct string_vector_vector *v);

// add a new string vector to the end of the vector vector and return its index
int append_vector(struct string_vector_vector *v);
// add a copy of a string to the vector at index
void add_string_copy_to_vector(struct string_vector_vector *v, int index, const char *s);

#endif
