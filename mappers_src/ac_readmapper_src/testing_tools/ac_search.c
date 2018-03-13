#include "../trie.h"
#include "../strings.h"
#include "../string_vector.h"
#include "../string_vector_vector.h"
#include "../aho_corasick.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE_SIZE 1024

static char * read_database(FILE *file)
{
    char buffer[MAX_LINE_SIZE];
    size_t seq_size = MAX_LINE_SIZE;
    size_t n = 0;
    char *seq = malloc(seq_size);
    
    while (fgets(buffer, MAX_LINE_SIZE, file) != 0) {
        
        for (char *c = buffer; *c; ++c) {
            if (!isalpha(*c)) continue;
            
            seq[n++] = *c;
            
            if (n == seq_size) {
                seq_size *= 2;
                seq = (char*)realloc(seq, seq_size);
            }
        }
    }
    
	return seq;
}

struct callback_data {
	struct string_vector *patterns;
	struct string_vector_vector *cigars;
	const char *database;
};

static void match_callback(int string_label, size_t index, void * data)
{
	struct callback_data *cb_data = (struct callback_data*)data;
	const char *str = cb_data->patterns->strings[string_label];
	size_t n = strlen(str);
	printf("string \"%s\" [%d] matches at index %lu [%s].\n", 
		str, string_label, index - n + 1, cb_data->database + index - n + 1);
	struct string_vector *cigars = cb_data->cigars->string_vectors[index];
	printf("\t%d matches:\n", cigars->used);
	for (size_t i = 0; i < cigars->used; i++) {
		printf("\t\t%s\n", cigars->strings[i]);
	}
	
}

int main(int argc, const char **argv)
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s database-file patterns-file\n", argv[0]);
		return EXIT_FAILURE;
	}

	FILE *database_file = fopen(argv[1], "r");
	if (!database_file) {
		fprintf(stderr, "Could not open file %s\n", argv[1]);
		return EXIT_FAILURE;
	}
	
	FILE *patterns_file = fopen(argv[2], "r");
	if (!patterns_file) {
		fprintf(stderr, "Could not open file %s\n", argv[2]);
		return EXIT_FAILURE;
	}
	
	char *database = read_database(database_file);
	fclose(database_file);

	struct string_vector *patterns = empty_string_vector(10);
	struct string_vector_vector *cigars = empty_string_vector_vector(10);
	struct trie *trie = empty_trie();

	char buffer[MAX_LINE_SIZE];
	while (fgets(buffer, MAX_LINE_SIZE, patterns_file) != 0) {
		char pattern[MAX_LINE_SIZE], cigar[MAX_LINE_SIZE];
		sscanf(buffer, "%s %s", (char*)&pattern, (char*)&cigar);
        printf("adding \"%s\" [%s]\n", pattern, cigar);

        if (string_in_trie(trie, pattern)) {
        	// The pattern is already in the tree, but if we are called here
        	// we have a new CIGAR for the same pattern.
        	struct trie *node = get_trie_node(trie, pattern);
        	add_string_copy_to_vector(cigars, node->string_label, cigar);

    	} else {
        	// NB: the order is important here -- info->patterns->used will be updated
        	// when we add the pattern to the vector, so we insert in the trie first.
        	int index = append_vector(cigars);
        	add_string_to_trie(trie, pattern, index);
        	add_string_copy(patterns, pattern);
        	add_string_copy_to_vector(cigars, index, cigar);
    	}
	}
    fclose(patterns_file);

    printf("computing failure links.\n");
    compute_failure_links(trie);

    printf("searching...\n");
    struct callback_data cb_data;
    cb_data.patterns = patterns;
    cb_data.cigars = cigars;
    cb_data.database = database;
    aho_corasick_match(database, strlen(database), trie, match_callback, &cb_data);

    free(database);
    delete_trie(trie);
    delete_string_vector(patterns);

	return EXIT_SUCCESS;
}
