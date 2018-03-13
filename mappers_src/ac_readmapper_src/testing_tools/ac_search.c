#include "../trie.h"
#include "../strings.h"
#include "../string_vector.h"
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
	const char *database;
};

static void match_callback(int string_label, size_t index, void * data)
{
	struct callback_data *cb_data = (struct callback_data*)data;
	const char *str = cb_data->patterns->strings[string_label];
	size_t n = strlen(str);
	printf("string \"%s\" [%d] matches at index %lu [\"%s\"].\n", 
		str, string_label, index - n + 1, cb_data->database + index - n + 1);
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
	struct trie *trie = empty_trie();
	char buffer[MAX_LINE_SIZE];
	int string_label = 0;
	while (fgets(buffer, MAX_LINE_SIZE, patterns_file) != 0) {
		char *str = string_copy(strtok(buffer, "\n"));
        printf("adding \"%s\"\n", str);
        patterns = add_string_copy(patterns, str);
        add_string_to_trie(trie, str, string_label++);
	}
    fclose(patterns_file);

    printf("computing failure links.\n");
    compute_failure_links(trie);

    printf("searching...\n");
    struct callback_data cb_data;
    cb_data.patterns = patterns;
    cb_data.database = database;
    aho_corasick_match(database, strlen(database), trie, match_callback, &cb_data);

    free(database);
    delete_trie(trie);
    delete_string_vector(patterns);

	return EXIT_SUCCESS;
}
