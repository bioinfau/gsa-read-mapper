#include "../trie.h"
#include "../strings.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_LINE_SIZE 1024

int main(int argc, const char **argv)
{
	struct trie *trie = empty_trie();

	if (argc != 2) {
		fprintf(stderr, "Usage: %s input-file\n", argv[0]);
		return EXIT_FAILURE;
	}

	// read lines from input file and put them in the trie
	FILE *infile = fopen(argv[1], "r");
	
	if (!infile) {
		fprintf(stderr, "Could not open file %s\n", argv[1]);
		return EXIT_FAILURE;
	}

	char buffer[MAX_LINE_SIZE];
	int string_label = 0;
	while (fgets(buffer, MAX_LINE_SIZE, infile) != 0) {
		char *str = string_copy(strtok(buffer, "\n"));
        printf("adding \"%s\"\n", str);
        add_string_to_trie(trie, str, string_label++);
	}
    fclose(infile);

    printf("computing failure links.\n");
    compute_failure_links(trie);

    printf("printing trie graph to \"trie.dot\"\n");
    print_dot(trie, "trie");

	return EXIT_SUCCESS;
}