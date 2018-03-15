#include "trie.h"
#include "strings.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_LINE_SIZE 1024

int main(int argc, const char **argv)
{
	struct trie *trie = empty_trie();

	if (argc != 2) {
		fprintf(stderr, "Usage: %s patterns-file\n", argv[0]);
		return EXIT_FAILURE;
	}

	// read lines from input file and put them in the trie
	FILE *infile = fopen(argv[1], "r");

	if (!infile) {
		fprintf(stderr, "Could not open file %s\n", argv[1]);
		return EXIT_FAILURE;
	}

	size_t string_label = 0;
	char buffer[MAX_LINE_SIZE];
	while (fgets(buffer, MAX_LINE_SIZE, infile) != 0) {
		char pattern[MAX_LINE_SIZE], cigar[MAX_LINE_SIZE];
		sscanf(buffer, "%s %s", (char*)&pattern, (char*)&cigar);
        printf("adding \"%s\"\n", pattern);

        if (string_in_trie(trie, pattern)) {
        	printf("...already in trie.\n");
        } else {
        	add_string_to_trie(trie, pattern, (int)string_label++);
        }

	}
    fclose(infile);

    printf("computing failure links.\n");
    compute_failure_links(trie);

    printf("printing trie graph to \"trie.dot\"\n");
    print_dot(trie, "trie");

	return EXIT_SUCCESS;
}
