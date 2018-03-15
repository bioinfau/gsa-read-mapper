#include <stdlib.h>
#include <stdio.h>
#include "edit_distance_generator.h"

static void callback(const char *string, const char *cigar, void * data)
{
	printf("%s %s\n", string, cigar);
}

int main(int argc, const char **argv)
{
	if (argc != 4) {
		fprintf(stderr, "Usage: %s distance alphabet string\n", argv[0]);
		return EXIT_FAILURE;
	}

	generate_all_neighbours(argv[3], argv[2], atoi(argv[1]), callback, 0);

	return EXIT_SUCCESS;
}
