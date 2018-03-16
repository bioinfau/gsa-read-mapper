#include "edit_distance_generator.h"
#include "options.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

static void callback(const char *string, const char *cigar, void *data) {
    printf("%s %s\n", string, cigar);
}

int main(int argc, char *argv[]) {
    struct options options;
    options.edit_distance = 0;
    options.extended_cigars = false;
    options.verbose = false;
    
    const char *progname = argv[0];
    int opt;
    static struct option longopts[] = {
        { "help",             no_argument,            NULL,           'h' },
        { "distance",         required_argument,      NULL,           'd' },
        { "extended-cigar",   no_argument,            NULL,           'x' },
        { "verbose",          no_argument,            NULL,           'v' },
        { NULL,               0,                      NULL,            0  }
    };
    while ((opt = getopt_long(argc, argv, "hvd:x", longopts, NULL)) != -1) {
        switch (opt) {
            case 'h':
                printf("Usage: %s [options] alphabet string\n\n", argv[0]);
                printf("Options:\n");
                printf("\t-h | --help:\t\t Show this message.\n");
                printf("\t-d | --distance:\t Maximum edit distance for the search.\n");
                printf("\t-v | --verbose:\t Verbose output.\n");
                printf("\t-x | --extended-cigar:\t Use extended CIGAR format in SAM output.\n");
                printf("\n\n");
                return EXIT_SUCCESS;
                
            case 'd':
                options.edit_distance = atoi(optarg);
                break;
                
            case 'x':
                options.extended_cigars = true;
                break;
                
            default:
                fprintf(stderr, "Usage: %s [options] ref.fa reads.fq\n", argv[0]);
                return EXIT_FAILURE;
        }
    }
    
    argc -= optind;
    argv += optind;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s [options] alphabet string\n", argv[0]);
        return EXIT_FAILURE;
    }

    generate_all_neighbours(argv[3], argv[2], atoi(argv[1]), callback, 0, &options);

    return EXIT_SUCCESS;
}
