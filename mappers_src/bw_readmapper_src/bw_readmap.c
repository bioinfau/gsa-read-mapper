

/*
 Readmapper based on Burrows-Wheeler transform search.
 */

#include "suffix_array.h"
#include "fasta.h"
#include "fastq.h"
#include "sam.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>

struct search_info {
    int edit_dist;
    struct fasta_records *records;
    FILE *sam_file;
};

static struct search_info *empty_search_info()
{
    struct search_info *info =
        (struct search_info*)malloc(sizeof(struct search_info));
    info->edit_dist = 0;
    info->records = empty_fasta_records();
    return info;
}

static void delete_search_info(struct search_info *info)
{
    delete_fasta_records(info->records);
    free(info);
}

struct read_search_info {
    const char *ref_name;
    const char *read_name;
    const char *read;
    const char *quality;
    struct search_info *search_info;
};

static struct read_search_info *empty_read_search_info()
{
    struct read_search_info *info =
        (struct read_search_info*)malloc(sizeof(struct read_search_info));
    
    info->ref_name = 0;
    info->read_name = 0;
    info->read = 0;
    info->quality = 0;
    info->search_info = 0;

    return info;
}

static void delete_read_search_info(struct read_search_info *info)
{
    free(info);
}

static void read_callback(const char *read_name,
                          const char *read,
                          const char *quality,
                          void * callback_data) {
    struct search_info *search_info = (struct search_info*)callback_data;
    
    // I allocate and deallocate the info all the time... I might
    // be able to save some time by not doing this, but compared to
    // building and removeing the trie, I don't think it will be much.
    struct read_search_info *info = empty_read_search_info();
    info->search_info = search_info;
    info->read = read;
    info->quality = quality;
    info->read_name = read_name;
    
    // FIXME: SEARCH HERE
    
    delete_read_search_info(info);
}

int main(int argc, char * argv[])
{
    const char *prog_name = argv[0];
#if 0
    const char *algorithm = "naive";
    int opt; int edit_dist = 0;
    static struct option longopts[] = {
        { "help",       no_argument,            NULL,           'h' },
        { "distance",   required_argument,      NULL,           'd' },
        { "algorithm",  required_argument,      NULL,           'a' },
        { NULL,         0,                      NULL,            0  }
    };
    while ((opt = getopt_long(argc, argv, "hd:", longopts, NULL)) != -1) {
        switch (opt) {
            case 'h':
                printf("Usage: %s [options] ref.fa reads.fq\n\n", prog_name);
                printf("Options:\n");
                printf("\t-h | --help:\t\t Show this message.\n");
                printf("\t-d | --distance:\t Maximum edit distance for the search.\n");
                printf("\t-a | --algorithm:\t Algorithm to use for the search.\n");
                printf("\t\t\t\t Choices are:\n");
                printf("\t\t\t\t\t\"naive\"\n");
                printf("\t\t\t\t\t\"bmh\" (Boyer-Moore-Horspool)\n");
                printf("\t\t\t\t\t\"kmp\" (Knuth-Morris-Pratt)\n");
                printf("\t\t\t\t\t\"bsearch\" (suffix array binary search)\n");
                printf("\n\n");
                return EXIT_SUCCESS;
                
            case 'd':
                edit_dist = atoi(optarg);
                break;
                
            case 'a':
                algorithm = optarg;
                break;
                
            default:
                fprintf(stderr, "Usage: %s [options] ref.fa reads.fq\n", prog_name);
                return EXIT_FAILURE;
        }
    }
    argc -= optind;
    argv += optind;
    
    if (argc != 2) {
        fprintf(stderr, "Usage: %s [options] ref.fa reads.fq\n", prog_name);
        return EXIT_FAILURE;
    }
    
    FILE *fasta_file = fopen(argv[0], "r");
    if (!fasta_file) {
        fprintf(stderr, "Could not open %s.\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    FILE *fastq_file = fopen(argv[1], "r");
    if (!fastq_file) {
        fprintf(stderr, "Could not open %s.\n", argv[1]);
        return EXIT_FAILURE;
    }
    
    struct search_info *search_info = empty_search_info();
    search_info->edit_dist = edit_dist;
    
    if (strcmp(algorithm, "naive") == 0) {
        search_info->match_func = naive_exact_match;
    } else if (strcmp(algorithm, "bmh") == 0) {
        search_info->match_func = boyer_moore_horspool;
    } else if (strcmp(algorithm, "kmp") == 0) {
        search_info->match_func = knuth_morris_pratt;
    } else if (strcmp(algorithm, "bsearch") == 0) {
        search_info->match_func = suffix_array_bsearch_match;
    } else {
        fprintf(stderr, "Unknown search algorithm %s.\n", algorithm);
        return EXIT_FAILURE;
    }
    
    read_fasta_records(search_info->records, fasta_file);
    fclose(fasta_file);
    
    search_info->sam_file = stdout;
    
    scan_fastq(fastq_file, read_callback, search_info);
    delete_search_info(search_info);
    fclose(fastq_file);
    
#endif
    
    return EXIT_SUCCESS;
}
