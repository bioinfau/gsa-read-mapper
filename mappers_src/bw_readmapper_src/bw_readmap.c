

/*
 Readmapper based on Burrows-Wheeler transform search.
 */

#include "suffix_array.h"
#include "suffix_array_records.h"
#include "fasta.h"
#include "fastq.h"
#include "sam.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>

struct search_info {
    int edit_dist;
    struct fasta_records *fasta_records;
    struct suffix_array_records *sa_records;
    FILE *sam_file;
};

static struct search_info *empty_search_info()
{
    struct search_info *info =
        (struct search_info*)malloc(sizeof(struct search_info));
    info->edit_dist = 0;
    info->fasta_records = empty_fasta_records();
    info->sa_records = empty_suffix_array_records();
    return info;
}

static void delete_search_info(struct search_info *info)
{
    delete_fasta_records(info->fasta_records);
    delete_suffix_array_records(info->sa_records);
    free(info);
}

static void read_callback(const char *read_name,
                          const char *read,
                          const char *quality,
                          void * callback_data) {
    struct search_info *search_info = (struct search_info*)callback_data;
    
    // FIXME: SEARCH HERE
    
}

int main(int argc, char * argv[])
{
    const char *prog_name = argv[0];

    int opt; bool preprocess = false; int edit_dist = 0;
    static struct option longopts[] = {
        { "help",        no_argument,            NULL,           'h' },
        { "preprocess",  required_argument,      NULL,           'p' },
        { "distance",    required_argument,      NULL,           'd' },
        { NULL,         0,                      NULL,            0  }
    };
    while ((opt = getopt_long(argc, argv, "hpd:", longopts, NULL)) != -1) {
        switch (opt) {
            case 'h':
                printf("Usage: %s -p | --preprocess ref.fa\n"
                       "       %s [-d distance] ref.fa reads.fq\n",
                       prog_name, prog_name);
                printf("Options:\n");
                printf("\t-h | --help:\t\t Show this message.\n");
                printf("\t-d | --distance:\t Maximum edit distance for the search.\n");
                printf("\t-p | --preprocess:\t Preprocess a reference genome.\n");
                printf("\n\n");
                return EXIT_SUCCESS;
                
            case 'p':
                preprocess = true;
                break;

            case 'd':
                edit_dist = atoi(optarg);
                break;
                
            default:
                fprintf(stderr,
                        "Usage: %s -p | --preprocess ref.fa\n"
                        "       %s [-d distance] ref.fa reads.fq\n",
                        prog_name, prog_name);
                return EXIT_FAILURE;
        }
    }
    argc -= optind;
    argv += optind;

    if (preprocess) {
        if (argc != 1) {
            fprintf(stderr,
                    "Usage: %s -p | --preprocess ref.fa\n"
                    "       %s [-d distance] ref.fa reads.fq\n",
                    prog_name, prog_name);
            return EXIT_FAILURE;
        }
        
        FILE *fasta_file = fopen(argv[0], "r");
        if (!fasta_file) {
            fprintf(stderr, "Could not open %s.\n", argv[0]);
            return EXIT_FAILURE;
        }
        
        struct fasta_records *records = empty_fasta_records();
        if (0 != read_fasta_records(records, fasta_file)) {
            fprintf(stderr, "Could not read FASTA file.\n");
            return EXIT_FAILURE;
        }
        fclose(fasta_file);
        
        struct suffix_array_records *sa_records =
            build_suffix_array_records(records);
        
        write_suffix_array_records(sa_records, records, argv[0]);
        
        delete_suffix_array_records(sa_records);
        delete_fasta_records(records);
        
    } else {
        if (argc != 2) {
            fprintf(stderr,
                    "Usage: %s -p | --preprocess ref.fa\n"
                    "       %s [-d distance] ref.fa reads.fq\n",
                    prog_name, prog_name);
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
        
        if (0 != read_fasta_records(search_info->fasta_records, fasta_file)) {
            fprintf(stderr, "Could not read FASTA file.\n");
            delete_search_info(search_info);
            return EXIT_FAILURE;
        }
        fclose(fasta_file);
        
        if (0 != read_suffix_array_records(search_info->sa_records,
                                           search_info->fasta_records,
                                           argv[0])) {
            fprintf(stderr, "Could not read suffix arrays.\n");
            delete_search_info(search_info);
            return EXIT_FAILURE;
        }
        
        write_suffix_array_records(
            search_info->sa_records, search_info->fasta_records, "test"
        );
        
        search_info->sam_file = stdout;
        scan_fastq(fastq_file, read_callback, search_info);

        delete_search_info(search_info);
        fclose(fastq_file);

    }
        
    return EXIT_SUCCESS;
}
