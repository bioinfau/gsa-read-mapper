

/*
 Readmapper based on Burrows-Wheeler transform search.
 */

#include "suffix_array.h"
#include "suffix_array_records.h"
#include "fasta.h"
#include "fastq.h"
#include "sam.h"
#include "search.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>

struct search_info {
    int edit_dist;
    struct fasta_records *fasta_records;
    struct suffix_array_records *sa_records;
    int max_edit_distance;
    FILE *samfile;
};

static struct search_info *empty_search_info(int max_edit_distance)
{
    struct search_info *info =
        (struct search_info*)malloc(sizeof(struct search_info));
    info->edit_dist = 0;
    info->fasta_records = empty_fasta_records();
    info->sa_records = empty_suffix_array_records();
    info->max_edit_distance = max_edit_distance;
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
    struct search_info *info = (struct search_info*)callback_data;
    
    size_t n = strlen(read) + (size_t)info->max_edit_distance;
    char cigar[n + 1], cigar_buffer[n + 1];
    cigar[n] = cigar_buffer[n] = '\0';
    
    size_t no_records = info->fasta_records->names->used;
    for (size_t seq_no = 0; seq_no < no_records; seq_no++) {
        char *ref_name = info->fasta_records->names->strings[seq_no];
        struct suffix_array *sa = info->sa_records->suffix_arrays[seq_no];
        
        search(read_name, read, strlen(read),
               quality,
               ref_name, 0, sa->length - 1, info->edit_dist,
               cigar, cigar_buffer + n - 1,
               sa, info->samfile);
        
    }
}

int main(int argc, char * argv[])
{
    const char *prog_name = argv[0];

    int opt;
    bool preprocess = false;
    bool validate = false;
    int edit_dist = 0;
    static struct option longopts[] = {
        { "help",        no_argument,            NULL,           'h' },
        { "preprocess",  no_argument,            NULL,           'p' },
        { "validate",    no_argument,            NULL,           'v' },
        { "distance",    required_argument,      NULL,           'd' },
        { NULL,         0,                       NULL,            0  }
    };
    while ((opt = getopt_long(argc, argv, "hpvd:", longopts, NULL)) != -1) {
        switch (opt) {
            case 'h':
                printf("Usage: %s -p | --preprocess ref.fa\n"
                       "       %s -v | --validate ref.fa\n"
                       "       %s [-d distance] ref.fa reads.fq\n",
                       prog_name, prog_name, prog_name);
                printf("Options:\n");
                printf("\t-h | --help:\t\t Show this message.\n");
                printf("\t-d | --distance:\t Maximum edit distance for the search.\n");
                printf("\t-p | --preprocess:\t Preprocess a reference genome.\n");
                printf("\t-v | --validate:\t Checks that preprocessed data matches the computed tables reference genome.\n");
                printf("\n\n");
                return EXIT_SUCCESS;
                
            case 'p':
                preprocess = true;
                break;

            case 'v':
                validate = true;
                break;
                
            case 'd':
                edit_dist = atoi(optarg);
                break;
                
            default:
                fprintf(stderr,
                        "Usage: %s -p | --preprocess ref.fa\n"
                        "       %s -v | --validate ref.fa\n"
                        "       %s [-d distance] ref.fa reads.fq\n",
                        prog_name, prog_name, prog_name);
                return EXIT_FAILURE;
        }
    }
    argc -= optind;
    argv += optind;

    if (preprocess) {
        if (argc != 1) {
            fprintf(stderr,
                    "Usage: %s -p | --preprocess ref.fa\n"
                    "       %s -v | --validate ref.fa\n"
                    "       %s [-d distance] ref.fa reads.fq\n",
                    prog_name, prog_name, prog_name);
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
        
#if 0
        fprintf(stderr, "suffix array structures after built.\n");
        for (int i = 0; i < records->names->used; i++) {
            struct suffix_array *sa = sa_records->suffix_arrays[i];
            const char *string = records->sequences->strings[i];
          
            fprintf(stderr, "suffix array:\n");
            for (int i = 0; i < sa->length; i++) {
                fprintf(stderr, "sa[%d] == %4lu %s\n",
                        i, sa->array[i], string + sa->array[i]);
            }

            for (size_t i = 0; i < sa->c_table_no_symbols; i++) {
                char symbol = sa->c_table_symbols[i];
                printf("O(%c,) =", (symbol == 0) ? '$' : symbol);
                for (size_t j = 0; j < sa->length; ++j) {
                    size_t idx = o_table_index(sa, symbol, j);
                    printf(" %lu", sa->o_table[idx]);
                }
                printf("\n");
            }

        }
#endif
        
        write_suffix_array_records(sa_records, records, argv[0]);
        
        fprintf(stderr, "Clean up suffic array records.\n");
        delete_suffix_array_records(sa_records);
        fprintf(stderr, "Clean up FASTA records.\n");
        delete_fasta_records(records);
        fprintf(stderr, "All done.\n");
       
    } else if (validate) {
        
        if (argc != 1) {
            fprintf(stderr,
                    "Usage: %s -p | --preprocess ref.fa\n"
                    "       %s -v | --validate ref.fa\n"
                    "       %s [-d distance] ref.fa reads.fq\n",
                    prog_name, prog_name, prog_name);
            return EXIT_FAILURE;
        }
        
        FILE *fasta_file = fopen(argv[0], "r");
        if (!fasta_file) {
            fprintf(stderr, "Could not open %s.\n", argv[0]);
            return EXIT_FAILURE;
        }
        
        struct fasta_records *fasta_records = empty_fasta_records();
        if (0 != read_fasta_records(fasta_records, fasta_file)) {
            fprintf(stderr, "Could not read FASTA file.\n");
            return EXIT_FAILURE;
        }
        fclose(fasta_file);
        
        struct suffix_array_records *built_sa_records =
            build_suffix_array_records(fasta_records);
        
        struct suffix_array_records *loaded_sa_records = empty_suffix_array_records();
        
        if (0 != read_suffix_array_records(loaded_sa_records,
                                           fasta_records,
                                           argv[0])) {
            fprintf(stderr, "Could not read suffix arrays.\n");
            delete_fasta_records(fasta_records);
            delete_suffix_array_records(built_sa_records);
            delete_suffix_array_records(loaded_sa_records);
            return EXIT_FAILURE;
        }
        
        fprintf(stdout, "Checking that suffix arrays are the same.\n");
        size_t no_records = fasta_records->names->used;
        for (size_t i = 0; i < no_records; ++i) {
            const char *seq_name = fasta_records->names->strings[i];
            fprintf(stdout, "...comparing for %s.\n", seq_name);
            struct suffix_array *built_sa = built_sa_records->suffix_arrays[i];
            struct suffix_array *loaded_sa = loaded_sa_records->suffix_arrays[i];
            if (built_sa->length != loaded_sa->length) {
                fprintf(stderr, "The suffix arrays disagree on length. %lu != %lu\n",
                        built_sa->length, loaded_sa->length);
                return EXIT_FAILURE;
            }
            for (size_t j = 0; j < built_sa->length; j++) {
                if (built_sa->array[j] != loaded_sa->array[j]) {
                    fprintf(stderr, "The suffix arrays disagree on index %lu : %lu != %lu\n",
                            j, built_sa->array[j], loaded_sa->array[j]);
                    return EXIT_FAILURE;
                }
            }
        }
        
        delete_fasta_records(fasta_records);
        delete_suffix_array_records(built_sa_records);
        delete_suffix_array_records(loaded_sa_records);
        return EXIT_SUCCESS;
        
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
    
        struct search_info *search_info = empty_search_info(edit_dist);
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
        
        search_info->samfile = stdout;
        scan_fastq(fastq_file, read_callback, search_info);

        delete_search_info(search_info);
        fclose(fastq_file);

    }
        
    return EXIT_SUCCESS;
}
