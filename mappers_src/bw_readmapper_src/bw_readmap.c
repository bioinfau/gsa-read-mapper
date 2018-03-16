

/*
 Readmapper based on Burrows-Wheeler transform search.
 */

#include "fasta.h"
#include "fastq.h"
#include "sam.h"
#include "search.h"
#include "suffix_array.h"
#include "suffix_array_records.h"
#include "options.h"

#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(const char *prog_name, FILE *file)
{
    fprintf(file, "Usage: %s -p | --preprocess ref.fa\n"
            "              %s [search options] ref.fa reads.fq\n",
            prog_name, prog_name);
    
    fprintf(file, "Options:\n");
    fprintf(file, "\t-h | --help:\t\t Show this message.\n");
    fprintf(file, "\t-p | --preprocess:\t Preprocess a reference genome.\n");
    fprintf(file, "\nSearch options:\n");
    fprintf(file, "\t-d | --distance:\t Maximum edit distance for the search.\n");
    fprintf(file, "\t-x | --extended-cigar:\t Use extended CIGAR notation in SAM output.\n");
    fprintf(file, "\n\n");
}

#define FASTQ_BUFFER_SIZE 1024

int main(int argc, char *argv[]) {
    
    int opt;
    struct options options;
    options.verbose = false;
    options.extended_cigars = false;
    options.edit_distance = 0;
    bool preprocess = false;
    
    static struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"preprocess", no_argument, NULL, 'p'},
        {"distance", required_argument, NULL, 'd'},
        {"extended-cigar", no_argument, NULL, 'x'},
        {NULL, 0, NULL, 0}};
    while ((opt = getopt_long(argc, argv, "hpd:x", longopts, NULL)) != -1) {
        switch (opt) {
            case 'h':
                print_usage(argv[0], stdout);
                return EXIT_SUCCESS;
                
            case 'p':
                preprocess = true;
                break;
                
            case 'd':
                options.edit_distance = atoi(optarg);
                break;
                
            case 'x':
                options.extended_cigars = true;
                break;
                
            default:
                print_usage(argv[0], stderr);
                return EXIT_FAILURE;
        }
    }
    argc -= optind;
    argv += optind;
    
    if (preprocess) {
        if (argc != 1) {
            print_usage(argv[0], stderr);
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
            print_usage(argv[0], stderr);
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
        
        struct fasta_records *fasta_records = empty_fasta_records();
        if (0 != read_fasta_records(fasta_records, fasta_file)) {
            fprintf(stderr, "Could not read FASTA file.\n");
            return EXIT_FAILURE;
        }
        fclose(fasta_file);
        
        struct suffix_array_records *sa_records = empty_suffix_array_records();
        if (0 != read_suffix_array_records(sa_records, fasta_records, argv[0])) {
            fprintf(stderr, "Could not read suffix arrays.\n");
            delete_fasta_records(fasta_records);
            return EXIT_FAILURE;
        }
        
        FILE *samfile = stdout;
        char read_name_buffer[FASTQ_BUFFER_SIZE];
        char read_buffer[FASTQ_BUFFER_SIZE];
        char quality_buffer[FASTQ_BUFFER_SIZE];
        while (fastq_parse_next_record(fastq_file, (char*)&read_name_buffer,
                                       (char*)&read_buffer, (char*)&quality_buffer)) {
            
            size_t n = strlen(read_buffer) + (size_t)options.edit_distance;
            char cigar[n + 1], cigar_buffer[n + 1];
            cigar[n] = cigar_buffer[n] = '\0';
            
            size_t no_records = fasta_records->names->used;
            for (size_t seq_no = 0; seq_no < no_records; seq_no++) {
                char *ref_name = fasta_records->names->strings[seq_no];
                struct suffix_array *sa = sa_records->suffix_arrays[seq_no];
                
                search(read_name_buffer, read_buffer, strlen(read_buffer), quality_buffer, ref_name, 0,
                       sa->length - 1, options.edit_distance, cigar, cigar_buffer + n - 1, sa,
                       samfile, &options);
            }
            
        }
        
        delete_fasta_records(fasta_records);
        delete_suffix_array_records(sa_records);
        fclose(fastq_file);
    }
    
    return EXIT_SUCCESS;
}
