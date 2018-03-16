
/*
 Readmapper based on explicitly generating approximate matches and
 using the Aho-Corasick algorithm for matching.
*/

#include "fasta.h"
#include "fastq.h"
#include "sam.h"
#include "string_vector_vector.h"
#include "aho_corasick.h"
#include "edit_distance_generator.h"
#include "options.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>

struct search_info {
    struct fasta_records *records;
    struct options *options;
    FILE *sam_file;
};

static struct search_info *empty_search_info(struct options *options)
{
    struct search_info *info =
        (struct search_info*)malloc(sizeof(struct search_info));
    info->options = options;
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
    FILE *sam_file;
    
    struct string_vector *patterns;
    struct string_vector_vector *cigars;
    struct trie *patterns_trie;
};

static struct read_search_info *empty_read_search_info()
{
    struct read_search_info *info =
        (struct read_search_info*)malloc(sizeof(struct read_search_info));
    
    info->ref_name = 0;
    info->read_name = 0;
    info->read = 0;
    info->quality = 0;
    
    info->patterns = empty_string_vector(256); // arbitrary start size...
    info->cigars = empty_string_vector_vector(256); // arbitrary start size...
    info->patterns_trie = empty_trie();
    
    return info;
}

static void delete_read_search_info(struct read_search_info *info)
{
    delete_string_vector(info->patterns);
    delete_string_vector_vector(info->cigars);
    delete_trie(info->patterns_trie);
    free(info);
}

static void build_trie_callback(const char *pattern, const char *cigar, void * data)
{
    struct read_search_info *info = (struct read_search_info*)data;
    
    // patterns generated when we explore the neighbourhood of a read are not unique
    // so we need to check if we have seen it before
    if (string_in_trie(info->patterns_trie, pattern)) {
        // The pattern is already in the tree, but if we are called here
        // we have a new CIGAR for the same pattern.
        struct trie *node = get_trie_node(info->patterns_trie, pattern);
        add_string_copy_to_vector(info->cigars, node->string_label, cigar);

    } else {
        // NB: the order is important here -- info->patterns->used will be updated
        // when we add the pattern to the vector, so we insert in the trie first.
        int index = append_vector(info->cigars);
        add_string_to_trie(info->patterns_trie, pattern, index);
        add_string_copy(info->patterns, pattern);
        add_string_copy_to_vector(info->cigars, index, cigar);
    }
}

static void match_callback(int string_label, size_t index, void * data)
{
    struct read_search_info *info = (struct read_search_info*)data;
    const char *str = info->patterns->strings[string_label];
    struct string_vector *cigars = info->cigars->string_vectors[string_label];
    size_t n = strlen(str);
    size_t start_index = index - n + 1 + 1; // +1 for start correction and +1 for 1-indexed
    for (int i = 0; i < cigars->used; i++) {
        sam_line(info->sam_file,
                 info->read_name, info->ref_name, start_index,
                 cigars->strings[i],
                 info->read,
                 info->quality);
    }
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
    info->sam_file = search_info->sam_file;
    info->read = read;
    info->quality = quality;
    
    generate_all_neighbours(read, "ACGT", search_info->options->edit_distance,
                            build_trie_callback, info, search_info->options);
    compute_failure_links(info->patterns_trie);
    
    info->read_name = read_name;
    for (int i = 0; i < search_info->records->names->used; ++i) {
        info->ref_name = search_info->records->names->strings[i];
        const char *ref = search_info->records->sequences->strings[i];
        size_t n = search_info->records->seq_sizes->sizes[i];
        aho_corasick_match(ref, n, info->patterns_trie, match_callback, info);
    }
    
    delete_read_search_info(info);
}

int main(int argc, char * argv[])
{
    const char *prog_name = argv[0];
    
    struct options options;
    options.edit_distance = 0;
    options.extended_cigars = false;
    options.verbose = false;
    
    int opt;
    static struct option longopts[] = {
        { "help",       no_argument,            NULL,           'h' },
        { "distance",   required_argument,      NULL,           'd' },
        { "extended-cigar",   no_argument,      NULL,           'x' },
        { NULL,         0,                      NULL,            0  }
    };
    while ((opt = getopt_long(argc, argv, "hd:x", longopts, NULL)) != -1) {
        switch (opt) {
            case 'h':
                printf("Usage: %s [options] ref.fa reads.fq\n\n", prog_name);
                printf("Options:\n");
                printf("\t-h | --help:\t\t Show this message.\n");
                printf("\t-x | --extended-cigar:\t Use extended CIGAR format in SAM output.\n");
                printf("\t-d | --distance:\t Maximum edit distance for the search.\n");
                printf("\n\n");
                return EXIT_SUCCESS;
                
            case 'd':
                options.edit_distance = atoi(optarg);
                break;

            case 'x':
                options.extended_cigars = true;
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
    
    struct search_info *search_info = empty_search_info(&options);
    read_fasta_records(search_info->records, fasta_file);
    fclose(fasta_file);
    
    search_info->sam_file = stdout;
    
    scan_fastq(fastq_file, read_callback, search_info);
    delete_search_info(search_info);
    fclose(fastq_file);
    
    return EXIT_SUCCESS;
}
