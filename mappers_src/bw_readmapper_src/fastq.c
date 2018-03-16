
#include "fastq.h"
#include "strings.h"

#include <stdlib.h>
#include <string.h>

#define MAX_LINE_SIZE 1024

void scan_fastq(FILE *file, fastq_read_callback_func callback, void * callback_data)
{
    char buffer[MAX_LINE_SIZE];
    
    while (fgets(buffer, MAX_LINE_SIZE, file) != 0) {
        char *name = string_copy(strtok(buffer+1, "\n"));
        fgets(buffer, MAX_LINE_SIZE, file);
        char *seq = string_copy(strtok(buffer, "\n"));
        fgets(buffer, MAX_LINE_SIZE, file);
        fgets(buffer, MAX_LINE_SIZE, file);
        char *qual = string_copy(strtok(buffer, "\n"));
        
        callback(name, seq, qual, callback_data);
        
        free(name);
        free(seq);
        free(qual);
    }
}

bool fastq_parse_next_record(FILE *file, char *read_name_buffer,
                             char *read_buffer, char *quality_buffer)
{
    char buffer[MAX_LINE_SIZE];
    
    if (fgets(buffer, MAX_LINE_SIZE, file) == 0) return false; // read name line
    strcpy(read_name_buffer, strtok(buffer+1, "\n"));
    if (fgets(buffer, MAX_LINE_SIZE, file) == 0) return false; // read line
    strcpy(read_buffer, strtok(buffer, "\n"));
    if (fgets(buffer, MAX_LINE_SIZE, file) == 0) return false; // '+' line
    if (fgets(buffer, MAX_LINE_SIZE, file) == 0) return false; // quality line
    strcpy(quality_buffer, strtok(buffer, "\n"));
    
    return true;
}
