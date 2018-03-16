
#ifndef FASTQ_H
#define FASTQ_H

#include <stdio.h>
#include <stdbool.h>

typedef void (*fastq_read_callback_func)(const char *read_name,
                                         const char *read,
                                         const char *quality,
                                         void * callback_data);

void scan_fastq(FILE *file, fastq_read_callback_func callback, void * callback_data);
bool fastq_parse_next_record(FILE *file, char *read_name_buffer,
                             char *read_buffer, char *quality_buffer);

#endif
