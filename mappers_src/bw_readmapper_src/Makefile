CFLAGS=-Wall -O3 -std=c99

source_files = $(wildcard *.c)
object_files = $(source_files:.c=.o)

bw_readmapper: $(object_files)
	cc -o bw_readmapper $(object_files)

clean:
	-rm bw_readmapper
	-rm *.o

depend:
	makedepend $(source_files)

# DO NOT DELETE

bw_readmap.o: fasta.h string_vector.h size_vector.h fastq.h sam.h search.h
bw_readmap.o: suffix_array_records.h suffix_array.h options.h
cigar.o: cigar.h
fasta.o: fasta.h string_vector.h size_vector.h strings.h
fastq.o: fastq.h strings.h
options.o: options.h
pair_stack.o: pair_stack.h
sam.o: sam.h
search.o: cigar.h sam.h search.h suffix_array_records.h fasta.h
search.o: string_vector.h size_vector.h suffix_array.h options.h strings.h
size_vector.o: size_vector.h
string_vector.o: string_vector.h strings.h
strings.o: strings.h
suffix_array.o: suffix_array.h strings.h pair_stack.h
suffix_array_records.o: suffix_array_records.h fasta.h string_vector.h
suffix_array_records.o: size_vector.h suffix_array.h
