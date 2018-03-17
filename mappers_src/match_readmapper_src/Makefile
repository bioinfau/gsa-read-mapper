CFLAGS=-Wall -O3 -std=c99

source_files = $(wildcard *.c)
object_files = $(source_files:.c=.o)

match_readmapper: $(object_files)
	cc -o match_readmapper $(object_files)

clean:
	-rm match_readmapper
	-rm *.o

depend:
	makedepend $(source_files)

# DO NOT DELETE

cigar.o: cigar.h
edit_distance_generator.o: edit_distance_generator.h options.h cigar.h
fasta.o: fasta.h string_vector.h size_vector.h strings.h
fastq.o: fastq.h strings.h
match.o: match.h
match_readmap.o: match.h suffix_array.h fasta.h string_vector.h size_vector.h
match_readmap.o: fastq.h sam.h edit_distance_generator.h options.h
options.o: options.h
pair_stack.o: pair_stack.h
queue.o: queue.h
sam.o: sam.h
size_vector.o: size_vector.h
string_vector.o: string_vector.h strings.h
strings.o: strings.h
suffix_array.o: suffix_array.h match.h strings.h pair_stack.h
trie.o: trie.h queue.h
