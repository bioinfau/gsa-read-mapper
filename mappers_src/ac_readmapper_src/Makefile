CFLAGS=-Wall -O3 -std=c99

source_files = $(wildcard *.c)
object_files = $(source_files:.c=.o)

ac_readmapper: $(object_files)
	cc -o ac_readmapper $(object_files)

clean:
	-rm ac_readmapper
	-rm *.o

depend:
	makedepend $(source_files)
# DO NOT DELETE

ac_readmap.o: fasta.h string_vector.h size_vector.h fastq.h sam.h
ac_readmap.o: string_vector_vector.h aho_corasick.h trie.h
ac_readmap.o: edit_distance_generator.h options.h
aho_corasick.o: aho_corasick.h trie.h
cigar.o: cigar.h
edit_distance_generator.o: edit_distance_generator.h options.h cigar.h
fasta.o: fasta.h string_vector.h size_vector.h strings.h
fastq.o: fastq.h strings.h
match.o: match.h
options.o: options.h
pair_stack.o: pair_stack.h
queue.o: queue.h
sam.o: sam.h
size_vector.o: size_vector.h
string_vector.o: string_vector.h strings.h
string_vector_vector.o: string_vector_vector.h string_vector.h
strings.o: strings.h
trie.o: trie.h queue.h
