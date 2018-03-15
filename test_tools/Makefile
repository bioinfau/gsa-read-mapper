CFLAGS=-Wall
LDFLAGS=-L.

SHARED_OBJ := \
	trie.o queue.o strings.o string_vector.o \
	string_vector_vector.o aho_corasick.o \
	fasta.o size_vector.o strings.o \
	queue.o edit_distance_generator.o cigar.o

all: display_match display_trie edit_cloud ac_search

display_match: display_match.o libgsa.a
	cc -o display_match display_match.o -L. -lgsa

display_trie: display_trie.o libgsa.a
		cc -o display_trie display_trie.o -L. -lgsa

edit_cloud: edit_cloud.o libgsa.a
		cc -o edit_cloud edit_cloud.o -L. -lgsa

ac_search: ac_search.o libgsa.a
		cc -o ac_search ac_search.o -L. -lgsa

libgsa.a: $(SHARED_OBJ)
		ar -csru libgsa.a $(SHARED_OBJ)

clean:
	-rm display_match
	-rm display_trie
	-rm edit_cloud
	-rm ac_search
	-rm *.o

depend:
	makedepend $(source_files)

# DO NOT DELETE

ac_search.o: trie.h strings.h string_vector.h string_vector_vector.h
ac_search.o: aho_corasick.h
aho_corasick.o: aho_corasick.h trie.h
cigar.o: cigar.h
display_match.o: fasta.h string_vector.h size_vector.h
display_trie.o: trie.h strings.h
edit_cloud.o: edit_distance_generator.h
edit_distance_generator.o: edit_distance_generator.h cigar.h
fasta.o: fasta.h string_vector.h size_vector.h strings.h
queue.o: queue.h
size_vector.o: size_vector.h
string_vector.o: string_vector.h strings.h
string_vector_vector.o: string_vector_vector.h string_vector.h
strings.o: strings.h
trie.o: trie.h queue.h
