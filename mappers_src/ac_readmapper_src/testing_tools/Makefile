CFLAGS=-Wall -DEXTENDED_CIGAR

all: display_trie ac_search edit_cloud

display: display_trie strings.txt
	./display_trie strings.txt && dot < trie.dot -Tpdf > trie.pdf

display_trie: display_trie.o ../trie.o ../queue.o ../strings.o
	cc -o display_trie display_trie.o ../trie.o ../queue.o ../strings.o

search: ac_search bar-database.txt AA-patterns.txt 
	./ac_search bar-database.txt AA-patterns.txt 

ac_search: ac_search.o ../trie.o ../queue.o ../strings.o ../string_vector.o ../string_vector_vector.o ../aho_corasick.o
	cc -o ac_search ac_search.o  ../trie.o ../queue.o ../strings.o ../string_vector.o ../string_vector_vector.o ../aho_corasick.o

edit_cloud: edit_cloud.o ../edit_distance_generator.o ../cigar.o
	cc -o edit_cloud edit_cloud.o ../edit_distance_generator.o ../cigar.o

depend:
	makedepend *.c

clean:
	-rm *.o
	-rm display_trie
	-rm ac_search
	-rm edit_cloud

# DO NOT DELETE

ac_search.o: ../trie.h ../strings.h ../string_vector.h
ac_search.o: ../string_vector_vector.h ../aho_corasick.h ../trie.h
display_trie.o: ../trie.h ../strings.h
edit_cloud.o: ../edit_distance_generator.h
