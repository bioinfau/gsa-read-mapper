
all: display_trie ac_search edit_cloud

display: display_trie strings.txt
	./display_trie strings.txt && dot < trie.dot -Tpdf > trie.pdf

display_trie: display_trie.o ../trie.o ../queue.o ../strings.o
	cc -o display_trie display_trie.o ../trie.o ../queue.o ../strings.o

search: ac_search strings.txt database.txt
	./ac_search database.txt strings.txt

ac_search: ac_search.o ../trie.o ../queue.o ../strings.o ../string_vector.o ../aho_corasick.o
	cc -o ac_search ac_search.o  ../trie.o ../queue.o ../strings.o ../string_vector.o ../aho_corasick.o

edit_cloud: edit_cloud.o ../edit_distance_generator.o ../cigar.o
	cc -o edit_cloud edit_cloud.o ../edit_distance_generator.o ../cigar.o
