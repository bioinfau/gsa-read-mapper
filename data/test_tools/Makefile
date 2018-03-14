CFLAGS=-Wall

source_files = $(wildcard *.c)
SHARED_OBJECTS := fasta.o size_vector.o string_vector.o strings.o

display_match: display_match.o $(SHARED_OBJECTS)
	cc -o display_match display_match.o $(SHARED_OBJECTS)

clean:
	-rm display_match
	-rm *.o

depend: $(source_files)
	makedepend $(source_files)

# DO NOT DELETE

fasta.o: fasta.h string_vector.h size_vector.h strings.h
size_vector.o: size_vector.h
string_vector.o: string_vector.h strings.h
strings.o: strings.h
