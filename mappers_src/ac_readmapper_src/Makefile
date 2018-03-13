source_files = $(wildcard *.c)
object_files = $(source_files:.c=.o)

ac_readmapper: $(object_files)
	cc -o ac_readmapper $(object_files)

depend:
	makedepend $(source_files)
# DO NOT DELETE
