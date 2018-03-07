
mappers=$(wildcard *_src)

all:
	echo $(mappers)
	for m in $(mappers); do       		 \
		mapper_name=`basename $$m _src`	;\
		echo $$mapper_name		 		;\
		(cd $$m && make)  		 		;\
		cp $$m/$$mapper_name $$mapper_name ;\
	done

depend:
	echo $(mappers)
	for m in $(mappers); do       \
		echo $$m 				 ;\
		(cd $$m && make depend)  ;\
	done

