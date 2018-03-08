
mappers=$(wildcard *_src)

all:
	echo $(mappers)
	for m in $(mappers); do       		                  \
		mapper_name=`basename $$m _src`	                 ;\
		echo $$mapper_name		 		                 ;\
		(cd $$m && make)  		 		                 ;\
		if [ $$m/$$mapper_name -nt $$mapper_name ]; then  \
			cp $$m/$$mapper_name $$mapper_name           ;\
		fi                                               ;\
	done

depend:
	echo $(mappers)
	for m in $(mappers); do       \
		echo $$m 				 ;\
		(cd $$m && make depend)  ;\
	done

