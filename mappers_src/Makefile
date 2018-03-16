
mappers=$(wildcard *_src)

all:
	echo $(mappers)
	for m in $(mappers); do       		                  \
		mapper_name=`basename $$m _src`	                 ;\
		echo $$mapper_name		 		                       ;\
		(cd $$m && make)  		 		                       ;\
		if [ ! -e $$mapper_name ]; then                   \
			cp $$m/$$mapper_name $$mapper_name             ;\
		fi                                               ;\
		if [ $$m/$$mapper_name -nt $$mapper_name ]; then  \
			cp $$m/$$mapper_name $$mapper_name             ;\
		fi                                               ;\
	done                                               ;\
	for m in $(mappers); do       		                  \
		mapper_name=`basename $$m _src`	                 ;\
		if [ ! -x $$mapper_name ]; then                   \
		  echo "Hey! where did $${mapper_name} go?"      ;\
			exit 1                                         ;\
		fi                                               ;\
	done


depend:
	echo $(mappers)
	for m in $(mappers); do       \
		echo $$m 				 ;\
		(cd $$m && make depend)  ;\
	done
