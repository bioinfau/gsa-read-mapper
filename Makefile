

all: mappers test evaluate

mappers:
	cd mappers_src && make

test: mappers
	(export PATH=mappers_src:${PATH} ; cd evaluation && ./test_mappers.sh)

evaluate: mappers
	(export PATH=mappers_src:${PATH} ; cd evaluation && ./evaluate_mappers.sh)


