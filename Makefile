

all: mappers test test evaluate

mappers:
	cd mappers_src && make

test: test_exact test_approximative

test_exact: mappers
	(export PATH=mappers_src:${PATH} ; cd evaluation && ./test_mappers_exact.sh)

test_approximative: mappers
	(export PATH=mappers_src:${PATH} ; cd evaluation && ./test_mappers_approximative.sh)

evaluate: mappers
	(export PATH=mappers_src:${PATH} ; cd evaluation && ./evaluate_mappers.sh)


