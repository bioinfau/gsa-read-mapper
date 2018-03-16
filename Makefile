

all: mappers test test evaluate

mappers:
	cd mappers_src && make

test: tests
tests: test_exact test_approximative

test_exact: mappers
	(export PATH=${PWD}/mappers_src:${PATH} ; cd evaluation && ./test_mappers_exact.sh)

test_approximative: mappers
	(export PATH=${PWD}/mappers_src:${PATH} ; cd evaluation && ./test_mappers_approximative.sh)

evaluate: evaluate_exact evaluate_approximative

evaluate_exact: mappers
	(export PATH=${PWD}/mappers_src:${PATH} ; cd evaluation && ./evaluate_mappers_exact.sh)

evaluate_approximative: mappers
	(export PATH=${PWD}/mappers_src:${PATH} ; cd evaluation && ./evaluate_mappers_approximative.sh)
