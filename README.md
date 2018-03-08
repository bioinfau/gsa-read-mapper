# Read-mapper project for Genome Scale Algorithms

## Downloading the testing and evaluation scripts

You can download the scripts as [a ZIP file](https://github.com/bioinfau/gsa-read-mapper/archive/master.zip), but if you are familiar with Git, then an easier solution is to clone this repository

```sh
git clone https://github.com/bioinfau/gsa-read-mapper.git
```

If you plan to use the clone as a shared repository to work with during the project, you might want to fork it before you clone your own copy.

## Building mappers

Once you have downloaded or cloned the code you can build the read-mappers I have included here using the command

```sh
make mappers
```

If all goes well, this should build three mappers in the directory `mappers_src`:
 * `bwa` — the Burrows-Wheeler based read-mapper you experimented with in the first week’s exercises.
 * `match_readmapper` — A read-mapper based on constructing the edit-distance cloud around reads and searching for them with exact pattern matching algorithms. You can pick the algorithm to use via the `--algorithm` option.
 * `ac_readmapper` — Another read-mapper that builds the edit-distance cloud and search for matches this way, but using the Aho-Corsick algorithm for the search.

### Adding your own mapper

To add your own read-mapper to the build setup, you can simply put it in a sub-directory of `mappers_src`. The name of your sub-directory most be the name of your executable followed by `_src` and inside the sub-directory you must have a Makefile for building your executable. When you invoke

```sh
make mappers
```

the build setup will enter all directories that ends in `_src` and call `make` inside them. It expects that this builds an executable that is named the same as the directory, except for the suffix `_src`. The executable is then copied into the `mappers` directory.

### Testing mappers

If you invoke 

```sh
make test
```

you will run the script [`evaluation/test_mapper.sh`](https://github.com/mailund/gsa-read-mapper/blob/master/evaluation/test_mappers.sh). This script uses a reference implementation to build a SAM file and then it tests that all the other mappers specified in the script produce the same SAM file.

You can modify the [header of the script](https://github.com/mailund/gsa-read-mapper/blob/a748068714fabeb8989382664c1dfea8e87fb79b/evaluation/test_mappers.sh#L3-L25) to configure how the tests are run.

The relevant variables you can modify are:
 * `ref_mapper` — this is the read-mapper the others will be compared against. If you do not change this variable, the reference read-mapper is `match_readmapper`. The performance of this mapper, however, is such that you cannot use for approximate pattern matching. It is simply too slow. So replace the reference when you test that your own mapper can find approximative matches. You can, for example, use the `ac_mapper` instead.
 * `mappers` — this is a list of the mappers to test against the reference. This is where you want to add your own read-mapper.
 * `report_file` and `log_file` determines where results and logging is written. There is no need to change these.
 * `d` — this is the maximum edit-distance to search in an approximate pattern matching. Unless you change it, it is set to zero, which means exact pattern matching. If you increase the distance, you probably want to use a different `ref_mapper`.
 * `reference` — this is the file that contains the reference genome. Unless you change it, it is a short prefix of the gorilla chromosome 1 where I have replaced ’N’ characters with random ‘A’, ‘C’, ‘G’, or ’T’, characters.
 * `reads` — this is the file containing the reads. By default it is a file that contains 10 reads of length 10 that I have copied from the reference string and modified up to distance d=2.

```sh
## Modify here to add or remove mappers or change options
## =============================================================

# The mapper we use as the goal to hit
ref_mapper=match_readmapper

# list of read-mappers to evaluate
mappers="ac_readmapper"

# file name for report
report_file=../test-report.txt
log_file=../test.log

# max edit distance to explore
d=0

# Reference genome
reference=../data/gorGor3-small-noN.fa

# Reads
reads=../data/sim-reads-d2-tiny.fq

## =============================================================
```

## Evaluating mapper performance



## The data files

## Evaluation

The current status of the read-mappers are:

![](evaluation-report.txt.png)

