# Compare different tools for filtering out human reads


## Command-line parameter description

For each human filtering tool `-t ` run 
filtering on all samples  `-s` and
save output `-o`:

```bash
anti_human.py -s samples.dat -t tools.dat -o results.dat
```

## Input files description

* tools file with list of tools for testing
* samples file with list of test cases, each case is a 
pair of FASTQ and annotation file for reads
* fastq files 
* annotation files for each fastq file

### Tools file

For each tool run the comparison, each tool
is on its own line:

```R
bmtagger
```

Here is list of tools supported:

1. bmtagger

### Test cases

Main input file list all of the test cases. Each test
case is a pair of FASTQ file and annotation file. The format is 
tab delimited with three columns name of the test case, FASTQ file
and annotation file:

```R
sample_name1    example1.fastq   read_annotation_for_example1.dat
sample_name2    example2.fastq   read_annotation_for_example2.dat
...
sample_nameM    exampleM.fastq   read_annotation_for_exampleM.dat
```

Each FASTQ file has reads from simulated or real dataset.

## Annotation file

Annotation file annotate each read from a FASTQ file. The format is
tab delimited with columns:

```R
#readId    isHumanRead
read_id1    Y
read_id2    N
...
read_idK    N
```

Header line is optional. All ids should match to ids from the corresponding 
FASTQ file.

## Output file

Creates tab-delimited file with the following field:

```R
#sample\_name   tool    TruePositive TrueNegative FalsePositive FalseNegative ...
```


