# Run different tools for filtering out human reads


## Command-line parameter description

For each human filtering tool `-t ` run 
filtering on all samples  `-s` and
save output `-o`:

```bash
anti_human.py -s samples.dat -t tools.dat -o results.dat
```

## Input files description

* tools file with list of tools for testing
* samples file with list of samples, each sample is a FASTQ file


### Tools file

For each tool run the comparison, each tool
is on its own line:

```R
bmtagger
```

Here is list of tools supported:

1. bmtagger
2. all_human 
3. none_human
4. bwa
5. bowtie2
6. blat

### Samples file

Main input file list all of the samples. Each sample
is a pair of name and FASTQ file. The format is 
tab delimited with ??? columns name, FASTQ file :

```R
sample_name1    example1.fastq
sample_name2    example2.fastq
...
sample_nameM    exampleM.fastq
```

Each FASTQ file has reads from simulated or real dataset.

## Output file

Creates tab-delimited file with the following field:

```R
#tool\_name    sample\_name   read\_id is\_human\_prediction
```


## Development
 
To add new human filtering tool, add entry to `toolname\_to\_runner` dictionary
and implement class that has `extract\_human\_reads` method for a FASTQ sample file.



## Testing

To test program run `nosetests` from directory where README.md is located.
