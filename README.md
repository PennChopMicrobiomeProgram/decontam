# Run different tools for filtering out human reads

## Command-line parameter description

For each human filtering tool `-t ` run 
filtering on all samples  `-s` and
save output `-o`:

```bash
anti _ human.py -s samples.dat -t tools.dat -o results.dat
```

## Input files description

* *tools* file with list of tools for testing
* *samples* file with list of samples, each sample is a paired-end run with 2 FASTQ files.


### Tools file

For each tool run the comparison, each tool
is on its own line:

```R
bmtagger
```

Here is list of tools supported:

1. bmtagger
2. all _ human 
3. none _ human
4. bwa
5. bowtie2
6. blat

### Samples file

Main input file list all of the samples. Each sample
is a pair of name and FASTQ file. The format is 
tab delimited with 3  columns(name, forward FASTQ file, reverse FASTQ file):

```R
sample _ name1    example1 _ R1.fastq    example1 _ R2.fastq
sample _ name2    example2 _ R1.fastq    example2 _ R2.fastq
...                                                    
sample _ nameM    exampleM _ R1.fastq    exampleM _ R2.fastq
```

Each pair of FASTQ files have reads from simulated or real dataset.
The pair of forward and reverse fastq files for a sample should have the same read _ ids.

## Output file

Creates tab-delimited file with the following field:

```R
#tool _ name    sample _ name   read _ id is _ human _  prediction
```


## Development
 
To add new human filtering tool, add entry to `toolname _ to _ runner` dictionary
and implement class that has `extract _ human _ reads` method for a FASTQ sample file.



## Testing

To test program run `nosetests` from directory where README.md is located.
