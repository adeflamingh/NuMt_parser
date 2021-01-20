# Numt Parser

`Numt Parser` a software for the detection and filtering of Nuclear Mitochondrial pseudogene (numt) contamination in Mitochondrial shotgun sequencing datasets.

Copyright 2020: Angel G. Rivera-Colon & Alida de Flamingh

## Installation and Usage Requirements

* Python 3

## Usage

### Pre-*Numt Parser*

Before processing, raw reads must be processed (trimmed to desired length and filtering) as desired.

Reads are then aligned using a short read aligner (such as `bwa` or `bowtie2`) to both a true mt sequence and a characterized numt reference.

Resulting alignments can then be additionally filtered (using tools such as `samtools`, `Picard`, etc.) to removed unmapped reads and low quality alignments. Final filtered alignments should be saved as `SAM` files for compatibility with `Numt Parser`.

### *Numt Parser* Analysis

`Numt Parser` requires four (4) inputs:

1. Cytoplasmic mitochondria reference file (FASTA)
2. numt reference file (FASTA)
3. Cytoplasmic mitochondria alignment file (SAM)
4. numt alignment file (SAM)

A fifth parameter naming the location and name of the resulting output table (TSV) is also required.

#### Usage example

```sh
$ numt_parser.py \
    --mt-fasta /path/to/cymt.fasta \
    --numt-fasta /path/to/numt.fasta \
    --mt-sam /path/to/mt_alignment.sam \
    --numt-sam /path/to/numt_alignment.sam \
    --outfile /path/to/numt_parser_output.tsv
```

### *Numt Parser* output

The output of `Numt Parser` is a table containing the identity statistics of each processed read.

```sh
#read_ID  mt_aln_bp  mt_mismatch  mt_identity  numt_aln_bp  numt_mismatch  numt_identity  candidate
Read_01   100        0            1.000000     100          0              1.000000       undetermined
Read_02   48         0            1.000000     48           0              1.000000       undetermined
Read_03   81         5            0.938272     None         None           None           cymt
Read_04   73         5            0.931507     76           1              0.986842       numt
Read_05   67         3            0.955224     67           0              1.000000       numt
Read_06   None       None         None         65           0              1.000000       numt
Read_07   100        0            1.000000     100          6              0.940000       cymt
Read_08   None       None         None         66           0              1.000000       numt
Read_09   80         3            0.962500     80           4              0.950000       cymt
```

### Post-*Numt Parser* processing

Using the resulting output table, raw reads files can be filtered to obtain specific reads originating from either cytoplasmic mitochondria or numt pseudogene templates:

1. Create a "Read ID list" text file with the names of the reads to include or exclude from the dataset. For example:

```sh
$ cat numt_parser_output.tsv | grep -v '^#' | grep 'cymt' | cut -f1 > Readlist.txt
```

The command above will filter reads tagged as `numt` or `undetermined`, retaining those of cytoplasmic mitochondrial origin (`cymt`).

```sh
$ cat ReadList.txt
read_03
read_07
read_09
```

2. Use the "Read ID list" to filter the original BAM alignment file for reads mapping to the mt Alignment file. E.g. in `Picard` (<https://github.com/broadinstitute/picard>) use the `FilterSamReads` function:

```sh
$ java -jar picard.jar FilterSamReads \
    I=Mt_alignment.bam \
    O=Mt_alignment_filtered.bam \
    READ_LIST_FILE=Readlist.txt \
    FILTER=includeReadList
 ```

3. The resulting output BAM file will contain only reads that have been identified by `Numt Parser` as being of putitive mt origin. Similar read filtering can be done in the program `seqtk` (<https://github.com/lh3/seqtk>) or `samtools` (<http://www.htslib.org/doc/samtools.html>).

## Authors

Alida de Flamingh & Angel G. Rivera-Colon

University of Illinos at Urbana-Champaign
