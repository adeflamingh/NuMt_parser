# numt parser

`numt parser` a software for the detection and filtering of Nuclear Mitochondrial pseudogene (numt) contamination in Mitochondrial shotgun sequencing datasets.

Copyright 2020: Angel G. Rivera-Colon & Alida de Flamingh

## Installation and Usage Requirements

* Python 3

## Usage

### Pre-*numt parser*

Before processing, raw reads must be processed (trimmed to desired length and filtering) as desired.

Reads are then aligned using a short read aligner (such as `bwa` or `bowtie2`) to both a true mt sequence and a characterized numt reference.

Resulting alignments can then be additionally filtered (using tools such as `samtools`, `Picard`, etc.) to removed unmapped reads and low quality alignments. Final filtered alignments should be saved as `SAM` files for compatibility with `numt parser`.

### *numt parser* Analysis

`numt parser` requires four (4) inputs:

1. mt Reference file (FASTA)
2. numt Reference file (FASTA)
3. mt Alignment file (SAM)
4. numt Alignment file (SAM)

A fifth parameter naming the location and name of the resulting output table (TSV) is also required.

#### Usage example

```sh
numt_parser.py \
    --mt-fasta /path/to/mitochondria.fasta \
    --numt-fasta /path/to/numt.fasta \
    --mt-sam /path/to/mt_alignment.sam \
    --numt-sam /path/to/numt_alignment.sam \
    --outfile /path/to/numt_parser_output.tsv
```

### *numt parser* output

The output of `numt parser` is a table containing the identity statistics of each processed read.

```sh
#read_ID  mt_aln_bp  mt_mismatch  mt_identity  numt_aln_bp  numt_mismatch  numt_identity  candidate
read_01   100        0            1.000000     100          0              1.000000       Unknown
read_02   48         0            1.000000     48           0              1.000000       Unknown
read_03   81         5            0.938272     None         None           None           mt
read_04   73         5            0.931507     76           1              0.986842       numt
read_05   67         3            0.955224     67           0              1.000000       numt
read_06   None       None         None         65           0              1.000000       numt
read_07   100        0            1.000000     100          6              0.940000       mt
read_08   None       None         None         66           0              1.000000       numt
read_09   80         3            0.962500     80           4              0.950000       mt
```

### Post-*numt parser* processing

Using the resulting output table, raw reads files can be filtered to obtain specific reads originating from either mt or numt templates:

1. Create a "Read ID list" text file with the names of the reads to include or exclude from the dataset. For example:

```sh
cat numt_parser_output.tsv | grep -v '^#' | grep -E -v 'numt|Und' | cut -f1 > Readlist.txt
```

The command above will filter reads tagged as `numt` or `Unknown`, retaining those of mitochondrial origin.

```sh
read_03
read_07
read_09
```

2. Use the "Read ID list" to filter the original BAM alignment file for reads mapping to the mt Alignment file. E.g. in `Picard` (<https://github.com/broadinstitute/picard>) use the `FilterSamReads` function:

```sh
java -jar picard.jar FilterSamReads \
    I=Mt_alignment.bam \
    O=Mt_alignment_filtered.bam \
    READ_LIST_FILE=Readlist.txt \
    FILTER=includeReadList
 ```

3. The resulting output BAM file will contain only reads that have been identified by `numt parser` as being of putitive mt origin. Similar read filtering can be done in the program `seqtk` (<https://github.com/lh3/seqtk>) or `samtools` (<http://www.htslib.org/doc/samtools.html>).

## Authors

Alida de Flamigh & Angel G. Rivera-Colon

University of Illinos at Urbana-Champaign
