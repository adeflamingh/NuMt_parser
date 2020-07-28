# NuMt Parser

`NuMt Parser` a software for the detection and filtering of Nuclear Mitochondrial pseudogene (NuMt) contamination in Mitochondrial shotgun sequencing datasets.

Copyright 2020: Angel G. Rivera-Colon & Alida de Flamingh

## Installation and Usage Requirements

* Python 3

## Usage

### Pre-NuMt Parser

Before processing, raw reads must be processed (trimmed to desired length and filtering) as desired.

Reads are then aligned using a short read aligner (such as `bwa` or `bowtie2`) to both a true Mt sequence and a characterized NuMt reference.

Resulting alignments can then be additionally filtered (using tools such as `samtools`, `Picard`, etc.) to removed unmapped reads and low quality alignments. Final filteresd alignments should be saved as `SAM` files for compatibility with `NuMt Parser`.

### NuMt Parser Analysis

`NuMt Parser` requires four (4) inputs:

1. Mt Reference file (FASTA)
2. NuMt Reference file (FASTA)
3. Mt Alignment file (SAM)
4. NuMt Alignment file (SAM)

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

### NuMt Parser output

The output of `NuMt Parser` is a table containing the identity statistics of each processed read.

```sh
#Read_ID  Mt_aln_bp  Mt_mismatch  Mt_identity  NuMt_aln_bp  NuMt_mismatch  NuMt_identity  Candidate
read_01   100        0            1.000000     100          0              1.000000       Unknown
read_02   48         0            1.000000     48           0              1.000000       Unknown
read_03   81         5            0.938272     None         None           None           Mt
read_04   73         5            0.931507     76           1              0.986842       NuMt
read_05   67         3            0.955224     67           0              1.000000       NuMt
read_06   None       None         None         65           0              1.000000       NuMt
read_07   100        0            1.000000     100          6              0.940000       Mt
read_08   None       None         None         66           0              1.000000       NuMt
read_00   80         3            0.962500     80           4              0.950000       Mt
```

### Post-NuMt Parser processing

Using the resulting output table, raw reads files can be filtered to obtain specific reads originating from either Mt or NuMt templates:

1)Create a "Read ID list" text file with the names of the reads to include or exclude from the dataset 
    E.g. Readlist.txt for the table above with reads originating from Mt:
    read_03
    read_07
    read_00
 
 2)Use the "Read ID list" to filter the original bam alignment file for reads mapping to the Mt Alignment file
    E.g. in Picard (https://github.com/broadinstitute/picard) use the "FilterSamReads" function
   
   ```sh
    picard FilterSamReads \
    I=Mt_alignment.bam \
    O=Mt_alignment_filtered.bam \
    READ_LIST_FILE=Readlist.txt \
    FILTER=includeReadList
    ```
 
 3)The resulting output bam file will contain only reads that have been identified by NuMt parser as being of putitive Mt origin.  
 Similar read filtering can be done in the program seqtk (https://github.com/lh3/seqtk)
 

## Authors

Alida de Flamigh & Angel G. Rivera-Colon

University of Illinos at Urbana-Champaign
