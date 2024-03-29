# Numt Parser

`Numt Parser` a software for the detection and filtering of Nuclear Mitochondrial pseudogene (numt) contamination in Mitochondrial shotgun sequencing datasets.

Copyright 2020: Angel G. Rivera-Colon & Alida de Flamingh

## Installation and Usage Requirements

* Python 3

## Usage

### Pre-*Numt Parser*

Before filtering, raw reads must be processed (trimmed to desired length and filtering) as desired.

Reads are then aligned using a short read aligner (such as `bwa` or `bowtie2`) to both a true mt sequence and a characterized numt reference.

Resulting alignments can then be additionally filtered (using tools such as `samtools`, `Picard`, etc.) to remove unmapped reads and low quality alignments. Final filtered alignments should be saved as `SAM` files for compatibility with `Numt Parser`.

### *Numt Parser* Analysis

`Numt Parser` requires four (4) inputs files:

1. Cytoplasmic mitochondria reference file (FASTA)
2. numt reference file (FASTA)
3. Cytoplasmic mitochondria alignment file (SAM)
4. numt alignment file (SAM)

A fifth parameter specifying the location and name of the resulting output table (TSV) is also required.

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
#read_ID  read_in_pair  mt_aln_bp  mt_mismatch  mt_identity  numt_aln_bp  numt_mismatch  numt_identity  candidate
Read_01   0             100        0            1.000000     100          0              1.000000       undetermined
Read_02   0             48         0            1.000000     48           0              1.000000       undetermined
Read_03   0             81         5            0.938272     None         None           None           cymt
Read_04   0             73         5            0.931507     76           1              0.986842       numt
Read_05   0             67         3            0.955224     67           0              1.000000       numt
Read_06   0             None       None         None         65           0              1.000000       numt
Read_07   0             100        0            1.000000     100          6              0.940000       cymt
Read_08   0             None       None         None         66           0              1.000000       numt
Read_09   0             80         3            0.962500     80           4              0.950000       cymt
```

When the input alignments originate from single-end reads, the `read_in_pair` column will show a value of `0`, since none of the reads are paired. When the alignments originate form paired-end reads, the output table will instead show a value of `1` or `2` in the `read_in_pair` column to denote the first or second read in the pair, respectively. In this case, the `read_id` will be seen up to two times, once per read in the pair. In cases for which singleton paired reads are present (e.g., a read `2` without its paired read `1` present in the file), `read_in_pair` information will be reported as described on the [SAM flag](https://broadinstitute.github.io/picard/explain-flags.html) for the alignment.

```sh
#read_ID  read_in_pair  mt_aln_bp  mt_mismatch  mt_identity  numt_aln_bp  numt_mismatch  numt_identity  candidate
Read_01   1             85         4            0.952941     85           0              1.000000       numt
Read_01   2             63         0            1.000000     63           0              1.000000       undetermined
Read_02   1             100        0            1.000000     None         None           None           cymt
Read_02   2             73         0            1.000000     None         None           None           cymt
Read_03   1             67         3            0.955224     67           0              1.000000       numt
Read_03   2             52         4            0.923077     52           0              1.000000       numt
Read_04   1             100        0            1.000000     100          6              0.940000       cymt
Read_04   2             85         4            0.952941     85           0              1.000000       numt
```

When analyzing paired-end data, each read in the pair is processed independently. Therefore, it is possible for each read in the pair to be tagged as having a different origin. In the example above, the first read in the `Read_01` pair is tagged as originating from a `numt` template, while the second read gets tagged as `undertermined`. A more drastic example is observed for the `Read_04` read pair, as the first and second reads in the pair are tagged as originating from different `cymt` and `numt` templates.

### Post-*Numt Parser* processing

Using the resulting output table, raw reads (or alignment) files can be filtered to obtain specific reads originating from either cytoplasmic mitochondria or numt pseudogene templates. An example of such filtering to obtain reads of `cymt` origin can be seen in the steps below:

1. Create a "Read ID list" text file with the names of the reads to include or exclude from the dataset. This process may be slightly different when processing single- or paired-end data.

* **For single-end (or unpaired) reads**

When processing single-end data, each `read_id` is only seen once. Therefore the output table can simply be filtered by rows.

```sh
$ cat numt_parser_output.tsv | grep -v '^#' | grep 'cymt' | cut -f1 > Readlist.txt
```

The command above will filter out reads tagged as `numt` or `undetermined`, retaining those of cytoplasmic mitochondrial origin (`cymt`).

```sh
$ cat ReadList.txt
read_03
read_07
read_09
```

* **For paired-end reads**

When processing alignments from paired-end data, each `read_id` might be seen up to two times in the output table (once per read in the pair). Since each read in the pair is processed independently, it is possible for them to have separate tags. A simple filtering alternative could be to select rows containing the `cymt` tag, and then only selecting `read_id`s seen twice:

```sh
$ cat numt_parser_output.tsv | grep -v '^#' | grep 'cymt' | cut -f1 | sort | \
    uniq -c | awk '$1 == 2 {print $2}' > Readlist.txt
```

In this example, after selecting only the rows contining the `cymt` tag, we tally the `read_id`s and select those seen two times. This guarantees that we are only keeping those read for which both the first and second reads in the paired are putatively of `cymt` origin. Note that more complex filtering criteria are also possible.

2. Use the "Read ID list" to filter the original BAM alignment file to retain reads of putative `cymt` origin, e.g., in [`Picard`](<https://github.com/broadinstitute/picard>) use the `FilterSamReads` function:

```sh
$ java -jar picard.jar FilterSamReads \
    I=Mt_alignment.bam \
    O=Mt_alignment_filtered.bam \
    READ_LIST_FILE=Readlist.txt \
    FILTER=includeReadList
 ```

3. The resulting output BAM file will contain only reads that have been identified by `Numt Parser` as being of putitive mt origin. Similar read filtering can be done in the program [`seqtk`](<https://github.com/lh3/seqtk>) or [`samtools`](<http://www.htslib.org/doc/samtools.html>).

It is imperative to mention that, while we have provided some filtering examples in this documentation, it is up to the user to decide a filtering criterion that is more suitable for their data. For example, when filtering paired-end data, users may elect to retain only reads pairs that are for which both reads are listed as having `cymt` origin, as described above. Alternatively, the users may choose to only discard individual reads of `numt` origin regardless of paired status or orientation, e.g., when data/endogenous reads are limited.

## Authors

Alida de Flamingh & Angel G. Rivera-Colon

University of Illinos at Urbana-Champaign
