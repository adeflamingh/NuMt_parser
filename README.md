# Numt Parser

*Numt Parser* a software for the detection and filtering of Nuclear Mitochondrial pseudogene (numt) contamination in Mitochondrial shotgun sequencing datasets.

## Installation and usage requirements

* Python >3.7

Other standard bioinformatic software are required for the processing and alignment of short-read sequencing data (e.g., `bwa`, `samtools`, etc.).

## Usage

### General data requirements

*Numt Parser* was designed to work on Illumina short-read sequencing data, most commonly focusing on the shotgun sequencing of mitochondrial genomes. Compatibility with other sequencing applications (e.g., long-read sequencing or RNAseq) has not been tested and thus we cannot guarantee the proper behavior of the software.

### Pre-*Numt Parser*

Before filtering, raw sequencing reads reads must be processed (trimmed to desired length and filtering) as desired.

Reads are then mapped using a short read aligner (such as `bwa` or `bowtie2`) to both a true mt sequence and a characterized numt reference.

Resulting alignments can then be additionally filtered (using tools such as `samtools`, `Picard`, etc.) to remove unmapped reads and low quality alignments. Final filtered alignments should be saved as `SAM` files for compatibility with *Numt Parser*.

### *Numt Parser* Analysis

The main *Numt Parser* executable, the `numt_parser.py` script, requires four (4) inputs files:

1. Cytoplasmic mitochondria reference file (FASTA)
2. numt reference file (FASTA)
3. Cytoplasmic mitochondria alignment file (SAM)
4. numt alignment file (SAM)

A fifth parameter specifying the location and name of the resulting output table (TSV) is also required.

#### Usage

```sh
numt_parser.py started on XXXX-XX-XX XX:XX:XX

usage: numt_parser.py [-h] --mt-fasta MT_FASTA --numt-fasta NUMT_FASTA --mt-sam MT_SAM --numt-sam
                      NUMT_SAM --outfile OUTFILE

Parse a set of mitochondrial reads and compare similarity between mitochondrial and NUMT references.

options:
  -h, --help                show this help message and exit
  --mt-fasta MT_FASTA       Mitochondrial Sequence FASTA
  --numt-fasta NUMT_FASTA   numt Sequence FASTA
  --mt-sam MT_SAM           Read alignments to the mt reference in SAM format
  --numt-sam NUMT_SAM       Read alignments to the numt reference in SAM format
  --outfile OUTFILE         Path and name to the output TSV file. Example: ./<sample_id>.tsv
```

#### Example

```sh
$ numt_parser.py \
    --mt-fasta /path/to/cymt.fasta \
    --numt-fasta /path/to/numt.fasta \
    --mt-sam /path/to/mt_alignment.sam \
    --numt-sam /path/to/numt_alignment.sam \
    --outfile /path/to/numt_parser_output.tsv
```

### *Numt Parser* output

The output of `numt_parser.py` is a table containing the identity statistics of each processed read.

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

## Citing *Numt Parser*

If you use the *Numt Parser* software for your research, please cite our *JoH* paper:

> de Flamingh A, Rivera-Colón AG, Gnoske TP, Kerbis Peterhans JC, Catchen JM, Malhi RS, Roca AL, ***Numt Parser*: Automated identification and removal of nuclear mitochondrial pseudogenes (numts) for accurate mitochondrial genome reconstruction in Panthera**, *Journal of Heredity*, Volume 114, Issue 2, March 2023, Pages 120–130, <https://doi.org/10.1093/jhered/esac065>

## On CIGAR strings

The Concise Idiosyncratic Gapped Alignment Report (CIGAR) string is a compact representation of an alignment [Li et al. 2009](https://doi.org/10.1093/bioinformatics/btp352). The CIGAR strings describes the status of each base in the alignment in form of several operations. These operations describe if a base in the alignment is "aligned" (i.e., matched, `M`), is part of an indel (`I` or `D` for insertions and deletions, respectively), or of it is clipped (`S` and `H` for soft and hard clips, respectively). Several other CIGAR operations are possible and are described in the [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

Since *Numt Parser* was designed to work on the alignments coming from the mapping of short-read sequencing data, we have limited its function to five CIGAR operations described by the most commonly used short-read aligners (e.g., `bwa`, `bowtie2`, `minimap2`). The five compatible operations are: `M`, `I`, `D`, `S`, and `H`. When *Numt Parser* detects an incompatible operation it will result in an error describing it. For example:

```sh
Error: 'P' is not a valid CIGAR opertation. Valid operations are: M,D,I,S,H. See README for more info.
```

This error describes the presence of a `P` or padding operation in a CIGAR string. The `P` operation can be commonly seen in multiple-sequence alignments (and in the alignment to a padded reference, as described in section 3 of the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)). In regards to *Numt Parser*, observing `P` operations might indicate that, instead of mapping the reads, the user multiple-sequence aligned the reads to genome. Similarly, the `N` operation is commonly used to represent introns in the alignment of RNAseq data and might describe that the user is mapping RNA instead of DNA sequences.

## Authors

**Alida de Flamingh**  
Carl R. Woese Institute for Genomic Biology  
University of Illinois at Urbana-Champaign 

**Angel G. Rivera-Colon**  
Department of Evolution, Ecology, and Behavior  
University of Illinos at Urbana-Champaign
