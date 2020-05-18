#!/bin/bash

src=/mnt/c/Users/AngelRC/Dropbox/lion_NUMT/data_angel/scripts
cd $src
# FASTA files
mt_fasta=/mnt/c/Users/AngelRC/Dropbox/lion_NUMT/data_angel/fasta_files/lion_ref.fasta
numt_fasta=/mnt/c/Users/AngelRC/Dropbox/lion_NUMT/data_angel/fasta_files/P_leo_KF907306_NUMT.fasta
# SAM files
mt_sam=/mnt/c/Users/AngelRC/Dropbox/lion_NUMT/data_angel/sam_testfiles/FMNH23970_real_15059reads.sam
numt_sam=/mnt/c/Users/AngelRC/Dropbox/lion_NUMT/data_angel/sam_testfiles/FMNH23970_NUMT_15666reads.sam
# Output file
out=/mnt/c/Users/AngelRC/Dropbox/lion_NUMT/data_angel/per_identity_stats/FMNH23970.identity_stats.tsv

# Run script
./compare_mt_numt.py \
    --mt-fasta $mt_fasta \
    --numt-fasta $numt_fasta \
    --mt-sam $mt_sam \
    --numt-sam $numt_sam \
    --outfile $out

