#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH -N 1
#SBATCH --mail-user=deflami2@illinoise.edu
#SBATCH --mail-type=ALL
#SBATCH -J aln_stat
# ----------------Load Modules--------------------
module load BEDTools/2.28.0-IGB-gcc-8.2.0    
module load   BCFtools/1.12-IGB-gcc-8.2.0  
module load SAMtools
# ----------------Commands------------------------
##FILES NEEDED
#bam files for sample
#reference genome fasta

#launch this code from the directory containing bams
sample_names=$(ls -1 | grep '^Singer' | cut -f1 -d '.' | sort -u)

for sample in $sample_names
do
    ##INDEX BAM FILES 
    #create .bai for each file
    samtools index ${sample}.bam
    
    ##CREATE BED FILE
    #bed file is used to mask low coverage regions of your mitogenome (after calling snps from it) 
    #specify a minimum coverage (here it is 3)
    bedtools genomecov -bga -ibam ${sample}.bam -g lion_ref.fasta | awk '$4 < 3' | bedtools merge -i - > ${sample}.lc_mask.bed
    
    ##CALL VARIANTS:    
    #The first mpileup part generates genotype likelihoods at each genomic position with coverage. Use the -Ou option when piping 
    #between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion. 
    #The second call part makes the actual calls. 
    #-c switch tells the program to use the consensus/biallelic caller, not the -m multiallelic default calling method
    #do not use -v option (asks to output only variant sites), all sites should be included 
    bcftools mpileup -Ou -f lion_ref.fasta ${sample}.bam | bcftools call -c -o ${sample}.calls.vcf.gz
    
    ##FILTER VARIANTS
    #use quality cut-off appropriate for your dataset
    bcftools view -i '%QUAL>=20' ${sample}.calls.vcf.gz -o ${sample}.filt.vcf.gz 
    
    ##INDEX VARIANT FILE
    bcftools index ${sample}.filt.vcf.gz
    
    ##CREATE CONSENSUS FASTA
    #cat reference, insert appropriate variants, mask areas of low coverage
    cat lion_ref.fasta | bcftools  consensus --mask ${sample}.lc_mask.bed ${sample}.filt.vcf.gz > ${sample}_consensus.fa 
done

