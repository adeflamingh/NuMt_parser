!#/usr/bin/end python

import pysam

fq= open('fastqHeader.txt').readlines()
fq= [x.strip() for x in fq]
fq= set(fq)

infile= pysam.AlignmentFile('in.bam')
outfile= pysam.AlignmentFile('out.bam', template= infile, mode= 'wb')
for aln in samfile:
    if aln.query_name in fq:
        outfile.write(aln)

infile.close()
outfile.close()
