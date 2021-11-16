#!/usr/bin/env python3
import argparse
PROG = 'parse_blast_table.py'
#
# Script to parse BLAST tabular format and subset reads based on top hits
#

# -----------
# Input files
# -----------

#
# Command Line Options
#
def parse_args():
    p = argparse.ArgumentParser(prog=PROG)
    p.add_argument('-b', '--blast-table', required=True, help='BLAST output table.')
    p.add_argument('-o', '--outfile',     required=True, help='Outfile path.')
    p.add_argument('-c', '--cymt-id',     required=True, help='Sequenece ID of the CYMT reference.')
    # Check input
    args = p.parse_args()
    assert os.path.exists(args.blast_table)
    return args

# ------------------------
# BLAST outfmt 6 structure
# ------------------------

#  1.   qseqid    query (e.g., unknown gene) sequence id
#  2.   sseqid    subject (e.g., ref genome) sequence id
#  3.   pident    percentage of identical matches
#  4.   length    alignment length (sequence overlap)
#  5.   mismatch  number of mismatches
#  6.   gapopen   number of gap openings
#  7.   qstart    start of alignment in query
#  8.   qend      end of alignment in query
#  9.   sstart    start of alignment in subject
#  10.  send      end of alignment in subject
#  11.  evalue    expect value
#  12.  bitscore  bit score

# -------
# Classes
# -------

#
# Keep elements of a blast hit
class BlastHit:
    def __init__(self, query_id, subject_id, pident, length, mismatch, evalue, bitscore):
        assert type(pident) in [int, float]
        assert type(length) in [int, float]
        assert type(mismatch) in [int, float]
        assert type(evalue) in [int, float]
        assert type(bitscore) in [int, float]
        self.quer = query_id
        self.subj = subject_id
        self.pidt = pident
        self.leng = length
        self.mist = mismatch
        self.eval = evalue
        self.bits = bitscore
    def __str__(self):
        return f'{self.quer} {self.subj} {self.pidt} {self.leng} {self.mist} {self.eval} {self.bits}'

# ---------
# Functions
# ---------

#
# Parse BLAST table and extract needed values
# Return a per-read dictionary of hits
def parse_blast_table(blast_file):
    # Check input
    # assert path.exist(blast_file) is True
    # Create output dictionary
    read_hit_dict = dict()
    # Open file
    for line in open(blast_file, 'r'):
        # Ignore empty lines or comments
        if len(line) == 0 or line[0] == '#':
            continue
        fields = line.strip('\n').split('\t')
        # Get the needed elements from the line
        query_id = fields[0]
        subject_id = fields[1]
        pident = float(fields[2])
        length = int(fields[3])
        mismatch = int(fields[4])
        evalue = float(fields[10])
        bitscore = float(fields[11])
        # Populate class
        blast_hit = BlastHit(query_id, subject_id, pident, length, mismatch, evalue, bitscore)
        # Set output dictionary
        read_hit_dict.setdefault(query_id, [])
        read_hit_dict[query_id].append(blast_hit)
    return read_hit_dict

#
# Function to find the top hit
def find_top_blast_hit(blast_hit_list):
    # Check inputs
    assert type(blast_hit_list) is list
    assert isinstance(blast_hit_list[0], BlastHit)
    # Output
    # This is a list in case there are more than one hit with the same top eval
    top_hits = list()
    # Find the top (lowest) evalue to use as comparison against all hits
    top_eval = min([ h.eval for h in blast_hit_list ])
    # Now, loop over the list and only extract hits that match the top (lower) evalue
    # In most case, there will only be a single top value
    for hit in blast_hit_list:
        if hit.eval <= top_eval:
            top_hits.append(hit)
    return top_hits


#
# Function to check alignments and print the corresponding CYMT reads from the top BLAST hit
def extract_top_hit_reads(read_hit_dict, cymt_id_name):
    # Check inputs
    assert type(read_hit_dict) is dict
    assert cymt_id_name is not None
    # CYMT read list output
    cymt_reads = []
    # Iterate over the dictionary
    for read_id in sorted(read_hit_dict.keys()):
        read_hits = read_hit_dict[read_id]
        assert isinstance(read_hits[0], BlastHit)
        # Find the top hit(s) in the list
        top_hits = find_top_blast_hit(read_hits)
        # Loop over the top hits and only retain those belonging to the CYMT reference
        for hit in top_hits:
            assert isinstance(hit, BlastHit)
            # This method includes reads if they have equal matches to CYMT and NUMT sequences.
            # In theory, if there is more than one top hit, one should belong to the CYMT and the other to the NUMT. Their e-values MUST be equal.
            if hit.subj == cymt_id_name:
                cymt_reads.append(hit.quer)
    return cymt_reads

#
# Function to print the output
def print_cymt_read_ids(cymt_reads_list, outfile):
    # Check inputs
    assert type(cymt_reads_list) is list
    out_f = open(outfile, 'w')
    # Loop over reads and save
    for read in cymt_reads_list:
        out_f.write(f'{read}\n')

#
# Main running function
def filter_blast_table(blast_table, outfile, cymt_id_name='KP202262.1_Ref_P_leo'):
    # 1. Read BLAST table and index all reads in dictionary
    read_hit_dict = parse_blast_table(blast_table)
    # 2. For each read, extract top hit and save read IDs
    cymt_reads_list = extract_top_hit_reads(read_hit_dict, cymt_id_name)
    # 3. Save output
    print_cymt_read_ids(cymt_reads_list, outfile)

# --------
# Run Code
# --------
args = parse_args()
filter_blast_table(args.blast_table, args.outfile, args.cymt_id)

