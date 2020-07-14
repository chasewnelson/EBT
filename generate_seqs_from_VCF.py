#! /usr/bin/env python

###############################################################################
## LICENSE
##
## Copyright (C) 2020
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# DATE: 2020
# AUTHOR: Chase W. Nelson
# AFFILIATION: Biodiversity Research Center, Academia Sinica, Taipei, Taiwan
# CONTACT: cnelson@gate.sinica.edu.tw
# FUNDING: Academia Sinica Postdoctoral Research Fellowship, P.I. Wen-Hsiung Li
# CITATION: Nelson CW, Evolutionary Bioinformatics Toolkit, https://github.com/chasewnelson/EBT

from Bio import SeqIO
import os
import random
import re
import sys

# Usage = "Script to generate a FASTA with randomly interspersed variants (from VCF)"

# Example:
# generate_seqs_from_VCF.py reference.fasta variants.vcf 1000

# Check argument number
if not len(sys.argv) == 4:
    raise SystemExit("\n### TERMINATED: need 3 unnamed arguments:\n" + \
                     "    # (1) FASTA reference (1 sequence)\n" + \
                     "    # (2) .VCF SNP report\n" + \
                     "    # (3) number of sequences to generate\n\n")

infile_FASTA = sys.argv[1]
infile_VCF = sys.argv[2]
num_seqs = int(sys.argv[3])

### Gather VCF files in working directory
#VCF_filenames = [name for name in os.listdir(".") if os.path.isfile(name) and name.endswith(".vcf")]

print("FASTA reference: " + str(infile_FASTA))
print("VCF with variants: " + str(infile_VCF))
print("Num sequences to generate: " + str(num_seqs))
#print("VCF files to examine, gathered from working directory: " + str(VCF_filenames))

# Check infiles exist
if not os.path.isfile(infile_FASTA):
    raise SystemExit("\n### TERMINATED: infile 1 (FASTA reference) does not exist\n")

if not os.path.isfile(infile_VCF):
    raise SystemExit("\n### TERMINATED: infile 2 (.VCF SNP report) does not exist\n")

# Create outfile_prefix name
outfile_prefix = str(infile_VCF)
outfile_prefix = outfile_prefix.replace(".vcf", "") # TODO: learn regex
outfile_prefix = outfile_prefix.replace(".fasta", "") # TODO: learn regex
outfile_prefix = outfile_prefix.replace(".txt", "") # TODO: learn regex
outfile_prefix = outfile_prefix.replace(".tsv", "") # TODO: learn regex
outfile_name = outfile_prefix + "_nSeqs" + str(num_seqs) + ".fasta"
print("Will write output to: " + outfile_name + "\n")


###############################################################################
### Extract SNP data from VCF

# Build a dictionary of site(int)->nt(str)->freqs(list)
# This will later be used to find mean frequency across all samples
site_nt_freqs = {}
site_nt_counts = {}

# Prepare regex
#VCF_line_pattern = r"([\w\_\.]+)\t(\d+)\t([^\s]+)\t([acgtuACGTU]+)\t([acgtuACGTU]+)\t(\d+)\t(\w+)\t([\w\=\;\.\,]+)"
#VCF_INFO_regex = re.compile(VCF_INFO_pattern)
#VCF_DP_pattern = r"DP=\d+" # DP=266
#VCF_DP_regex = re.compile(VCF_DP_pattern)
VCF_AF_pattern = r"AF=[\d\.]+" # AF=0.015038
VCF_AF_regex = re.compile(VCF_AF_pattern)
#VCF_SB_pattern = r"SB=\d+" # SB=0
#VCF_SB_regex = re.compile(VCF_SB_pattern)
#VCF_DP4_grppattern = r"DP4=(\d+),(\d+),(\d+),(\d+)" # DP4=27,235,0,4
#VCF_DP4_grpregex = re.compile(VCF_DP4_pattern)
#sample_count = 0


#for this_VCF in VCF_filenames:
if os.path.isfile(infile_VCF): # this_VCF
    with open(infile_VCF, "r") as f:
        lines = (line.rstrip() for line in f)
        this_VCF_ID = str(infile_VCF).rstrip(".vcf")
        #sample_count += 1
        variant_count = 0

        for line in lines:
            if not line.startswith("#"):  # skip headers
                variant_count += 1
                line_list = line.split("\t")

                # Extract information for this SNP
                # this_CHROM = line_list[0]
                this_POS = int(line_list[1])
                # this_ID = line_list[2]
                # this_REF = line_list[3]
                this_ALT = str(line_list[4])
                # this_QUAL = line_list[5]
                # this_FILTER = line_list[6]
                this_INFO = line_list[7]
                # this_INFO_regexgroups = VCF_INFO_regex.match(this_INFO)
                #this_DP = VCF_DP_regex.match(this_INFO)
                #this_DP.replace("DP=", "")
                this_AF = VCF_AF_regex.search(this_INFO)
                this_AF = this_INFO[this_AF.start():this_AF.end()]  # AF=0.367003
                this_AF = float(this_AF.replace("AF=", ""))
                # this_SB = VCF_SB_regex.match(this_INFO)
                # this_DP4_grps = VCF_DP4_grpregex.match(this_INFO)

                # Add data for this record
                if this_POS not in site_nt_freqs.keys(): # this_POS not in site_nt_freqs:
                    site_nt_freqs[this_POS] = {}
                    site_nt_counts[this_POS] = {}
                    site_nt_freqs[this_POS][this_ALT] = this_AF
                    site_nt_counts[this_POS][this_ALT] = 1
                elif this_ALT not in site_nt_freqs[this_POS]:
                    site_nt_freqs[this_POS][this_ALT] = {}
                    site_nt_counts[this_POS][this_ALT] = {}
                    site_nt_freqs[this_POS][this_ALT] = this_AF
                    site_nt_counts[this_POS][this_ALT] = 1
                else:
                    site_nt_freqs[this_POS][this_ALT] += this_AF
                    site_nt_counts[this_POS][this_ALT] += 1

else:
    raise SystemExit("\n### TERMINATED: file doesn't exist: " + str(infile_VCF) + "\n")

###############################################################################
### Record reference sequence
# Open FASTA for reading
rec = SeqIO.read(infile_FASTA, "fasta") # will throw error if more than one record!
ref_seq = str(rec.seq)
ref_seq_len = len(ref_seq)


###############################################################################
### Build lists of nucleotides to observe at each site, from which we will draw without replacement
site_nt_lists = {}
for this_site in site_nt_freqs.keys():
    this_nt_REF = ref_seq[this_site - 1]

    for this_nt_ALT in site_nt_freqs[this_site].keys():
        this_nt_ALT_freq = float(site_nt_freqs[this_site][this_nt_ALT]) / float(site_nt_counts[this_site][this_nt_ALT])
        num_ALT = int(round(this_nt_ALT_freq * num_seqs)) # rounded to nearest int
        num_REF = num_seqs - num_ALT

        if num_REF < 0:
            raise SystemExit("\n### TERMINATED: negative number of nucleotides implied at site: " + str(this_site) + "\n")
        else:
            # Add reference nucleotides
            for ref_i in range(num_REF):
                if this_site in site_nt_lists.keys():
                    site_nt_lists[this_site].append(str(this_nt_REF))
                else:
                    site_nt_lists[this_site] = [str(this_nt_REF)]

            # Add alternate nucleotides
            for alt_i in range(num_ALT):
                if this_site in site_nt_lists:
                    site_nt_lists[this_site].append(str(this_nt_ALT))
                else:
                    site_nt_lists[this_site] = [str(this_nt_ALT)]


###############################################################################
### Randomize the order of nucleotides at each site
this_seed = random.randrange(sys.maxsize)
random.seed(this_seed)
print("\nRandom number seed for nucleotide incorporation: " + str(this_seed) + "\n")
for this_site in site_nt_lists.keys():
    random.shuffle(site_nt_lists[this_site])


###############################################################################
### Build FASTA alignment with variants!

# Open FASTA for writing
outfile_hdl = open(outfile_name, "w")

for this_seq_idx in range(num_seqs):
    this_seq = str(ref_seq)
    this_seq_ID = str(outfile_prefix) + "_pseudoseq_" + str(this_seq_idx + 1)

    # Make sure we used all the nucleotides

    for this_site_idx in range(len(this_seq)):
        this_site_num = this_site_idx + 1

        if this_site_num in site_nt_lists:
            this_seq = this_seq[:this_site_idx] + site_nt_lists[this_site_num].pop() + this_seq[(this_site_idx + 1):]
            #this_seq[this_site_idx] = site_nt_lists[this_site_num].pop() # strings immutable
            # else it's the REF for sure; leave it

    # Write to FASTA file
    outfile_hdl.write(">" + this_seq_ID + "\n")
    outfile_hdl.write(str(this_seq) + "\n")

# Close FASTA output file
outfile_hdl.close()

# Make sure we used all the nucleotides were used
for this_seq_idx in range(num_seqs):
    this_site_num = this_seq_idx + 1

    if this_site_num in site_nt_lists and len(site_nt_lists[this_site_num]) > 0:
        raise SystemExit("\n### WARNING: sequences remain at site " + str(this_site_num) + ". BIG PROBLEM.\n")


###############################################################################
### FINISH
print("\n")
raise SystemExit("\n### All done; we stopped here, dear researcher.\n")

