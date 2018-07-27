#! /usr/bin/perl

# Takes in two arguments: (1) aligned amino acid sequence; and (2) nucleotide sequence.
#	Prints the codon-aligned nucleotide sequence to screen.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# align_codon2aa.pl MSAARKRTL-L ATGTCGGCGGCTCGTAAGCGCACCTTGTTG
#########################################################################################

# Copyright (C) 2018 Chase W. Nelson
# Date created: July 27, 2018
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use warnings;

# This example should NOT work
#my $this_aln_aa_seq = 'MSAARKRTLLKVIILGDSGVGKTSLINQYVNNKFSNQHKATIGADFLTKEVVEDRLVTMQIWDTAGQERLQGLGVAFYRGADCCVLVYDVYVMKSFDNLDYWRDEFVIQAAPSDQEYFPFVVLGNKVDVDGGNSRVVSEKKAKAWCAAKGGIPYFESSAKEDFNVDAAFQCIAKNALKNETEEEIYLPDTIDVNASRPQKTSGCEC';
#my $this_nt_seq = 'ATGTCGGCGGCTCGTAAGCGCACCTTGTTGAAGGTCATCATCCTCGGCGATAGCGGGGTTGGAAAGACATCGCTGATAAACCAATATGTTAATAACAAGTTCAGCAATCAACACAAGGCCACTATTGGAGCTGATTTCCTAACTAAAGAAGTATAAGTTGAAGACAGACTTGTTACAATGCAGATTTGGGATACGGCTGGACAAGAACGGCTTCAGGGTCTCGGTGTTGCCTTCTACAGAGGTGCAGATTGTTGTGTTCTTGTTTACGATGTTTATGTCATGAAGTCATTCGATAACCTTGACTACTGGAGGGATGAGTTTGTGATCCAGGCAGCCCCATCCGATCAGGAGTACTTCCCCTTTGTAGTACTTGGAAATAAAGTGGACGTGGATGGTGGCAATAGTCGAGTGGTCTCCGAGAAGAAGGCCAAAGCATGGTGTGCAGCGAAAGGAGGCATCCCCTACTTTGAGTCATCAGCCAAGGAAGACTTCAACGTGGATGCTGCATTCCAGTGTATTGCCAAGAACGCATTGAAGAACGAGACGGAGGAGGAAATCTACCTGCCTGATACGATCGACGTGAACGCCAGCAGGCCACAGAAAACTTCCGGATGCGAGTGTTAA';

# This example SHOULD work
#my $this_aln_aa_seq = 'MSAARKRTL-L';
#my $this_nt_seq = 'ATGTCGGCGGCTCGTAAGCGCACCTTGTTG'

my $this_aln_aa_seq = $ARGV[0];
my $this_nt_seq = $ARGV[1];

print "\n\n##########################################################################################\n";
print "### Performing codon alignment...\n";

unless($this_nt_seq) { 
	die "\n## Must submit the following arguments (in order):\n".
		"##\t(1) aligned amino acid sequence\n##\t(2) nucleotide sequence\n## TERMINATED.\n\n"; 
}

my $this_nt_seq_aligned = &align_codon2aa($this_aln_aa_seq, $this_nt_seq);

print "\n### CODON ALIGNMENT: $this_nt_seq_aligned\n";

exit;

#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################


##########################################################################################
sub align_codon2aa {
	my ($curr_aa_seq, $curr_nt_seq) = @_;
	
	# RECORD positions of GAPS in the amino acid, and also check for correct lengths
	# Store gap indices for each seq in an array of arrays
	my @aa_gap_indices;

	my $curr_aa_seq_length = length($curr_aa_seq);
	my $ungapped_curr_aa_seq_length = $curr_aa_seq_length;
	my $nt_seq_length = length($curr_nt_seq);
		
	for(my $pos_index = 0; $pos_index < $curr_aa_seq_length; $pos_index++) {
		my $curr_aa = substr($curr_aa_seq, $pos_index, 1);
		
		if($curr_aa =~ /-/) { # It's a gap
			$ungapped_curr_aa_seq_length -= 1;
			push(@aa_gap_indices, $pos_index);
		}
	}
	
	if(($ungapped_curr_aa_seq_length * 3) != $nt_seq_length) {
		# If nucleotide sequence contains a STOP codon not represented in the amino acid
		# sequence or vice versa, modify the nucleotide sequence
		my $last_codon = substr($curr_nt_seq, -3);
		my $last_aa = substr($curr_aa_seq, -1);
		
		#print "\nlast_codon = $last_codon\n";
		
		if(($last_codon eq 'TAA' || $last_codon eq 'TAG' || $last_codon eq 'TGA') && ($last_aa ne '*')) {
			# Remove the STOP codon to nt seq
			#warn "\n### STOP CODON REMOVED.\n";
			print "\n### REMOVING STOP codon:\nBEFOR: $curr_nt_seq\n";
			substr($curr_nt_seq, -3, 3, '');
			print "AFTER: $curr_nt_seq\n";
			
			
		} elsif(($last_codon ne 'TAA' && $last_codon ne 'TAG' && $last_codon ne 'TGA') && ($last_aa eq '*')) {
			# Add the STOP codon to nt seq
			#warn "\n### STOP CODON ADDED.\n";
			print "\n### ADDING STOP codon:\nBEFOR: $curr_nt_seq\n";
			$curr_nt_seq .= 'TAA';
			print "AFTER: $curr_nt_seq\n";
		}
		
		$nt_seq_length = length($curr_nt_seq);
		#print "\nungapped_curr_aa_seq_length X 3 = " . $ungapped_curr_aa_seq_length * 3 . "\n" .
		#	"nt_seq_length = $nt_seq_length\n";
		
		
		# If that didn't fix it
		if(($ungapped_curr_aa_seq_length * 3) != $nt_seq_length) {
			return '<ERRONEOUS INPUT; SEQUENCE LENGTHS NOT CONSISTENT>'; # DONE
		}
	}

	# INSERT the appropriate gaps into the nucleotide sequences
	if(scalar(@aa_gap_indices) >= 1) {
		
		foreach(@aa_gap_indices) {
			#print "$_ and ";
			substr($curr_nt_seq, ($_ * 3), 0, '---')
		}
		
		#print "that's it.\n";
		
		#print "New nt sequence length of <THIS> is " . length($curr_nt_seq) . "\n";
		
	}
	
	# WRITE the new nucleotide sequence with amino acid alignment imposed
	my $curr_length = length($curr_nt_seq);
	
	if ($curr_length < 3) {
		die "\n\n### WARNING: Must be at least 3 nucleotides.\n\n";
	}
	
	if(($curr_length % 3) > 0) {
		warn "\n\n### WARNING: Length is not a multiple of 3 (complete codons).\n\n";
	}
	
	#print "len is ".length($curr_nt_seq)."\n$curr_nt_seq";
	
	#print "\nAligned Nucleotide Sequence:\n$curr_seq\n\n";

	return "$curr_nt_seq";

}
