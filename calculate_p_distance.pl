#! /usr/bin/perl

# Takes in as arguments:
#	[0] one FASTA file containing one sequence;
#	[1] a second FASTA file containing a sequence aligned to the first.
# OUTPUTS: a p-distance AFTER removing positions containing gaps in both sequences.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# calculate_p_distance.pl <aligned_seq_1.fasta> <aligned_seq_2.fasta>
#########################################################################################

# Copyright (C) 2017 Chase W. Nelson
# Date created: November 11, 2017
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use warnings;


if(! $ARGV[1]) {
	die "\n### WARNING: Two fasta files must be supplied as arguments.\n\n";
}

my $fasta_file_1 = $ARGV[0];
my $fasta_file_2 = $ARGV[1];

my $seq1;
my $seq2;

open(FASTA_FILE_1,"$fasta_file_1"); # First FASTA sequence
while (<FASTA_FILE_1>) { # For each line of the FASTA
	unless($_ =~/^>/) {
		chomp;
		$seq1 .= $_;
	}
}
close FASTA_FILE_1;

open(FASTA_FILE_2,"$fasta_file_2"); # First FASTA sequence
while (<FASTA_FILE_2>) { # For each line of the FASTA
	unless($_ =~/^>/) {
		chomp;
		$seq2 .= $_;
	}
}
close FASTA_FILE_2;

my $p_distance = &calculate_p_distance($seq1, $seq2);

print "p=$p_distance\n\n";


#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################


#########################################################################################
sub calculate_p_distance {
	my ($seq1,$seq2)=@_;
	
	my $seq_length = length($seq1);
	
	print "\nOriginal sequence length: $seq_length\n";
	
#	my @indices_both_gaps;
	my %indices_both_gaps;
	
	# First, find positions that contain gaps in both and eliminate them
	for(my $i=0; $i<$seq_length; $i++) {
		
		# If they're both gaps
		if(substr($seq1, $i, 1) eq '-' && substr($seq2, $i, 1) eq '-') {
#			push(@indices_both_gaps, $i);
			$indices_both_gaps{$i}=1;
		}
		
		# If they're both N
		if(substr($seq1, $i, 1) eq 'N' && substr($seq2, $i, 1) eq 'N') {
#			push(@indices_both_gaps, $i);
			$indices_both_gaps{$i}=1;
		}
		
	}
	
	# Cut out positions with gaps in both sequences
	my $running_num_deletions = 0;
	
	for(my $i=0; $i<$seq_length; $i++) {
		if(exists $indices_both_gaps{$i}) {
			if($indices_both_gaps{$i} == 1) {
				substr($seq1, ($i - $running_num_deletions), 1, '');
				substr($seq2, ($i - $running_num_deletions), 1, '');
				$running_num_deletions++;
			}
		}
	}
	
	my $new_seq1_length = length ($seq1);
	my $new_seq2_length = length ($seq2);
	
	print "Sequence length after excluding gaps: $new_seq1_length\n";
	
	if($new_seq1_length != $new_seq2_length) {
		die "### SEQUENCES DIFFER IN LENGTH AFTER GAP PROCESSES: Chase's fault. Email him.\n\n";
	}
	
	# We find out new positions that match (1) and don't (0)
	my $match_string = &match_two_seqs($seq1,$seq2); # these are the new seqs with deletions
	#print "$match_string\n";
	# RESULT: 111111110010111010000000000101110000000000111111000011010010110000110110110010110111111110111010010101110110110111110110010000110110110010111110010000000010110011111110110010010110111010100110111111111010111111110110110110010110011100110110111010010111000110110010110110111010110111000010110110110110100000111010000000111000110110010110000011110110110110110110110010111111110110110110111110010110110110000110111110110110000111100111010100111100110010010010010101001010100111011110011110010010100110001011010110110110110111111010111111110110110110010111110110110111111111110110110111110110111110110110110110010110111110110010111110010110011110011110000000010000000000000000000000000101
	
	my $match_string_length = length($match_string);
	
	my $num_matching_positions;
	for(my $i=0; $i<$match_string_length; $i++) {
		if(substr($match_string, $i, 1) == 1) {
			$num_matching_positions++;
		}
	}
	
	my $num_mismatch_positions = $match_string_length - $num_matching_positions;
	my $p_distance = ($num_mismatch_positions / $match_string_length);
	
	return $p_distance;
}


#########################################################################################
sub match_two_seqs {
	my ($seq1,$seq2)=@_;
	
	my $match_string = '';
	
	if(length($seq1) != length($seq2)) {
		$match_string = 'ERROR';
	} else {
		while(length($seq1)>0 && length($seq2)>0) {
			my $seq1_char = substr($seq1,0,1,'');
			my $seq2_char = substr($seq2,0,1,'');
			
			if($seq1_char eq $seq2_char) {
				$match_string .= '1';
			} else {
				$match_string .= '0';
			}
		}
	}
	
	return $match_string;
}