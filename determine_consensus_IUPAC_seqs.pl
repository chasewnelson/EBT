#! /usr/bin/perl

# INPUT: Takes in one FASTA file with multiple aligned NUCLEOTIDE seqs as an argument.
# OUTPUTS: two files to the working directory, one with a consensus sequence 
#	(*_consensus.fa) and one with an IUPAC sequence (_IUPAC.fa); and brief summary 
#	statistics to the Terminal.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# determine_consensus_IUPAC_seqs.pl <aligned_seqs.fasta>
#########################################################################################

# Copyright (C) 2016 Chase W. Nelson
# Date created: November 2016
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use warnings;

my $filename = $ARGV[0];

# Read in the sequence from the file
my $seq = '';
my @seqs_arr;
#my $header = '';
#my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN, "$filename") or die "Could not open aligned FASTA file $filename\n";

print "\n\nReading in FASTA sequences from files...\n\n";

while(<IN>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			#$header = $_;
			$seq_num ++;
		} else {
			push(@seqs_arr,$seq);
			#push(@headers_arr,$header);
			#$header = $_;
			$seq_num ++;
			
			my $this_seq_length = length($seq);
			
			#print "\nseq $seq_num is of length $this_seq_length\n";
			
			if($last_seq_length && ($last_seq_length != $this_seq_length)) {
				die "\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
			} else {
				$last_seq_length = $this_seq_length;
				#print "\nseq: $seq\n";
				$seq = '';
			}
		}
	} else {
		$seq .= $_;
	}
}

push(@seqs_arr,$seq);
#push(@headers_arr,$header);

close IN;

print "\n\nTabulating nucleotides...\n\n";

$seq_num = scalar(@seqs_arr);

my @A_counts;
my @C_counts;
my @G_counts;
my @T_counts;

# INITIALIZE ARRAY VALUES
for(my $site_id=1; $site_id<=length($seq); $site_id++) { # for each site
#	my $curr_site_index = $site_id - 1;
	push(@A_counts,0);
	push(@C_counts,0);
	push(@G_counts,0);
	push(@T_counts,0);
}

for(my $seq_id=1; $seq_id <= $seq_num; $seq_id++) { # for each sequence
	my $curr_seq_index = $seq_id - 1;
	my $curr_seq = $seqs_arr[$curr_seq_index];
	#print "\ncurr_seq: $curr_seq\n";
	
	for(my $site_id=1; $site_id<=length($curr_seq); $site_id++) { # for each site
		my $curr_site_index = $site_id - 1;
		my $curr_nt = substr($curr_seq, $curr_site_index, 1); # substr is 0-based indexing. VAR,OFFSET,LEN
		#print "\ncurr_nt: $curr_nt";
		
		if($curr_nt eq 'A') {
			$A_counts[$curr_site_index]++;
			#$C_counts[$curr_site_index]+=0;
			#$G_counts[$curr_site_index]+=0;
			#$T_counts[$curr_site_index]+=0;
		} elsif($curr_nt eq 'C') {
			#$A_counts[$curr_site_index]+=0;
			$C_counts[$curr_site_index]++;
			#$G_counts[$curr_site_index]+=0;
			#$T_counts[$curr_site_index]+=0;
		} elsif($curr_nt eq 'G') {
			#$A_counts[$curr_site_index]+=0;
			#$C_counts[$curr_site_index]+=0;
			$G_counts[$curr_site_index]++;
			#$T_counts[$curr_site_index]+=0;
		} elsif($curr_nt eq 'T') {
			#$A_counts[$curr_site_index]+=0;
			#$C_counts[$curr_site_index]+=0;
			#$G_counts[$curr_site_index]+=0;
			$T_counts[$curr_site_index]++;
		}
	}	
}

print "\n\nGenerating consensus and IUPAC files...\n\n";

# Do calculations and output to summary file
my $consensus_file_name;
if($filename =~ '.fasta') {
	$consensus_file_name = $` . "_consensus.fa";
} elsif($filename =~ '.fa') {
	$consensus_file_name = $` . "_consensus.fa";
} else {
	$consensus_file_name = "fasta_consensus.fa";
}

my $IUPAC_file_name;
if($filename =~ '.fasta') {
	$IUPAC_file_name = $` . "_IUPAC.fa";
} elsif($filename =~ '.fa') {
	$IUPAC_file_name = $` . "_IUPAC.fa";
} else {
	$IUPAC_file_name = "fasta_IUPAC.fa";
}

open(OUT_CONSENSUS, ">>$consensus_file_name");
print OUT_CONSENSUS "\>$filename\_CONSENSUS\n";

open(OUT_IUPAC, ">>$IUPAC_file_name");
print OUT_IUPAC "\>$filename\_IUPAC\n";

my $nt_counter = 0;
my $total_aligned_sites = scalar(@A_counts);
my $determinate_sites = 0;
my $determinate_poly_sites = 0;
my $indeterminate_sites = 0;

for(my $site_id=1; $site_id<=$total_aligned_sites; $site_id++) { # for each site
	my $site_index = $site_id - 1;
	
	# The following method should automatically ignore gaps (-)
	my $A = $A_counts[$site_index];
	my $C = $C_counts[$site_index];
	my $G = $G_counts[$site_index];
	my $T = $T_counts[$site_index];
	
	my $total = ($A + $C + $G + $T);
	
#	print "My total is $total\n";
	
#	if($total != $seq_num) {
#		#die "\n\nNucleotide sum ($total) at site $site_id does not match number of sequences ($seq_num).\n".
#		#	"A=$A\nC=$C\nG=$G\nT=$T\nTERMINATED\n\n";
#		
##		print "\n\nNucleotide sum ($total) at site $site_id does not match number of sequences ($seq_num).\n".
##			"A=$A\nC=$C\nG=$G\nT=$T\nALIGNMENT GAP?\n\n";
#	}
	
	if($total > 0) { # there are some determinate nucleotides
		$determinate_sites++;
		# Which is the dominant nucleotide?
		my $maj_nt; # a variant nucleotide may have fixed
		my $maj_non_N_nt;
		my $maj_nt_count = 0;
		if($A > $maj_nt_count) {
			$maj_nt_count = $A;
			$maj_nt = 'A';
			$maj_non_N_nt = 'A';
		} elsif($A == $maj_nt_count) {
			$maj_nt = 'N';
		}
		if($C > $maj_nt_count) {
			$maj_nt_count = $C;
			$maj_nt = 'C';
			$maj_non_N_nt = 'C';
		} elsif($C == $maj_nt_count) {
			$maj_nt = 'N';
		}
		if($G > $maj_nt_count) {
			$maj_nt_count = $G;
			$maj_nt = 'G';
			$maj_non_N_nt = 'G';
		} elsif($G == $maj_nt_count) {
			$maj_nt = 'N';
		}
		if($T > $maj_nt_count) {
			$maj_nt_count = $T;
			$maj_nt = 'T';
			$maj_non_N_nt = 'T';
		} elsif($T == $maj_nt_count) {
			$maj_nt = 'N';
		}
		
		my $min_sum = $total - $maj_nt_count;
		
#		if($min_sum > 1) {
#			$non_singleton_sites{$site_id} = 1;
#		}
		
		# Go through variable sites and determine IUPAC code
		# A, C, G, T as expected
		# Y: pYrimidine; C or T
		# R: puRine; A or G
		# S: Strong bond: C or G
		# W: Weak bond: A or T
		# K: Keto: G or T [done]
		# M: aMino: A or C [done]
		# B: all except A, after which comes B: C, G, T [done]
		# D: all except C, after which comes D: A, G, T [done]
		# H: all except G, after which comes H: A, C, T [done]
		# V: all except T/U, after which comes V: A, C, G [done]
		# N: aNy base [done]
		
		my $IUPAC_nt = $maj_nt; # initialize as consensus
		
		if($A && $C && $G && $T) {  # N
		
			$determinate_poly_sites++;
			$IUPAC_nt = 'N';
			
		} elsif($A && $C && $G) { # V
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'V';
			
		} elsif($A && $C && $T) { # H
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'H';
			
		} elsif($A && $G && $T) { # D
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'D';
			
		} elsif($C && $G && $T) { # B
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'B';
			
		} elsif($A && $C) { # M
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'M';
			
		} elsif($G && $T) { # K
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'K';
			
		} elsif($A && $T) { # W
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'W';
			
		} elsif($C && $G) { # S
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'S';
			
		} elsif($A && $G) { # R
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'R';
			
		} elsif($C && $T) { # Y
			
			$determinate_poly_sites++;
			$IUPAC_nt = 'Y';
			
		} 
		
		if($nt_counter == 60) {
#			print OUT_CONSENSUS "$maj_nt\n";
			print OUT_CONSENSUS "$maj_non_N_nt\n";
			print OUT_IUPAC "$IUPAC_nt\n";
			$nt_counter = 0;
		} else {
#			print OUT_CONSENSUS "$maj_nt";
			print OUT_CONSENSUS "$maj_non_N_nt";
			print OUT_IUPAC "$IUPAC_nt";
			$nt_counter++; 
		}
		
	} else { # there are no determinate nucleotides
		$indeterminate_sites++;
		if($nt_counter == 60) {
			print OUT_CONSENSUS "N\n";
			print OUT_IUPAC "N\n";
			$nt_counter = 0;
		} else {
			print OUT_CONSENSUS "N";
			print OUT_IUPAC "N";
			$nt_counter++; 
		}
	}
}

close OUT_CONSENSUS;
close OUT_IUPAC;

print "SUMMARY:\nNum. sequences examined: $seq_num\n".
	"Num. determinate sites: $determinate_sites\n".
	"Num. determinate sites polymorphic: $determinate_poly_sites\n".
	"Num. indeterminate sites: $indeterminate_sites\n".
	"Sum determinate and indeterminate: ".($determinate_sites+$indeterminate_sites)."\n".
	"Total aligned sites: $total_aligned_sites\n\n";

print "\n******\nCOMPLETED.\n******\n\n";
