#! /usr/bin/perl

# Takes in as arguments:
#	[0] one FASTA file containing MULTIPLE ALIGNED SEQUENCES to be stored by the first
#	string to follow the '>' symbol, until the first space--that is, the "ID"; and
#	[1] a file with one list (one column) of ID's to output.
# OUTPUTS: a fasta file containing ONLY the sequences with headers matching the IDs given
#   in the file given in argument [1] (second argument).

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# store_fasta_by_ID.pl all_my_sequences.fasta just_wanted_seqs_headers.txt > just_wanted_seqs.fa
#
# where <all_my_sequences.fasta> might begin:
# >My_FASTA_Header1
# ATGCCAGGGTTT...
# >My_FASTA_Header2
# ATGCCAGGGTTC...
# ...
#
# and where <ust_wanted_seqs_headers.txt> might begin:
# My_FASTA_Header2
# My_FASTA_Header6
# ...
#########################################################################################

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

# We have co-opted by natural selection. That's what this is.

use strict;
#use warnings;
use Data::Dumper;
use List::Util qw(max);

my $fasta_file_name = $ARGV[0];
my $list_ID_file_name = $ARGV[1];
unless($fasta_file_name =~ /.fa/) { die "\n\n# FASTA file must contain .fa or .fasta extension. TERMINATED\n\n"; }
unless($list_ID_file_name =~ /.txt/) { die "\n\n# File with ID list file must contain .txt extension. TERMINATED\n\n"; }

# Later, flesh this out, and in the next file read check to see if it's present, to save
# time; for the moment, quick and dirty approach
#my %ids_present;
#
#open(IN_IDS, "$list_ID_file_name") or die "Could not open file $list_ID_file_name\n";
#
#print "\nReading in sequences IDs to select from $list_ID_file_name...\n";
#
#while(<IN_IDS>) {
#	chomp;
#	my $curr_id = $_;
#	$ids_present{$curr_id} = 1;
#}

# Read in the group of sequences from the fasta file
my $seq = '';
my %seqs_hash;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file_name") or die "Could not open file $fasta_file_name\n";

#print "\nRecording coding sequence data for $fasta_file_name...\n";

my $curr_header = '';
my $curr_seq_ID = '';

# DURING THIS STEP, let's automate choosing the sequence with fewest N's
while(<IN_FASTA>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$curr_header = $_;
			
#			if($curr_header =~/>([\w\-\:']+)/) { # USED FOR HPV # whatever follows '>' until first underscore or space # origin />([\w\.\-\:']+)/
			if($curr_header =~/>([\w\-\:\.\_']+)/) { # whatever follows '>' until first underscore or space # origin />([\w\.\-\:']+)/
				$curr_seq_ID = $1;
#				print "\nMatched $curr_seq_ID\n";
			}
			
			$seq_num ++;
		} else {
			# READY TO STORE
			# If there are multiples, it'll contain an underscore
			## HPV ONLY
			##if($curr_seq_ID =~ /_/) {
			##	$curr_seq_ID = $`;
			##}
			
			# Store sequence
			if($seqs_hash{$curr_seq_ID}) { # if MATCHES ANOTHER previous...
				my $this_N_count = () = $seq =~ /N/g;
				my $other_N_count = () = $seqs_hash{$curr_seq_ID} =~ /N/g;
				
				if($this_N_count < $other_N_count) { # ... choose one with fewest N's
					$seqs_hash{$curr_seq_ID} = $seq;
				} # else just keep it
				
			} else { # just store it
				$seqs_hash{$curr_seq_ID} = $seq;
			}
			
			# NEW HEADER & SEQUENCE!
			$curr_header = $_;
			
#			if($curr_header =~/>([\w\-\:']+)/) { # whatever follows '>' until first space # origin />([\w\.\-\:']+)/
			if($curr_header =~/>([\w\-\:\.\_']+)/) { # whatever follows '>' until first space # origin />([\w\.\-\:']+)/
				$curr_seq_ID = $1;
#				print "\nMatched $curr_seq_ID\n";
			}
			
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

close IN_FASTA;

# READY TO STORE FINAL SEQUENCE
# If there are multiples, it'll contain an underscore 
## HPV ONLY
##if($curr_seq_ID =~ /_/) {
##	$curr_seq_ID = $`;
##}

# Store sequence
if($seqs_hash{$curr_seq_ID}) { # if MATCHES ANOTHER previous...
	my $this_N_count = () = $seq =~ /N/g;
	my $other_N_count = () = $seqs_hash{$curr_seq_ID} =~ /N/g;
	
	if($this_N_count < $other_N_count) { # ... choose one with fewest N's
		$seqs_hash{$curr_seq_ID} = $seq;
	} # else just keep it
	
} else { # just store it
	$seqs_hash{$curr_seq_ID} = $seq;
}
			
# Simpler before			
#$seqs_hash{$curr_seq_ID} = $seq;
#push(@headers_arr,$header);

#print "\n";

my %ids_present;

# Finds IDS to print
open(IN_IDS, "$list_ID_file_name") or die "Could not open file $list_ID_file_name\n";

#print "\nReading in sequences IDs to select from $list_ID_file_name...\n";

# THIS FILE CONTAINS ONLY ID's, NO UNDERSCORES (_)
while(<IN_IDS>) {
	chomp;
	my $curr_id = $_;
	
	$ids_present{$curr_id} = 1;
	print ">$curr_id\n" . $seqs_hash{$curr_id} . "\n";
}

# This is for choosing all NON-INDICATED ID's... but could just change input ID file and use above
###### BEFORE ######
# THIS WILL PRINT every sequence that is in the original fasta but NOT in the list of IDs
# Thus, if there are inexact IDs, they WILL BE PRINTED. So we will have multiples in this
# file of some sequences which shouldn't be present at all; namely, things with "_" in them.
# keys %seqs_hash) { # these are the EXACT ID's in the ORIGINAL FASTA
###### NOW ######
# We have already trimmed the underscore (_) and what follows from the sequence ID
# Thus, we have already eliminated redundancies and chosen the sequences with fewest N's

#foreach(keys %seqs_hash) {
#	my $curr_id = $_;
#	
#	if(! $ids_present{$curr_id}) {
#		print ">$curr_id\n" . $seqs_hash{$curr_id} . "\n";
#	}
#}

# DONE