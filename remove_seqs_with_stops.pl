#! /usr/bin/perl

# You want to remove all sequences containing mid-sequence STOP codons.

# Takes in as arguments:
#	[0] one FASTA file containing multiple coding sequences that all begin at the first 
#	    position of the first codon (they need not be aligned to each other).
# OUTPUTS: a new file in the working directory without the STOP-containing sequences.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# remove_seqs_with_stops.pl <my_fasta_file.fasta>
#########################################################################################

# Copyright (C) 2017 Chase W. Nelson

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# DATE CREATED: November 14, 2017
# AUTHOR: Chase W. Nelson
# CONTACT: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use warnings;

my $filename = $ARGV[0];

# Read in the sequence from the file
my $seq = '';
my @seqs_arr;
my $header = '';
my @headers_arr;
my $seq_num = 0;

open(IN, "$filename") or die "Could not open file $filename\n";

print "\n### Reading in sequences...\n";

while(<IN>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$header = $_;
			$seq_num ++;
		} else {
			push(@seqs_arr,$seq);
			push(@headers_arr,$header);
			$header = $_;
			$seq_num ++;
			#print "\nseq: $seq\n";
			$seq = '';
		}
	} else {
		$seq .= $_;
	}
}

push(@seqs_arr,$seq);
push(@headers_arr,$header);

unless(@seqs_arr == @headers_arr) { die "\n### Num seqs does not equal num headers.\n\n"; }

close IN;

my %seqs_with_stops;

SEQUENCES: for(my $i=0; $i<scalar(@headers_arr); $i++) { # for each sequence
	my $curr_seq = $seqs_arr[$i];
	
	for(my $j=0; $j<=length($curr_seq); $j+=3) { # for each codon
		my $curr_codon = substr($curr_seq,$j,3);
		$curr_codon = uc($curr_codon); # uc returns uppercase
		$curr_codon =~ tr/U/T/;
		
		if($curr_codon eq 'TAA' || $curr_codon eq 'TAG' || $curr_codon eq 'TGA') { # if STOP
			if($j != (length($curr_seq) - 3)) { # if not the last codon in sequence
				$seqs_with_stops{$headers_arr[$i]} = 1;
				next SEQUENCES;
			}
		}
		
	}
	
}

# Loop through FASTA again and print only if no mid-sequence STOPS were found

# Output file name
my $output_filename;
if($filename =~ '.fasta') {
	$output_filename = $` . "_wo_stops.fasta";
} elsif($filename =~ '.fa') {
	$output_filename = $` . "_wo_stops.fa";
} else {
	$output_filename = "input_wo_stops.fasta";
}

open(OUT, ">>$output_filename");

open(IN, "$filename") or die "Could not open file $filename\n";

print "\n### Printing sequences to $output_filename\...\n";

my $print_seq_flag = 0;

while(<IN>) {
	chomp;
	if(/>/) {
		$header = $_;
		
		if(exists $seqs_with_stops{$header}) {
			$print_seq_flag = 0;
		} else {
			$print_seq_flag = 1;
			print OUT "$_\n";
		}
		
	} elsif($print_seq_flag == 1) {
		print OUT "$_\n";
	}
}

close IN;
close OUT;

my @excluded_headers = sort keys %seqs_with_stops;
my $num_excluded_seqs = scalar(@excluded_headers);

print "\n##########################################################################################\n";
print "### COMPLETED ###\n";
print "### The following $num_excluded_seqs sequences were excluded:\n";

foreach my $stop_headers (@excluded_headers) {
	print "$stop_headers\n";
}

print "##########################################################################################\n\n\n";
