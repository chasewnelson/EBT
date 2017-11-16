#! /usr/bin/perl

# Takes in as arguments:
#	[0] a TAB-DELIMITED FILE;
#	[1] the column number to examine.
# OUTPUTS: the unique values in the column, i.e., duplicate values eliminated. Results 
#   WILL INCLUDE the header! ALSO places a summary statistics file in the working
#   directory, by the name *_unique_values_summary.txt.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# get_unique_values.pl <my_table.txt> <column_number>
#########################################################################################

# Copyright (C) 2017 Chase W. Nelson
# Date created: October 4, 2013
# AUTHOR: Chase W. Nelson
# CONTACT: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a National Science Foundation
# Graduate Research Fellowship (DGE-0929297), a National Science Foundation East Asian 
# and Pacific Summer Institutes Fellowship, and a University of South Carolina 
# Presidential Fellowship.

use strict;
#use warnings;

my $infile = $ARGV[0]; # TAB-DELIM FILE
my $column = $ARGV[1]; # COLUMN NUMBER

my %h_seen_keys;
my $col_header;

#open(INFILE,"$infile");
#my $line = 0;
#while(<INFILE>) {
#	if($line == 0) {
#		chomp;
#		my @a_curr_line = split("\t",$_);
#		$col_header = $a_curr_line[($column-1)];
#		$line++;
#	} else {
#		chomp;
#		my @a_curr_line = split("\t",$_);
#		$h_seen_keys{$a_curr_line[($column-1)]} += 1;
#	}
#}

my $line = 0;

open(INFILE,"$infile");
while(<INFILE>) {
	chomp;
	my @a_curr_line = split("\t",$_);
	$h_seen_keys{$a_curr_line[($column-1)]} += 1;
	$line++;
}

close INFILE;


my $max_hits = 1;

print "unique_values_including_header:\n";

#print "$col_header\n";

foreach (keys %h_seen_keys) {
	print "$_\n";
	
	if($h_seen_keys{$_} > $max_hits) {
		$max_hits = $h_seen_keys{$_};
	}
}

my %numHits2count;

# Print summary statistics files
# Generate new prefix
my $new_file_prefix;
if($infile =~/\.txt/) {
	$new_file_prefix = $`;
} else {
	$new_file_prefix = "infile";
}

open(OUT_UNIQUE_SUMMARY,">>$new_file_prefix\_unique_values_summary.txt");

for(my $i=1; $i<=$max_hits; $i++) {
	my @hits;	
		
	foreach (keys %h_seen_keys) {
		
		if($h_seen_keys{$_} == $i) {
			#print "$_\n";
			push(@hits, $_);
			$numHits2count{$i}++;
		}
	}
	
	if($numHits2count{$i} > 0) {
		print OUT_UNIQUE_SUMMARY "VALUES WITH $i HIT(S): ";
		
		foreach (@hits) {
			print OUT_UNIQUE_SUMMARY "$_ ";
		}
		print OUT_UNIQUE_SUMMARY "\n";
	}
}

print OUT_UNIQUE_SUMMARY "\n\nRESULTS:\nNUM RECORDS (incl. possible header): $line\n".
	"NUM UNIQUE VALUES: " . scalar(keys %h_seen_keys) .
	"\nNUM HITS PER VALUE:\n";

my $num_multi_hits = 0;
for(my $i=1; $i<=$max_hits; $i++) {
	if(! exists $numHits2count{$i}) {
		$numHits2count{$i} = 0;
	}
	
	print OUT_UNIQUE_SUMMARY "$i HIT(S): $numHits2count{$i}\n";

	if($i>1) {
		$num_multi_hits += $numHits2count{$i};
	}
	
}

print OUT_UNIQUE_SUMMARY "\nNUM VALUES WITH ≥2 HITS: $num_multi_hits\n\n";

close OUT_UNIQUE_SUMMARY;

#print "\n### ALL DONE ###\n\n";
