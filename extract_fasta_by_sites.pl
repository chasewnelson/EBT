#! /usr/bin/perl

# Takes in as arguments:
#	[0] one FASTA file containing one or more aligned sequences; and
#	[1] a GTF file containing the CDS products to extract.
# OUTPUTS: one FASTA file for each CDS record in the GTF file. Just the DNA segment
#   corresponding to the CDS coordinates of each record will be present in the resulting
#   FASTA files.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# extract_fasta_by_sites.pl <multiple_aligned_seqs.fasta> <gene_coordinates_to_extract.gtf>
#########################################################################################

# Copyright (C) 2014 Chase W. Nelson
# Date created: 07 December 2014
# AUTHOR: Chase W. Nelson
# Affiliation1: Austin L. Hughes lab, University of South Carolina (Columbia, SC 29208, USA)
# Affiliation2: Wen-Hsiung Li lab, Academia Sinica (Taipei, Taiwan)
# Contact1: cwnelson88@gmail.com
# Contact2: nelsoncw@email.sc.edu

# ACKNOWLEDGMENTS: Written by C.W.N. with support from a National Science Foundation
# Graduate Research Fellowship (DGE-0929297), a National Science Foundation East Asian 
# and Pacific Summer Institutes Fellowship, and a University of South Carolina 
# Presidential Fellowship.

use strict;
use warnings;

# ONLY WORKS FOR + STRAND RIGHT NOW

# Check that an argument is given
if(! $ARGV[1]) {
	die "\n###WARNING: Two arguments must be supplied: (1) <multiple_aligned_seqs.fasta>; ".
		"(2) <gene_coordinates_to_extract.gtf>\n\n";
}

#my $fasta_file = glob "*.fasta";
#if(! $fasta_file) {
#	$fasta_file = glob "*.fa";
#}
my $fasta_file = $ARGV[0];
#my $gtf_file = glob "*.gtf"; #$ARGV[1];
my $gtf_file = $ARGV[1];

# Get product names first
my %products_hash;
open (GTF_FILE, $gtf_file);
while (<GTF_FILE>) {
	if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
		if($_ =~ /\s*gene_id\s*\"gene\:([\w\s\.\-\:']+)\"/) {
			$products_hash{$1} = 1;
		} elsif($_ =~ /\s*gene_id\s*\"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
			$products_hash{$1} = 1;
		} elsif($_ =~/\s*gene_id\s*\"([\w\s\.\-\:']+)\"/) {
			$products_hash{$1} = 1;
		} 
	}
}
close GTF_FILE;

my @product_names = keys %products_hash;

# Now get product coordinates
my %products_sites_hh;
foreach my $product (@product_names) {
	
	#print "\nproduct: $product\n";
	
	my $start_site_1;
	my $stop_site_1;
	
	my $start_site_2;
	my $stop_site_2;
	
	open (CURRINFILE, $gtf_file);
	while (<CURRINFILE>) {
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			chomp;
			
			if (! $start_site_1) {
				if ($_ =~/gene_id "$product";/) {
					if ($_ =~/CDS\t(\d+)\t(\d+)/) {
						$start_site_1 = $1;
						$stop_site_1 = $2;
					}
				} elsif($_ =~/gene_id "gene\:$product";/) {
					if ($_ =~/CDS\t(\d+)\t(\d+)/) {
						$start_site_1 = $1;
						$stop_site_1 = $2;
					}
				}
			} else {
				if ($_ =~/gene_id "$product";/) {
					if ($_ =~/CDS\t(\d+)\t(\d+)/) {
						$start_site_2 = $1;
						$stop_site_2 = $2;
						last; # This might be changed if we go on to add more segments
					}
				} elsif($_ =~/gene_id "gene\:$product";/) {
					if ($_ =~/CDS\t(\d+)\t(\d+)/) {
						$start_site_2 = $1;
						$stop_site_2 = $2;
						last; # This might be changed if we go on to add more segments
					}
				}
			}
		}
	}
	close CURRINFILE;
	
	if ((!$start_site_2) && (($stop_site_1 - $start_site_1 + 1) % 3) != 0) {
		die "\n\n## WARNING: The CDS coordinates for gene $product in the gtf file ".
			"do not yield a set of complete codons,\n".
			"## or are absent from the file. The number of nucleotides must ".
			"be a multiple of 3.\n## SNPGenie terminated.\n\n";
	}
	
	if ($start_site_2 && (((($stop_site_1 - $start_site_1 + 1) + ($stop_site_2 - $start_site_2 + 1)) % 3) != 0)) {
		die "\n\n## WARNING: The CDS coordinates for gene $product in the gtf file ".
			"do not yield a set of complete codons,\n".
			"## or are absent from the file. The number of nucleotides must ".
			"be a multiple of 3.\n## SNPGenie terminated.\n\n";
	}
	
	my @coord_arr;
	
	if (! $start_site_2) {
		$products_sites_hh{$product}->{start}=$start_site_1;
		$products_sites_hh{$product}->{stop}=$stop_site_1;
	} else {
		if ($start_site_1 < $start_site_2) {
			$products_sites_hh{$product}->{start}=$start_site_1;
			$products_sites_hh{$product}->{stop}=$stop_site_1;
			$products_sites_hh{$product}->{start2}=$start_site_2;
			$products_sites_hh{$product}->{stop2}=$stop_site_2;
		} else {
			$products_sites_hh{$product}->{start}=$start_site_2;
			$products_sites_hh{$product}->{stop}=$stop_site_2;
			$products_sites_hh{$product}->{start2}=$start_site_1;
			$products_sites_hh{$product}->{stop2}=$stop_site_1;
		}
	}
}

# Generate new prefix
my $new_file_prefix;
if($fasta_file =~/\.fasta/) { 
	$new_file_prefix = $`;
} elsif($fasta_file =~/\.txt/) {
	$new_file_prefix = $`;
} elsif($fasta_file =~/\.fa/) {
	$new_file_prefix = $`;
} else {
	$new_file_prefix = "new_file";
}


my $seq_counter = 0;
my $motif_counter = 0;

my $seq = '';
my $header = '';
my @seqs_arr;
my @headers_arr;
#my @id_arr;
my @seq_lengths;
my $last_seq_length;

open(FASTA_FILE, "$fasta_file") or die "Could not open file $fasta_file\n";

print "\n\nBuilding sequences...\n\n";

while(<FASTA_FILE>) {
	chomp;
	if(/>/) {
		if($seq_counter == 0) {
			$header = $_;
			$seq_counter ++;
		} else {
			push(@seqs_arr,$seq);
			push(@headers_arr,$header);
			$header = $_;
			$seq_counter ++;
			
			#print "\nThis sequence is:\n\n$seq\n\n";
			
			my $this_seq_length = length($seq);
			
			push(@seq_lengths,$this_seq_length);
			#print "\nseq $seq_num is of length $this_seq_length\n";
			
			if($last_seq_length && ($last_seq_length != $this_seq_length)) {
				die "\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
			} else {
				$last_seq_length = $this_seq_length;
				#print "\nseq: $seq\n";
				$seq = '';
				#print "\nGot here\n";
			}
		}
	} else {
		$seq .= $_;
	}
}

push(@seqs_arr,$seq);
push(@seq_lengths,$last_seq_length); # FIX THIS FOR ONE SEQ
push(@headers_arr,$header);

my $this_seq_length = length($seq);
			
push(@seq_lengths,$this_seq_length);
#print "\nseq $seq_num is of length $this_seq_length\n";

if($last_seq_length && ($last_seq_length != $this_seq_length)) {
	die "\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
} else {
	$last_seq_length = $this_seq_length;
	#print "\nseq: $seq\n";
	$seq = '';
}

close FASTA_FILE;


foreach my $product (@product_names) {
	
	my $this_file_name = "$new_file_prefix\_$product\.fasta";
	
	print "Printing $this_file_name\...\n\n";
	
	open(OUTPUT_FILE, ">>$this_file_name");
	
	if($products_sites_hh{$product}->{start2}) { # if there's a second segment
		print "\nGot here 2\n";
		my $start = $products_sites_hh{$product}->{start};
		my $start_index = ($start - 1);
		my $stop = $products_sites_hh{$product}->{stop};
		
		my $start2 = $products_sites_hh{$product}->{start2};
		my $start2_index = ($start2 - 1);
		my $stop2 = $products_sites_hh{$product}->{stop2};
		
		my $extracted_length1 = ($stop - $start + 1);
		my $extracted_length2 = ($stop2 - $start2 + 1);
		my $extracted_length_total = $extracted_length1 + $extracted_length2;
		
		if($start > $stop) {
			die "\nThe start site must occur before the stop site\n\n";
		} elsif($start2 > $stop2) {
			die "\nThe start site must occur before the stop site\n\n";
		}
		
		if($extracted_length_total % 3 != 0) {
			print "\n\nWARNING: extracted region is not a multiple of 3 (will not form whole \n".
				"codons if this is a coding region.\n\n";
		}
		
		# MUST DEFINE A PRODUCT-SPECIFIC FILE
		
		# FOREACH SEQ
		for(my $seq_index = 0; $seq_index < @seqs_arr; $seq_index++) {
			my $seq = $seqs_arr[$seq_index];
			my $extract_seq1 = substr($seq, $start_index, $extracted_length1); # substr( variable, start [0 index], length)
			my $extract_seq2 = substr($seq, $start2_index, $extracted_length2); # substr( variable, start [0 index], length)
			my $extract_seq_total = $extract_seq1 . $extract_seq2;
			
			#my $new_file_name = $new_file_prefix . "_$product\.fasta";
			#$this_file_name .= $new_file_prefix . $this_file_name;
			#open(OUTPUT_FILE, ">>$this_file_name");
			my $new_header = $headers_arr[$seq_index] . "product $product: $start\.\.$stop,$start2\.\.$stop2 ".
				"length $extracted_length_total";
			print OUTPUT_FILE "$new_header\n";
	
			for(my $i = 0; $i < length($extract_seq_total); $i+=60) {
				my $seq_line_to_print = substr($extract_seq_total, $i, 60);
				print OUTPUT_FILE "$seq_line_to_print\n";
			}
		}
		
	} else {
		#print "\nGot here 3\n";
		my $start = $products_sites_hh{$product}->{start};
		my $start_index = ($start - 1);
		my $stop = $products_sites_hh{$product}->{stop};
		
		my $extracted_length = ($stop - $start + 1);
		
		if($start > $stop) {
			die "\nThe start site must occur before the stop site\n\n";
		}
		
		if($extracted_length % 3 != 0) {
			print "\n\nWARNING: extracted region is not a multiple of 3 (will not form whole \n".
				"codons if this is a coding region.\n\n";
		}
		
		# MUST DEFINE A PROUCT-SPECIFIC FILE
		
		# FOREACH SEQ
		for(my $seq_index = 0; $seq_index < @seqs_arr; $seq_index++) {
#			print "my seq index: $seq_index\n";
			my $seq = $seqs_arr[$seq_index];
			my $extract_seq = substr($seq, $start_index, $extracted_length); # substr( variable, start [0 index], length)
		
			#my $new_file_name = $new_file_prefix . "_$product\.fasta";
			#$this_file_name .= $new_file_prefix . $this_file_name;
			#open(OUTPUT_FILE, ">>$this_file_name");
			my $new_header = $headers_arr[$seq_index] . "product $product: $start\.\.$stop length $extracted_length";
			print OUTPUT_FILE "$new_header\n";
	
			for(my $i = 0; $i < length($extract_seq); $i+=60) {
				my $seq_line_to_print = substr($extract_seq, $i, 60);
				print OUTPUT_FILE "$seq_line_to_print\n";
			}
		}
	}
	close OUTPUT_FILE;
}

print "\n### DONE ###\n\n";
