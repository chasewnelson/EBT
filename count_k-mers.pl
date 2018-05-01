#! /usr/bin/perl

# INPUT: Takes as arguments (1) ONE FASTA file with multiple aligned NUCLEOTIDE seqs;
#	(2) a GTF file with CDS annotations; and (3) max k (max k-mer length).
# OUTPUTS: command-line.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# count_k-mers.pl <aligned.fasta> <gene_annotations.gtf> <max_k>
# count_k-mers.pl HPV16_A1.fa HPV_annotations_plusE8.gtf 6
#########################################################################################

# Copyright (C) 2018 Chase W. Nelson
# DATE CREATED: May 1, 2018
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural 
#	History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# 	the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
#use warnings;
use List::Util qw(max sum);

my $fasta_filename = $ARGV[0];
my $gtf_filename = $ARGV[1];
my $max_k = 2;

#########################################################################################
# Read in the sequences from the file
my $seq = '';
my @seqs_arr;
#my $header = '';
#my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_filename") or die "Could not open FASTA file $fasta_filename\n";

print "\n\nReading in sequences from FASTA...\n\n";

while(<IN_FASTA>) {
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

close IN_FASTA;

$seq_num = scalar(@seqs_arr);


#########################################################################################
print "\n\nReading in coding annotations from GTF...\n\n";

# All products should be in GTF file, so don't need to get them from SNP reports as before
my @product_names_arr = &get_product_names_from_gtf($gtf_filename);
#@product_names_arr = sort {$a <=> $b} @product_names_arr; # take given order

# Determine all product coordinates from the start
my %product_coordinates_harr;
foreach my $product_name (@product_names_arr) {
	#print "\nproduct name is: $product_name\n";
	my @product_coord_arr = &get_product_coordinates($product_name);
	
	# Save to a hash of arrays
	push(@{$product_coordinates_harr{$product_name}->{product_coord_arr}},@product_coord_arr);
	#print "\n$product_name product_coord arr in harr: @{$product_coordinates_harr{$product_name}->{product_coord_arr}} \n\n";
}

# LABEL AS CODING IF WITHIN A CODING SEGMENT (from GTF)
my %hh_product_position_info;

# New segments approach
foreach (@product_names_arr) {
	my $product_name = $_;
	
	if(! exists $hh_product_position_info{$product_name}) {
		#print "product name is: $product_name\n";
		
		# Retrieve coordinates
		my @product_coord_arr = @{$product_coordinates_harr{$product_name}->{product_coord_arr}};
		
		my %product_starts;
		my %product_stops;
		
		my $num_segments = (@product_coord_arr / 2);
		for(my $i=1; $i<=scalar(@product_coord_arr); $i++) {
			$product_starts{$i} = $product_coord_arr[2*$i-2];
			$product_stops{$i} = $product_coord_arr[2*$i-1];
		}
		
		# Store start and stop for all segments
		foreach(sort {$a <=> $b} (keys %product_starts)) {
			my $this_start_key = 'start_' . $_;
			my $this_stop_key = 'stop_' . $_;
			$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
			$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
		}
		# Store number of segments
		$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
	}
}

#my @hh_keys = keys %hh_product_position_info;
#print "\n\n @hh_keys \n\n";

my @curr_products = sort(keys %hh_product_position_info);
#print "\n@curr_products";

# Ordering sorting products by the first segment's start site in the genome
my @curr_products_ordered_by_start = 
	sort { $hh_product_position_info{$a}->{start_1} <=> $hh_product_position_info{$b}->{start_1} } keys %hh_product_position_info;






##########################################################################################
# BEGIN K-MER ANALYSIS
my %k_motif_count_leading;
my %k_totals_leading = 0;

my %product_k_motif_count_leading;
my %product_k_totals_leading;

my %k_motif_count_lagging;
my %k_totals_lagging = 0;

my %product_k_motif_count_lagging;
my %product_k_totals_lagging;

for (my $k = 1; $k <= $max_k; $k++) {
	
	print "\n\nTabulating $k\-mers leading strand...\n\n";
	
	for(my $seq_id = 1; $seq_id <= scalar(@seqs_arr); $seq_id++) { # for each sequence
		
		my $curr_seq_index = $seq_id - 1;
		my $curr_seq = $seqs_arr[$curr_seq_index];
		#print "\ncurr_seq: $curr_seq\n";
		
		#################
		# LEADING STRAND
		for(my $site_id = 1; $site_id <= $last_seq_length; $site_id++) { # for each site
			if(($site_id + $k) <= $last_seq_length) { # we've got enough sequence left to record a k-mer of this size
				my $curr_site_index = $site_id - 1;
				
				my $curr_k_mer = substr($curr_seq, $curr_site_index, $k);
				#print "\ncurr_k_mer: $curr_k_mer";
				
				if(! ($curr_k_mer =~ 'N' || $curr_k_mer =~ '-')) {
					$k_motif_count_leading{$k}->{$curr_k_mer}++;
					$k_totals_leading{$k}++;
				}
			}
		}
		
		#################
		# LAGGING STRAND
		my $rev_seq = reverse($curr_seq);
		my $rev_compl = $rev_seq;
		$rev_compl =~ tr/ACGT/TGCA/;
		
		for(my $site_id = 1; $site_id <= $last_seq_length; $site_id++) { # for each site
			if(($site_id + $k) <= $last_seq_length) { # we've got enough sequence left to record a k-mer of this size
				my $curr_site_index = $site_id - 1;
				
				my $curr_k_mer = substr($rev_compl, $curr_site_index, $k);
				#print "\ncurr_k_mer: $curr_k_mer";
				
				if(! ($curr_k_mer =~ 'N' || $curr_k_mer =~ '-')) {
					$k_motif_count_lagging{$k}->{$curr_k_mer}++;
					$k_totals_lagging{$k}++;
				}
			}
		}
		
		#################
		# BY PRODUCT
		foreach my $curr_product (@curr_products_ordered_by_start) {
			#print "\nProduct: $curr_product, ";
			
			my $this_product_sequence = '';
			
			# SAVE a num_segments in the hh!
			my $num_segments = $hh_product_position_info{$curr_product}->{num_segments}; 
			#print "\nFor product $curr_product, the num_segments is: $num_segments\n";
			
			# Concatenate segments
			for(my $i = 1; $i <= $num_segments; $i++) {
				my $this_start_key = 'start_' . $i;
				my $this_stop_key = 'stop_' . $i;
				
				my $start = $hh_product_position_info{$curr_product}->{$this_start_key};
				my $stop = $hh_product_position_info{$curr_product}->{$this_stop_key};
				my $segment_length = ($stop - $start + 1);
				
				$this_product_sequence .= substr($curr_seq, $start - 1, $segment_length);
				
			}
			
			my $product_length = length($this_product_sequence);
			
			########################
			# PRODUCT LEADING STRAND
			for(my $site_id = 1; $site_id <= $product_length; $site_id++) { # for each site
				
				if(($site_id + $k) <= $product_length) { # we've got enough sequence left to record a k-mer of this size
					my $curr_site_index = $site_id - 1;
					
					my $curr_k_mer = substr($this_product_sequence, $curr_site_index, $k);
					#print "\ncurr_k_mer: $curr_k_mer";
					
					if(! ($curr_k_mer =~ 'N' || $curr_k_mer =~ '-')) {
						$product_k_motif_count_leading{$curr_product}->{$k}->{$curr_k_mer}++;
						$product_k_totals_leading{$curr_product}->{$k}++;
					}
				}
				
			}
			
			########################
			# PRODUCT LAGGING STRAND
			my $rev_product_seq = reverse($this_product_sequence);
			my $rev_compl_product_seq = $rev_product_seq;
			$rev_compl_product_seq =~ tr/ACGT/TGCA/;
			
			for(my $site_id = 1; $site_id <= $product_length; $site_id++) { # for each site
				
				if(($site_id + $k) <= $product_length) { # we've got enough sequence left to record a k-mer of this size
					my $curr_site_index = $site_id - 1;
					
					my $curr_k_mer = substr($rev_compl_product_seq, $curr_site_index, $k);
					#print "\ncurr_k_mer: $curr_k_mer";
					
					if(! ($curr_k_mer =~ 'N' || $curr_k_mer =~ '-')) {
						$product_k_motif_count_lagging{$curr_product}->{$k}->{$curr_k_mer}++;
						$product_k_totals_lagging{$curr_product}->{$k}++;
					}
				}
				
			}
			
		}
		
	} # done with last sequence
	
}

# PRINT LEADING GENOME
foreach my $k (sort {$a <=> $b} keys %k_motif_count_leading) {
	my $this_k_total = $k_totals_leading{$k};
	
	foreach my $mer (sort keys %{$k_motif_count_leading{$k}}) {
		my $mer_count = $k_motif_count_leading{$k}->{$mer};
		my $mer_prop = '';
		
		if($this_k_total > 0) {
			$mer_prop = $mer_count / $this_k_total;
		}
		
		print "genome\t+\t$mer\t$mer_prop\t$mer_count\n";
	}
}

# PRINT LEADING PRODUCTS
foreach my $product (@curr_products_ordered_by_start) {	
	foreach my $k (sort {$a <=> $b} keys %{$product_k_motif_count_leading{$product}}) {
		my $this_k_total = $product_k_totals_leading{$product}->{$k};
		
		foreach my $mer (sort keys %{$product_k_motif_count_leading{$product}->{$k}}) {
			my $mer_count = $product_k_motif_count_leading{$product}->{$k}->{$mer};
			my $mer_prop = '';
			
			if($this_k_total > 0) {
				$mer_prop = $mer_count / $this_k_total;
			}
			
			print "$product\t+\t$mer\t$mer_prop\t$mer_count\n";
		}
	}
}

# PRINT LAGGING GENOME
foreach my $k (sort {$a <=> $b} keys %k_motif_count_lagging) {
	my $this_k_total = $k_totals_lagging{$k};
	
	foreach my $mer (sort keys %{$k_motif_count_lagging{$k}}) {
		my $mer_count = $k_motif_count_lagging{$k}->{$mer};
		my $mer_prop = '';
		
		if($this_k_total > 0) {
			$mer_prop = $mer_count / $this_k_total;
		}
		
		print "genome\t-\t$mer\t$mer_prop\t$mer_count\n";
	}
}

# PRINT LAGGING PRODUCTS
foreach my $product (@curr_products_ordered_by_start) {	
	foreach my $k (sort {$a <=> $b} keys %{$product_k_motif_count_lagging{$product}}) {
		my $this_k_total = $product_k_totals_lagging{$product}->{$k};
		
		foreach my $mer (sort keys %{$product_k_motif_count_lagging{$product}->{$k}}) {
			my $mer_count = $product_k_motif_count_lagging{$product}->{$k}->{$mer};
			my $mer_prop = '';
			
			if($this_k_total > 0) {
				$mer_prop = $mer_count / $this_k_total;
			}
			
			print "$product\t-\t$mer\t$mer_prop\t$mer_count\n";
		}
	}
}

print "\n\nCOMPLETED.\n\n";

exit;


#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################

#########################################################################################
sub get_product_names_from_gtf {
	my ($gtf_filename) = @_;
	#print "\n\n$gtf_filename\n\n";
	my %products_hash;
	open (CURRINFILE, $gtf_filename);
	while (<CURRINFILE>) {
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			#print "this_line: $_";
			if($_ =~/\s*gene_id\s*\"gene\:([\w\s\.\-\:']+)\"/) { # transcript_id not a problem
				$products_hash{$1} = 1;
#				$seen_sense_strand_products = 1;
			} elsif($_ =~ /\s*gene_id\s*\"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
#				$seen_sense_strand_products = 1;
			} elsif($_ =~/\s*gene_id\s*\"([\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
#				$seen_sense_strand_products = 1;
			} else {
				#unlink $curr_snp_report_filename;
				
				die "\n\n## WARNING: CDS annotation(s) in $gtf_filename does not have a ".
					"gene_id. SNPGenie terminated.\n\n";
			}
		}
	}
	close CURRINFILE;
	
	my @product_names = keys %products_hash;
	#print "\n@product_names\n\n";
	return @product_names;
}

#########################################################################################
# Get an array with the product's start and stop sites, found in the .gtf file.
# Cause the program to DIE if the coordinates do not yield a multiple of three (i.e., 
# there is an incomplete codon)
# Returns an array with:
#	returned[0] = starting site
#	returned[1] = stop site
#	IF MORE SEGMENTS:
#	returned[2] = starting site 2
#	returned[3] = stop site 2
#	etc.
sub get_product_coordinates {
	my $product = $_[0];
	
	#print "\nThis time called for product: $product\n";
	
	my %start_site_h; # {segment #}->{start coordinate}
	my %stop_site_h;
	my %segment_lengths_h;
	
	open (CURRINFILE, $gtf_filename);
	while (<CURRINFILE>) { # go through the GTF file
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand #COMEBACK to \d+]
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			my $this_line = $_;
			my $product_present = 0;
			my $this_start;
			my $this_stop;
			
			#print "\nPRODUCT: $product\n";
			
			# incomplete segments
			if ($this_line =~/\s*gene_id\s+"$product\"/) { 
				if ($this_line =~/CDS\t(\d+)\t(\d+)/) {
					$product_present = 1;
					$this_start = $1;
					$this_stop = $2;
					
					#print "\n## $product\: $this_start\-$this_stop\n";
				}
			} elsif($this_line =~/\s*gene_id\s+"gene\:$product\"/) {
				if ($this_line =~/CDS\t(\d+)\t(\d+)/) {
					$product_present = 1;
					$this_start = $1;
					$this_stop = $2;
					#print "\n## $product\: $this_start\-$this_stop\n";
				}
			}
			
			# ADD SOME ERROR IF THE GENE NAME WASN'T RECOGNIZED BUT OUGHT TO HAVE BEEN?
			
			if($product_present == 1) {
				my $curr_max_key = max(keys %start_site_h); # this is the number of segments so far
				#print "curr_max_key is: $curr_max_key\n";
				my $curr_key = ($curr_max_key + 1);
				#print "curr_key is: $curr_key\n";
				$start_site_h{$curr_key} = $this_start;
				$stop_site_h{$curr_key} = $this_stop;
			}
		}
	}
	close CURRINFILE;
	
	#new incomplete segments
	my $num_segments = max(keys %start_site_h); # this is the number of segments total
	#print "\nAfter while, $product num segments: $num_segments\n";
	my $sum_of_lengths;
	
	for(my $i=1; $i<=$num_segments; $i++) {
		$segment_lengths_h{$i} = ($stop_site_h{$i} - $start_site_h{$i} + 1);
		$sum_of_lengths += $segment_lengths_h{$i};
	}
	
	# Make sure the sum of all nucleotides for this coding product add to a multiple of 3
	# incomplete segments
	if (($sum_of_lengths % 3) != 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		# FILE | PRODUCT | SITE | CODON | WARNING
		print ERROR_FILE "N/A\t$product\tN/A\t".
			"Coordinates do not specify complete codon set (nucleotides a multiple of 3). SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: The CDS coordinates for gene $product in the gtf file ".
			"do not yield a set of complete codons,\n".
			"## or are absent from the file. The number of nucleotides must ".
			"be a multiple of 3.\n## SNPGenie terminated.\n\n";
	}
	
	my @coord_arr;
	
	# Sort the segments
	my @sorted_segment_numbers;
	
	@sorted_segment_numbers = sort {($start_site_h{$a} <=> $start_site_h{$b}) || ($a <=> $b) } keys %start_site_h;
	
	foreach(@sorted_segment_numbers) {
		push(@coord_arr,$start_site_h{$_});
		push(@coord_arr,$stop_site_h{$_});
	}
	
#	for(my $i=1; $i<=$num_segments; $i++) {
#		push(@coord_arr,$start_site_h{$i});
#		push(@coord_arr,$stop_site_h{$i});
#	}
	
	return @coord_arr;
}


