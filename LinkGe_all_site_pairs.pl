#! /usr/bin/perl

# PROGRAM: Runs LinkGe for all pairs of SNPs provided in a CLC-style SNP report (below). 
#	This script then tallies the following
#	for each pair of sites: num WT/WT reads; num WT/mut reads; num mut/WT reads; num
#	mut/mut reads; num OTHER reads.

# INPUT ARGUMENTS:
# snp_file: REQUIRED. ONE SNP report in CLC format containing the headers: 
#	(1) 'Reference Position', the site number; (2) 'Reference', the reference nucleotide; 
#	and (3) 'Allele', the variant nucleotide.
# BAM_file: REQUIRED. Standard format, aligned BAM file.
# max_dist: the maximum distance at which to search for pairs of sites in the same read.
#	DEFAULT: 1000.

# EXAMPLE CALLS:
# LinkGe_all_site_pairs.pl --snp_file=<CLC_SNP_report>.txt --BAM_file=<aligned_BAM>.bam
# LinkGe_all_site_pairs.pl --snp_file=high_freq_SNPs.txt --BAM_file=sample.bam --max_dist=1000

# Copyright (C) 2018 Chase W. Nelson

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

# DATE CREATED: March 30, 2018
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com
# AFFILIATION1: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA
# AFFILIATION2: Special Volunteer, Division of Cancer Epidemiology & Genetics, National
#     Cancer Institute, National Institutes of Health, Rockville, MD 20850, USA
# AFFILIATION3: BigPlant Consortium, Center for Genomics and Systems Biology, New York 
#     University, New York, NY 10003, USA

use strict;
#use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(max);

STDOUT->autoflush(1);

#########################################################################################
# Get the time and begin
my $time1 = time;
my $local_time1 = localtime;

print "\n##########################################################################################";
print "\n## LinkGe analysis initiated at local time $local_time1";
print "\n##########################################################################################\n";

#########################################################################################
# INITIALIZE VARIABLES (all optional and/or defaulted)
#my $genome_length; # can figure out
my $snp_file; # containing the variant sites of interest (and any others, to be ignored)
my $BAM_file;
my $max_dist = 1000;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "snp_file=s" => \$snp_file,
			"BAM_file=s" => \$BAM_file,
			 "max_dist=i" => \$max_dist)
			
			or die "\n### WARNING: Error in command line arguments. SNPGenie for OID terminated.\n\n";

unless(-f "$snp_file") {
	die "\n### WARNING: An existing --snp_file option must be provided\n".
		"### Script terminated.\n\n";
}

unless(-f "$BAM_file") {
	die "\n### WARNING: An existing --BAM_file option must be provided\n".
		"### Script terminated.\n\n";
}

if($max_dist < 2) {
	warn "\n### WARNING: max_dist must an integer (number of sites) greater than 1.\n" . 
		"### RE-SETTING max_dist=1000 and proceeding.\n\n";
}

# Get header names
my @header_names_arr = &get_header_names($snp_file,$snp_file);
#print "@header_names_arr";

my %header_indices;

# Determine the index of each column
for (my $i=0; $i<scalar(@header_names_arr); $i++) {
	my $curr_header = $header_names_arr[$i];
	$header_indices{$curr_header} = $i;
}

# Make sure we got all the expected headers
#my @required_headers = ('Reference Position', 'Reference', 'Allele', 'Count', 'Coverage', 
#						'Frequency');
my @required_headers = ('Reference Position', 'Reference', 'Allele');
						
foreach(@required_headers) {
	unless(exists $header_indices{$_}) {
		die "### DIE: the header name $_ is not present in the epitope database file.\n\n";
	}
}

# Test
#print "\n";
#foreach my $keys (keys %header_indices) {
#	print "key: $keys\nvalue: " . $header_indices{$keys} . "\n";
#}

##########################################################################################
# Store SNP data from provided $snp_file
my %site_nt_data;

open(SNP_FILE, "$snp_file") or die "Could not open file $snp_file\n";

my $snp_file_line = 0;

SNP: while(<SNP_FILE>) {
	chomp;
	
	if($snp_file_line == 0) { # skip the header
		$snp_file_line++;
#		next EPITOPE;
	} else { # store the record
		
		my @line_arr = split(/\t/, $_, -1);
		
#		print "\nline_arr is:\n";
#		foreach (@line_arr) {
#			print "$_\n";
#		}
		
		# STORE DATA HERE
		my $ref_position;
		if($line_arr[$header_indices{'Reference Position'}] =~ /(\d+)/) {
			$ref_position = $1;
		} else {
			warn "\n### WARNING: missing 'Reference Position'! Skipping.\n";
			next SNP;
		}
		
		my $WT = '';
		if($line_arr[$header_indices{'Reference'}] =~ /(\w+)/) {
			$WT = $1;
		} else {
			warn "\n### WARNING: missing 'Reference'! Skipping.\n";
			next SNP;
		}
		
		my $Mut = '';
		if($line_arr[$header_indices{'Allele'}] =~ /(\w+)/) {
			$Mut = $1;
		} else {
			warn "\n### WARNING: missing 'Allele'! Skipping.\n";
			next SNP;
		}
		
#		my $count;
#		if($line_arr[$header_indices{'Count'}] =~ /(\d+)/) {
#			$count = $1; 
#		} else {
#			warn "\n### WARNING: missing 'Count'! Skipping.\n";
#			next SNP;
#		}
#		
#		my $coverage;
#		if($line_arr[$header_indices{'Coverage'}] =~ /(\d+)/) {
#			$coverage = $1; 
#		} else {
#			warn "\n### WARNING: missing 'Coverage'! Skipping.\n";
#			next SNP;
#		}
#		
#		my $frequency;
#		if($line_arr[$header_indices{'Frequency'}] =~ /([\d\.]+)/) {
#			$frequency = $1; 
#		} else {
#			warn "\n### WARNING: missing 'Frequency'! Skipping.\n";
#			next SNP;
#		}
		
		# Store the data: in future, if slow, determine only sites we NEED to store first
		$site_nt_data{$ref_position}->{WT} = $WT;
		$site_nt_data{$ref_position}->{Mut} = $Mut;
#		$site_nt_data{$ref_position}->{count} = $count;
#		$site_nt_data{$ref_position}->{coverage} = $coverage;
#		$site_nt_data{$ref_position}->{frequency} = $frequency;
		
#		print "\n### Finished storing position $ref_position: WT=$WT Mut=$Mut count=$count coverage=$coverage frequency=$frequency\n";
		
	}
	
}

close SNP_FILE;

##########################################################################################
# Determine all unique site pairs falling within a certain distance of one another
my @all_sites_ordered = sort {$a <=> $b} keys %site_nt_data;
#print "\n@all_sites_ordered\n";

my %all_site_pairs;

SITE_QUERY: for(my $i = 0; $i < scalar(@all_sites_ordered); $i++) {
	my $curr_site = $all_sites_ordered[$i];
	
	for(my $j = $i + 1; $j < scalar(@all_sites_ordered); $j++) {
		my $curr_next_site = $all_sites_ordered[$j];
		
		#print "\ncurr_site=$curr_site\nnext_site=$curr_next_site";
		#print "curr_site + max_dist is " . ($curr_site + $max_dist) . "\n";
		
		if($curr_next_site <= ($curr_site + $max_dist)) { # it's a viable pair
			#print "site $curr_next_site is within range\n";
			push(@{$all_site_pairs{$curr_site}}, $curr_next_site);
		} else { # no more sites in range of $curr_site
			next SITE_QUERY;
		}
	}
}

##########################################################################################
# Run LinkGe for all unique site pairs

# Create new directory to store output files; follow LinkGe format
##my %file_to_call;

# collect time information to provide a unique filename to each output text file
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $theTime = $hour . $minute . $second . "_" . $year . "-" . $month . "-" . $dayOfMonth;

# Remove extension from SNP file to name the output file
my $snp_file_prefix = "SNPs";
if ($snp_file =~ /(\w+).txt/) {
	$snp_file_prefix = $1;
}
my $prefix = $snp_file_prefix . "_". $theTime;

my $output_file_name = "LinkGe\_$prefix\.txt";
my $temp_directory_name = "$prefix\_temp";

mkdir($temp_directory_name);
chdir($temp_directory_name);

print "\n#### CALLING: ####\n";

my $counter = 1;

foreach my $site_1 (sort {$a <=> $b} keys %all_site_pairs) {
	
	#print "\nsite $site_1 is within range of: ";
	
	foreach my $site_2 ((sort {$a <=> $b} @{$all_site_pairs{$site_1}})) {
		#print "$site_2 ";
		
		# RUN LinkGe
		my $LinkGe_path = `which LinkGe.pl`;
		chomp($LinkGe_path);
		
		my $this_call = "$LinkGe_path -l $site_1\,$site_2 ../$BAM_file -o $site_1\_$site_2";
		
		print "$this_call\n";
		`$this_call`;
		
##		$file_to_call{"$site_1\_$site_2\.txt"} = $this_call;
		
#		print "\nHere's what I'll call:\nLinkGe.pl -l $site_1\,$site_2 $BAM_file -o $site_1\_$site_2\n\n";
#		`LinkGe.pl -l $site_1\,$site_2 $BAM_file -o $site_1\_$site_2`;
		
		$counter++;
		
		if($counter == 16) { # don't wait to overwhelm the computer
#?			sleep(1); # turns out not necessary
			$counter = 1;
		}
		
	}
	
	#print "\n";
}


##########################################################################################
# Retrieve all text files in temp working directory
my @working_directory_files = glob "*.txt";

# Store data from files (ALL LinkGe output here), which have sites in their names
my @LinkGe_files;
my %read_linkage_data_hh;
my %sum_linkage_data;

#print "\nworking_directory_files: @working_directory_files\n";

foreach(@working_directory_files) {
	if($_ =~ /(\d+)_(\d+).txt/) {
		my $LinkGe_file = $_;
		
		my $this_site_1 = $1;
		my $this_site_2 = $2;
		
		# Order them numerically
		if($this_site_2 < $this_site_1) {
			my $temp_site = $this_site_1;
			$this_site_1 = $this_site_2;
			$this_site_2 = $temp_site;
		}
		
		#print "\nthis_site_1 = $this_site_1\n";
		#print "this_site_2 = $this_site_2\n";
		
		my $this_site_1_WT = $site_nt_data{$this_site_1}->{WT};
		my $this_site_1_Mut = $site_nt_data{$this_site_1}->{Mut};
		my $this_site_2_WT = $site_nt_data{$this_site_2}->{WT};
		my $this_site_2_Mut = $site_nt_data{$this_site_2}->{Mut};
		
		# Open file and read data
		open(LINKGE_FILE, "$LinkGe_file") or die "Could not open file $LinkGe_file\n";
		
		LINKGE: while(<LINKGE_FILE>) {
			chomp;
			
			if($_ =~ /([ACGTU-])(\d+)\:([ACGTU-])(\d+)\:\t(\d+)\t([\d\.]+)/) { # for example, C711:T933:	65	0.4676259
				#print "\ndata: $_";
				my $record_site_1_nt = $1;
				my $record_site_1_position = $2;
				#print "\nrecord_site_1_position = $record_site_1_position\n";
				
				my $record_site_2_nt = $3;
				my $record_site_2_position = $4;
				#print "\nrecord_site_2_position = $record_site_2_position\n";
				
				# Order them numerically
				if($record_site_2_position < $record_site_1_position) {
					my $temp_site = $record_site_1_position;
					$record_site_1_position = $record_site_2_position;
					$record_site_2_position = $temp_site;
				}
				
				#print "record_site_1_position = $record_site_1_position\n";
				#print "record_site_2_position = $record_site_2_position\n";
				
				my $record_reads = $5;
				my $record_freq = $6;
				
				#print "\nrecord_reads = $record_reads\n";
				
				# STORAGE FORMAT:
				# $read_linkage_data_hh{site1}->{site2}->{WT_WT}; # or WT_Mut/WT_unk/Mut_WT/Mut_Mut/Mut_unk/unk_WT/unk_Mut/unk_unk
				# NINE POSSIBLE COMBOS
				if($record_site_1_nt eq $this_site_1_WT) { # site 1 is the reference (WT)
					if($record_site_2_nt eq $this_site_2_WT) { # WT/WT
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{WT_WT} += $record_reads;
						$sum_linkage_data{WT_WT} += $record_reads;
					} elsif($record_site_2_nt eq $this_site_2_Mut) { # WT/Mut
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{WT_Mut} += $record_reads;
						$sum_linkage_data{WT_Mut} += $record_reads;
					} else { # WT/?
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{WT_unk} += $record_reads;
						$sum_linkage_data{WT_unk} += $record_reads;
					}
				} elsif($record_site_1_nt eq $this_site_1_Mut) { # site 1 is the allele (Mut)
					if($record_site_2_nt eq $this_site_2_WT) { # Mut/WT
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{Mut_WT} += $record_reads;
						$sum_linkage_data{Mut_WT} += $record_reads;
					} elsif($record_site_2_nt eq $this_site_2_Mut) { # Mut/Mut
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{Mut_Mut} += $record_reads;
						$sum_linkage_data{Mut_Mut} += $record_reads;
					} else { # Mut/?
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{Mut_unk} += $record_reads;
						$sum_linkage_data{Mut_unk} += $record_reads;
					}
				} else { # Mut/?
					if($record_site_2_nt eq $this_site_2_WT) { # ?/WT
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{unk_WT} += $record_reads;
						$sum_linkage_data{unk_WT} += $record_reads;
					} elsif($record_site_2_nt eq $this_site_2_Mut) { # ?/Mut
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{unk_Mut} += $record_reads;
						$sum_linkage_data{unk_Mut} += $record_reads;
					} else { # ?/?
						$read_linkage_data_hh{$record_site_1_position}->{$record_site_2_position}->{unk_unk} += $record_reads;
						$sum_linkage_data{unk_unk} += $record_reads;
					}
				}
			}
		}
		
		close LINKGE_FILE;
	}
}

##########################################################################################
# Exit temp directory; print individual totals to output file; delete temp directory
@working_directory_files = glob "*.txt"; # don't know why necessary
@working_directory_files = sort {$a <=> $b} @working_directory_files;

chdir("..");

open(OUTPUT_FILE, ">>$output_file_name");

foreach my $curr_LinkGe_file (@working_directory_files) {
#	print "\ncurr_LinkGe_file=$curr_LinkGe_file\n";
	if($curr_LinkGe_file =~ /(\d+)_(\d+).txt/) {
		my $this_site_1 = $1;
		my $this_site_2 = $2;
		
		my $curr_LinkGe_file_path = "$temp_directory_name\/$curr_LinkGe_file";
		
	#	print "\ncurr_LinkGe_file_path=$curr_LinkGe_file_path\n";
		print OUTPUT_FILE "#############################\n" . 
			"### TESTING SITES $this_site_1 and $this_site_2\n";
		
		# Open file and read data
		open(CURR_LINKGE_FILE, $curr_LinkGe_file_path) or die "Could not open LinkGe file $curr_LinkGe_file_path\n";
		
		while(<CURR_LINKGE_FILE>) {
			print OUTPUT_FILE "$_";
		}
		
		print OUTPUT_FILE "\n";
		
		close CURR_LINKGE_FILE;
		
		unlink $curr_LinkGe_file_path;
	}
}

close OUTPUT_FILE;
rmdir($temp_directory_name);

##########################################################################################
# Print result totals
print "\n########################\n##### READ TOTALS ######\n########################\n" .
	"WT_WT: $sum_linkage_data{WT_WT}\n" .
	"WT_Mut: $sum_linkage_data{WT_Mut}\n" .
	"WT_unk: $sum_linkage_data{WT_unk}\n" .
	"Mut_WT: $sum_linkage_data{Mut_WT}\n" .
	"Mut_Mut: $sum_linkage_data{Mut_Mut}\n" .
	"Mut_unk: $sum_linkage_data{Mut_unk}\n" .
	"unk_WT: $sum_linkage_data{unk_WT}\n" .
	"unk_Mut: $sum_linkage_data{unk_Mut}\n" .
	"unk_unk: $sum_linkage_data{unk_unk}\n" . 
	"########################\n";

my $total_defined_reads = $sum_linkage_data{WT_WT} + $sum_linkage_data{WT_Mut} + 
							$sum_linkage_data{Mut_WT} + $sum_linkage_data{Mut_Mut};
							
my $WT_WT_prop = $sum_linkage_data{WT_WT} / $total_defined_reads;
my $WT_Mut_prop = $sum_linkage_data{WT_Mut} / $total_defined_reads;
my $Mut_WT_prop = $sum_linkage_data{Mut_WT} / $total_defined_reads;
my $Mut_Mut_prop = $sum_linkage_data{Mut_Mut} / $total_defined_reads;

print "\n########################\n### WT & Mut SUMMARY ###\n########################\n" .
	"READS: $total_defined_reads\n" .
	"WT_WT READS: $WT_WT_prop\n" .
	"WT_Mut READS: $WT_Mut_prop\n" .
	"Mut_WT READS: $Mut_WT_prop\n" .
	"Mut_Mut READS: $Mut_Mut_prop\n" .
	"########################\n";

print "\n";

# Print a completion message to screen
&end_the_program;

#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################

#########################################################################################
sub get_header_names {
	# Originally assumed that we've received a tempfile ending in "_snpg9temp.txt"
	# However, we're now calling it at least once before creating the tempfile to 
	# see what kind of processing (e.g., Geneious to CLC) is needed prior to tempfile
	# creation. Must include capability to get headers for .CSV file
	my ($curr_snp_report_filename,$filename) = @_;
	#print "\n$curr_snp_report_filename\n";
	
	#my $newline_char = &detect_newline_char($curr_snp_report_filename);
	#my $old_newline = $/;
	#$/ = $newline_char;
	
	my $seen_tab_delimited = 0;
	my $seen_comma_delimited = 0;
	my $seen_vcf_tab_delimited = 0;
	my @line_arr;
	
	my $line = 0;
	open (CURRINFILE, $curr_snp_report_filename);
	#seek(CURRINFILE,0,0);
	while (<CURRINFILE>) {
		#print "$_";
		if($line == 0) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			if($_ =~/\t\w+\t/) { # it's TAB-delimited
				@line_arr = split(/\t/,$_);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,\w+,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_);
				#print "COMMA!!!!!";
				last;
			}

			$line++;
		} elsif($line > 0 && $_ =~ /^##/) {
			$line++;
		} elsif($line > 0 && ($_ =~ /^#CHROM/)) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			if($_ =~/\t/) { # it's TAB-delimited
				@line_arr = split(/\t/,$_);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_);
				#print "COMMA!!!!!";
				last;
			} else {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | CODON | WARNING
				
				# No change OR error should occur if the file does not, in fact, end
				# with this SUFFIX
				my $file_nm = $curr_snp_report_filename;
				#$file_nm =~ s/_snpg9temp.txt/.txt/;
				$file_nm =~ s/_\w\w\w\w.txt/.txt/;
				
				print ERROR_FILE "$filename\tN/A\tN/A\t".
					"File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
				close ERROR_FILE;
				chdir('..');
				
				#unlink $curr_snp_report_filename;
				
				die "\n\n## WARNING: The SNP Report $filename is ".
					"not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie ".
					"terminated\n\n";
			}
		} else {
			chdir('SNPGenie_Results');
			open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
			# FILE | PRODUCT | SITE | CODON | WARNING
			
			# No change OR error should occur if the file does not, in fact, end
			# with this SUFFIX
			my $file_nm = $curr_snp_report_filename;
			#$file_nm =~ s/_snpg9temp.txt/.txt/;
			$file_nm =~ s/_\w\w\w\w.txt/.txt/;
			
			print ERROR_FILE "$filename\tN/A\tN/A\t".
				"File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
			close ERROR_FILE;
			chdir('..');
			
			#unlink $curr_snp_report_filename;
			
			die "\n\n## WARNING: The SNP Report $filename is ".
				"not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie ".
				"terminated\n\n";
		}
	}
	seek(CURRINFILE,0,0);
	close CURRINFILE;
	#$/ = $old_newline;
	return @line_arr;
}

#########################################################################################
# End the program by notifying the screen at command line
sub end_the_program {
	my $time2 = time;
	my $local_time2 = localtime;
	
	my $time_diff = ($time2 - $time1);
	my $time_diff_rounded = sprintf("%.2f",$time_diff);
	my $mins_elapsed = ($time_diff / 60);
	my $whole_mins_elapsed = int($mins_elapsed);
	my $whole_mins_in_secs = ($whole_mins_elapsed * 60);
	my $secs_remaining = ($time_diff - $whole_mins_in_secs);
	my $secs_remaining_rounded = sprintf("%.2f",$secs_remaining);
	
	print "LinkGe analysis completed at local time $local_time2. The process took $time_diff_rounded secs, i.e., ".
			"$whole_mins_elapsed mins and $secs_remaining_rounded secs\n";

	print "\n################################################################################".
		"\n##                   LinkGe analysis completed successfully.                  ##".
		"\n## Please find individual site pair results in $output_file_name\n".
		"################################################################################".
		"\n\n\n"; 
}


