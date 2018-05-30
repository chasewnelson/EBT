#! /usr/bin/perl

# PROGRAM: determines all haplotypes in a given FASTA alignment for the sites provided.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# haplotypes_provided_sites.pl --sites_file=PBS_neg_sites.txt --fasta_file=A1.fa > PBShaplo_A1.txt
#########################################################################################

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

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

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
print "\n## Haplotype analysis initiated at local time $local_time1";
print "\n##########################################################################################\n";


#########################################################################################
# INITIALIZE VARIABLES
my $fasta_file; # containing MULTIPLE ALIGNED sequences, e.g., of a sublineage
my $sites_file;
#my $procs_per_node;
my $num_bootstraps = 10000; # default 10000

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "fasta_file=s" => \$fasta_file,
			"sites_file=s" => \$sites_file,
#			"procs_per_node=i" => \$procs_per_node,
			"num_bootstraps=i" => \$num_bootstraps)
			
			or die "\n### WARNING: Error in command line arguments. Script terminated.\n\n";

unless(-f "$fasta_file") {
	die "\n### WARNING: An existing --fasta_file option must be provided\n".
		"### Script terminated.\n\n";
}
unless(-f "$sites_file") {
	die "\n### WARNING: An existing --sites_file option must be provided\n".
		"### Script terminated.\n\n";
}

#unless($procs_per_node >= 1) {
#	die "\n### WARNING: A positive --procs_per_node must be provided for parallelism.\n".
#		"### Script terminated.\n\n";
#}


#########################################################################################
# DETERMINE OUTPUT FILE PREFIX
my $new_file_prefix;
if($sites_file =~ '.txt') {
	$new_file_prefix = $`;
} else {
	$new_file_prefix = "INPUT";
}


#########################################################################################
# STORE FASTA SEQUENCE & FULL-SEQUENCE HAPLOTYPES

my $seq = '';
my @seqs_arr;
#my $header = '';
#my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file") or die "Could not open FASTA file $fasta_file\n";

print "\nReading in sequences from FASTA...\n";

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

push(@seqs_arr, $seq);
#push(@headers_arr,$header);

close IN_FASTA;

$seq_num = scalar(@seqs_arr);
my $genome_length = length($seq);


#########################################################################################
# READ IN AND STORE SITES
my @sites;

open(IN_SITES, "$sites_file") or die "Could not open sites file $sites_file\n";

print "\n\nReading in sites from file...\n\n";

while(<IN_SITES>) {
	chomp;
	
	my @sites_line_arr = split(/\s+/, $_, -1);
	
	foreach (@sites_line_arr) {
		if($_ =~ /\d+/) {
			push(@sites, $_);
		}
	}
}

close IN_SITES;

@sites = sort {$a <=> $b} @sites;
my $num_sites = scalar(@sites);

print "\nANALYZING $num_sites SITES:\n@sites\n\n";

my $max_site = max(@sites);

if($max_site > $genome_length) {
	die "\n### WARNING: sites fall outside of genome length $genome_length\.\n".
		"### SCRIPT TERMINATED.\n\n";
}

foreach (@sites) {
	print "site_$_\t";
}

print "\n";


##########################################################################################
# Loop through sequences and record haplotypes for sites
my %haplotypes;
foreach my $seq (@seqs_arr) {
	my $haplotype;
	
	foreach my $site (@sites) {
		$haplotype .= substr($seq, $site - 1, 1);
	}
	
	$haplotypes{$haplotype}++;
	
}

## ONLY DEFINED HAPLOTYPES
# Determine the most common haplotype
my $maj_haplotype = '';
my $maj_haplotype_count = 0;

foreach my $haplotype (keys %haplotypes) {
	my $curr_haplotype_count = $haplotypes{$haplotype};
	
	unless($haplotype =~ /N/) { # only fully defined ones
		if($curr_haplotype_count > $maj_haplotype_count) {
			$maj_haplotype_count = $curr_haplotype_count;
			$maj_haplotype = $haplotype;
		}
	}
}

print "\nMAJOR HAPLOTYPE (FULLY-DEFINED ONLY): $maj_haplotype\n\n";

print "haplotype\tcount\tdiffs_from_maj_haplotype\n";

foreach my $haplotype (keys %haplotypes) {
	unless($haplotype =~ /N/) { # only fully defined ones
		# Count differences from major haplotype
		my $diffs_from_maj_haplotype = 0;
		
		for(my $i = 0; $i < length($haplotype); $i++) {
			if(substr($haplotype, $i, 1) ne substr($maj_haplotype, $i, 1)) {
				$diffs_from_maj_haplotype++;
			}
		}
		
		print "$haplotype\t$haplotypes{$haplotype}\t$diffs_from_maj_haplotype\n";
	}
}

## ALL HAPLOTYPES
# Determine the most common haplotype_undef
my $maj_haplotype_undef = '';
my $maj_haplotype_undef_count = 0;

foreach my $haplotype (keys %haplotypes) {
	my $curr_haplotype_undef_count = $haplotypes{$haplotype};
	
	if($curr_haplotype_undef_count > $maj_haplotype_undef_count) {
		$maj_haplotype_undef_count = $curr_haplotype_undef_count;
		$maj_haplotype_undef = $haplotype;
	}
}

print "\nMAJOR HAPLOTYPE (ALL): $maj_haplotype_undef\n\n";

print "haplotype\tcount\tdiffs_from_maj_haplotype\n";

foreach my $haplotype (keys %haplotypes) {
	# Count differences from major haplotype
	my $diffs_from_maj_haplotype_undef = 0;
	
	for(my $i = 0; $i < length($haplotype); $i++) {
		if(substr($haplotype, $i, 1) ne substr($maj_haplotype_undef, $i, 1)) {
			$diffs_from_maj_haplotype_undef++;
		}
	}
	
	print "$haplotype\t$haplotypes{$haplotype}\t$diffs_from_maj_haplotype_undef\n";

}

# Print a completion message to screen
&end_the_program;

exit; 

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
		"\n##                 Haplotype analysis completed successfully.                 ##".
		"\n##                           Please have fun now.                             ##".
		"\n################################################################################".
		"\n\n\n"; 
}


