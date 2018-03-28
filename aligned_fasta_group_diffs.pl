#! /usr/bin/perl

# INPUT: takes in all .fa* files in the working directory, which must be aligned to one 
#    another. Considers each fasta file a meaningful "group" of sequences.
# OUTPUTS: a file describing positions at which the groups differ in their major
#    nucleotide. The user may provide the following arguments to control how these sites
#    are determined: <to be added>

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# aligned_fasta_group_diffs.pl
# aligned_fasta_group_diffs.pl --min_variant_maj_nt_freq=.9 --min_site_coverage=8
#########################################################################################

# Copyright (C) 2017 Chase W. Nelson
# Date created: November 15, 2017
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.


use strict;
#use warnings;
use Getopt::Long;
use Data::Dumper;

STDOUT->autoflush(1); # requires Data::Dumper, but needed

#########################################################################################
# INITIALIZE (OPTIONAL) INPUT VARIABLES
my $min_variant_maj_nt_freq;
my $min_site_coverage;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "min_variant_maj_nt_freq:f" => \$min_variant_maj_nt_freq, # optional floating point parameter
			"min_site_coverage:i" => \$min_site_coverage) # optional integer parameter
			
			or die "\n### WARNING: Error in command line arguments. Script terminated.\n\n";
			# If an argument is called as a flag, its value is 0; if not called, it's null
			
if(! $min_variant_maj_nt_freq) { # null or 0
	$min_variant_maj_nt_freq = 0.5; # default behavior
} elsif($min_variant_maj_nt_freq < 0.5) {
	die "\n### WARNING: The --min_variant_maj_nt_freq option must ≥0.5\n".
		"### Script terminated.\n\n";
}

if(! $min_site_coverage) { # null or 0
	$min_site_coverage = 5; # default behavior
} elsif($min_site_coverage < 1) {
	die "\n### WARNING: The --min_site_coverage option must ≥1\n".
		"### Script terminated.\n\n";
}

#########################################################################################
# PREPARE TO STORE FASTA DATA
my @fasta_file_names_arr = &get_fasta_file_names;
my %group2seqs_ha;

#my $filename = $ARGV[0];
#my $coding_start_site = $ARGV[1];
#my $coding_stop_site = $ARGV[2];
my $seq_length = 0;

foreach my $fasta_file_name (@fasta_file_names_arr) {
	# Read in the sequence from the file
	my $seq = '';
	my @seqs_arr;
	#my $header = '';
	#my @headers_arr;
	my $seq_num = 0;
	my $last_seq_length;
	
	open(IN, "$fasta_file_name") or die "Could not open file $fasta_file_name\n";
	
#	print "\nGenerating summary nucleotide statistics for $fasta_file_name\...\n";
	
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
	
	if($seq_length == 0) {
		$seq_length = $last_seq_length;
	} elsif($seq_length != $last_seq_length) {
		die "\n\nDIE: The sequences from $fasta_file_name must be the same length as previous files. TERMINATED.\n\n";
	}
	
	# Get name for group from $fasta_file_name
	my $group_name;
	if($fasta_file_name =~ '.fa') {
		$group_name = $`;
	} else {
		die "\n### NEED A FASTA FILE\n\n";
	}
	
	# Add this group of sequences to %group2seqs_ha
	$group2seqs_ha{$group_name} = \@seqs_arr;
	#push(@{$product_coordinates_harr{$product_name}->{product_coord_arr}},@product_coord_arr);
}

#print "\nSequence length is $seq_length\n";

# Now I want to go site-by-site and (1) determine the number of defined nucleotides at
# each sites; (2) calculate the consensus/majority nucleotide among those defined using 
# some threshold like 90%; (3) in cases where majority nucleotides differ between groups,
# make note of the sites, nucleotide identities, and coverage.


my @groups = sort keys %group2seqs_ha;

# Print header to output file
print "site";
foreach my $group (@groups) {
	print "\t";
	print "$group\_maj_nt\t";
	print "$group\_maj_nt_count\t";
	print "$group\_defined_count\t";
	print "$group\_maj_nt_freq";
}
print "\n";

# For each nucleotide
for(my $i = 0; $i < $seq_length; $i++) {
	my %group2maj_nt;
	
	# For each group (FASTA file)
	foreach my $group (@groups) {
		my @this_group_seqs = @{$group2seqs_ha{$group}};
		my %nt_count;
		
		# Count number of each nucleotide (A, C, G, T) at this site in this group
		foreach my $this_seq (@this_group_seqs) {
			my $this_nt = substr($this_seq,$i,1);
#			print "\nNT IS: $this_nt\n";
			$nt_count{$this_nt}++;
		}
		
		# Determine majority nucleotide
		my $maj_nt = 'N';
		my $maj_count = 0;
		my $defined_count = 0;
		
		# Check if the two highest have the same counts?	
		
		# Count the number of nucleotides that are defined, i.e., not 'N' or '-'	
		foreach my $observed_nt (keys %nt_count) {
			if($observed_nt ne 'N' && $observed_nt ne '-') {
				$defined_count += $nt_count{$observed_nt};
				
				# If it's got the highest count observed so far, make it the major nucleotide
				if($nt_count{$observed_nt} > $maj_count) {
					$maj_nt = $observed_nt;
					$maj_count = $nt_count{$observed_nt};
				}
			}
		}
		
		# Define frequency of major nucleotide in this group if there are defined nucleotides
		my $maj_nt_freq = 0;
		if($defined_count > 0) {
			$maj_nt_freq = $maj_count / $defined_count;
		}
		
		# If a group was all N's, this will remain N and 0
		$group2maj_nt{$group}->{maj_nt} = $maj_nt;
		$group2maj_nt{$group}->{maj_count} = $maj_count;
		$group2maj_nt{$group}->{defined_count} = $defined_count;
		$group2maj_nt{$group}->{maj_nt_freq} = $maj_nt_freq; # use this later to filter
		
#		print "\nsite $i group $group maj_nt $maj_nt with $maj_count / $defined_count = $maj_nt_freq\n";
		
	}
	
	my $prev_group_maj_nt = '';
	my $maj_nt_diff_flag = 0;
	
	# Find out if there's a difference, but don't consider if it's an N or gap (-)
	foreach my $group (@groups) {
		my $this_group_maj_nt = $group2maj_nt{$group}->{maj_nt};
		if($this_group_maj_nt ne 'N' && $this_group_maj_nt ne '-') {
			
			my $this_group_defined_count = $group2maj_nt{$group}->{defined_count};
			my $this_group_maj_nt_freq = $group2maj_nt{$group}->{maj_nt_freq};
			
			# Make sure it fits our minimum freq and seq coverage criteria
			if($this_group_maj_nt_freq >= $min_variant_maj_nt_freq && $this_group_defined_count >= $min_site_coverage) {
				if($prev_group_maj_nt eq '') { # haven't seen a defined nt yet
					$prev_group_maj_nt = $group2maj_nt{$group}->{maj_nt};
				} elsif($prev_group_maj_nt ne $group2maj_nt{$group}->{maj_nt}) {
		#			print "\nsite index $i has disparate major nucleotides among groups.\n";
					$maj_nt_diff_flag = 1;
					last;
				}
			}# else {
			#	print "\n### excluded a variant in $group index $i because criteria not satisfied\n\n";
			#}
		}
	}
	
	# If there are disparate major nucleotides among groups, print out site info
	if($maj_nt_diff_flag == 1) {
		my $site_num = $i+1;
		
		print "$site_num";
		
		foreach my $group (@groups) {
##			print "\t" . $group2maj_nt{$group}->{maj_nt};
##			print "\t" . $group2maj_nt{$group}->{maj_count};
##			print "\t" . $group2maj_nt{$group}->{defined_count};
##			print "\t" . $group2maj_nt{$group}->{maj_nt_freq};
			
#			if($site_num == 83) {
#				print "\n### GROUP $group\n";
#				print "maj_nt=" . $group2maj_nt{$group}->{maj_nt} . "\n";
#				print "maj_count=" . $group2maj_nt{$group}->{maj_count} . "\n";
#				print "defined_count=" . $group2maj_nt{$group}->{defined_count} . "\n";
#				print "maj_nt_freq=" . $group2maj_nt{$group}->{maj_nt_freq} . "\n";
#			}
			
			# ONLY PRINT if THIS group's frequency matches the frequency criterion
			if($group2maj_nt{$group}->{maj_nt_freq} >= $min_variant_maj_nt_freq && $group2maj_nt{$group}->{defined_count} >= $min_site_coverage) {
				print "\t" . $group2maj_nt{$group}->{maj_nt};
				print "\t" . $group2maj_nt{$group}->{maj_count};
				print "\t" . $group2maj_nt{$group}->{defined_count};
				print "\t" . $group2maj_nt{$group}->{maj_nt_freq};
			} else {
				print "\t";
				print "\t";
				print "\t";
				print "\t";
			}
			
		}
		print "\n";
		
	} # there's a difference in the major nucleotides
	
} # end last site


## my @A_counts;
## my @C_counts;
## my @G_counts;
## my @T_counts;
## 
## # INITIALIZE HASH VALUES
## for(my $site_id=1; $site_id<=length($seq); $site_id++) { # for each site
## 	my $curr_site_index = $site_id - 1;
## 	push(@A_counts,0);
## 	push(@C_counts,0);
## 	push(@G_counts,0);
## 	push(@T_counts,0);
## }
## 
## for(my $seq_id=1; $seq_id <= scalar(@seqs_arr); $seq_id++) { # for each sequence
## 	my $curr_seq_index = $seq_id - 1;
## 	my $curr_seq = $seqs_arr[$curr_seq_index];
## 	#print "\ncurr_seq: $curr_seq\n";
## 	
## 	for(my $site_id=1; $site_id<=length($curr_seq); $site_id++) { # for each site
## 		my $curr_site_index = $site_id - 1;
## 		my $curr_nt = substr($curr_seq,$curr_site_index,1); # substr is 0-based indexing. VAR,OFFSET,LEN
## 		#print "\ncurr_nt: $curr_nt";
## 		
## 		if($curr_nt eq 'A') {
## 			$A_counts[$curr_site_index]++;
## 			#$C_counts[$curr_site_index]+=0;
## 			#$G_counts[$curr_site_index]+=0;
## 			#$T_counts[$curr_site_index]+=0;
## 		} elsif($curr_nt eq 'C') {
## 			#$A_counts[$curr_site_index]+=0;
## 			$C_counts[$curr_site_index]++;
## 			#$G_counts[$curr_site_index]+=0;
## 			#$T_counts[$curr_site_index]+=0;
## 		} elsif($curr_nt eq 'G') {
## 			#$A_counts[$curr_site_index]+=0;
## 			#$C_counts[$curr_site_index]+=0;
## 			$G_counts[$curr_site_index]++;
## 			#$T_counts[$curr_site_index]+=0;
## 		} elsif($curr_nt eq 'T') {
## 			#$A_counts[$curr_site_index]+=0;
## 			#$C_counts[$curr_site_index]+=0;
## 			#$G_counts[$curr_site_index]+=0;
## 			$T_counts[$curr_site_index]++;
## 		}
## 	}	
## }
## 
## # Do calculations and output to summary file
## my $summary_file_name;
## if($filename =~ '.fasta') {
## 	$summary_file_name = $` . "_site_summary.txt";
## } elsif($filename =~ '.fa') {
## 	$summary_file_name = $` . "_site_summary.txt";
## } else {
## 	$summary_file_name = "fasta_site_summary.txt";
## }
## 
## my $total_transitions;
## my $total_transversions;
## 
## my %num_vars2num_nts; # keep track of how many variant nucleotides
## 
## #my %non_singleton_sites;
## my $non_singleton_transitions;
## my $non_singleton_transversions;
## 
## # For non-singleton polymorphic sites with only 2 nucleotides present
## my $non_singleton_2_transitions;
## my $non_singleton_2_transversions;
## my %non_singleton_2_codon_pos_counts;
## my %non_singleton_multiV_codon_pos_counts;
## 
## # For singleton polymorphic sites with only 2 nucleotides present
## my $singleton_2_transitions;
## my $singleton_2_transversions;
## my %singleton_2_codon_pos_counts;
## my %singleton_multiV_codon_pos_counts;
## 
## open(OUT, ">>$summary_file_name");
## print OUT "site\tA\tA_prop\tC\tC_prop\tG\tG_prop\tT\tT_prop\n";
## for(my $site_id=1; $site_id<=scalar(@A_counts); $site_id++) { # for each site
## 	my $site_index = $site_id - 1;
## 	
## 	# The following method should automatically ignore gaps (-)
## 	my $A = $A_counts[$site_index];
## 	my $C = $C_counts[$site_index];
## 	my $G = $G_counts[$site_index];
## 	my $T = $T_counts[$site_index];
## 	
## 	my $total = ($A + $C + $G + $T);
## 	
## #	print "My total is $total\n";
## 	
## 	if($total != $seq_num) {
## 		#die "\n\nNucleotide sum ($total) at site $site_id does not match number of sequences ($seq_num).\n".
## 		#	"A=$A\nC=$C\nG=$G\nT=$T\nTERMINATED\n\n";
## 		
## #		print "\n\nNucleotide sum ($total) at site $site_id does not match number of sequences ($seq_num).\n".
## #			"A=$A\nC=$C\nG=$G\nT=$T\nALIGNMENT GAP?\n\n";
## 	}
## 	
## 	if($total > 0) {
## 		
## #		die "\nMAKE SURE SITES WITH ONLY N's WON'T OFFPUT THE ALIGNMENT ANALYSIS\n".
## #			"(I think we're okay, but check)\n\n";
## 		
## 		# Which is the dominant nucleotide?
## 		my $maj_nt; # a variant nucleotide may have fixed
## 		my $maj_nt_count = 0;
## 		if($A > $maj_nt_count) {
## 			$maj_nt_count = $A;
## 			$maj_nt = 'A';
## 		}
## 		if($C > $maj_nt_count) {
## 			$maj_nt_count = $C;
## 			$maj_nt = 'C';
## 		} 
## 		if($G > $maj_nt_count) {
## 			$maj_nt_count = $G;
## 			$maj_nt = 'G';
## 		} 
## 		if($T > $maj_nt_count) {
## 			$maj_nt_count = $T;
## 			$maj_nt = 'T';
## 		}
## 		
## 		my $min_sum = $total - $maj_nt_count;
## 		
## #		if($min_sum > 1) {
## #			$non_singleton_sites{$site_id} = 1;
## #		}
## 
## 		# How many different variable nucleotides? Want to find only 2 (1 shared variant)
## 		my $num_alt_nts = 0;
## 		if($A > 0) {
## 			$num_alt_nts++;
## 		}
## 		if($C > 0) {
## 			$num_alt_nts++;
## 		}
## 		if($G > 0) {
## 			$num_alt_nts++;
## 		}
## 		if($T > 0) {
## 			$num_alt_nts++;
## 		}
## 		
## 		$num_vars2num_nts{$min_sum}->{$num_alt_nts}++;
## 		
## 		my $A_prop = ($A / $total);
## 		my $C_prop = ($C / $total);
## 		my $G_prop = ($G / $total);
## 		my $T_prop = ($T / $total);
## 		
## 		# Compute transitions and transversions
## 		my $num_pw_comps = (($total ** 2) - $total) / 2;
## 		my $AC = $A * $C;
## 		my $AG = $A * $G;
## 		my $AT = $A * $T;
## 		my $CG = $C * $G;
## 		my $CT = $C * $T;
## 		my $GT = $G * $T;
## 		
## 		my $transitions = $AG + $CT;
## 		my $transversions = $AC + $AT + $CG + $GT;
## 		
## 		$total_transitions += $transitions;
## 		$total_transversions += $transversions;
## 		
## 		
## 		
## 		
## 		if($min_sum == 1) { # SINGLETONS	
## 			if($num_alt_nts == 2) { # singleton AND one-variant
## 				$singleton_2_transitions += $transitions;
## 				$singleton_2_transversions += $transversions;
## 				
## 				if($coding_start_site && $coding_stop_site) { # second and third argts
## 					if($site_id < $coding_start_site) {
## 						$singleton_2_codon_pos_counts{UTR_5}++;
## 					} elsif($site_id <= $coding_stop_site) {
## 						my $codon_position = (($site_id - $coding_start_site) % 3) + 1;
## 						$singleton_2_codon_pos_counts{$codon_position}++;
## 					} else {
## 						$singleton_2_codon_pos_counts{UTR_3}++;
## 					}
## 				}
## 				
## 			} else { # singleton, multi-variant
## 				die "\n\nThis should never occur\n\n";
## 			}
## 		}
## 		
## 		
## 		
## 		
## 		
## 		if($min_sum > 1) { # non-singleton (more than one variable sequence)
## 			$non_singleton_transitions += $transitions;
## 			$non_singleton_transversions += $transversions;
## 			
## 			if($num_alt_nts == 2) { # non-singleton AND one-variant
## 				$non_singleton_2_transitions += $transitions;
## 				$non_singleton_2_transversions += $transversions;
## 				
## 				if($coding_start_site && $coding_stop_site) { # second and third argts
## 					if($site_id < $coding_start_site) {
## 						$non_singleton_2_codon_pos_counts{UTR_5}++;
## 					} elsif($site_id <= $coding_stop_site) { # CODING SITE
## 						my $codon_position = (($site_id - $coding_start_site) % 3) + 1;
## 						$non_singleton_2_codon_pos_counts{$codon_position}++;
## 						
## 						# DETERMINE CODON HERE FOR SYNONYMOUS COUNT?
## 						
## 					} else {
## 						$non_singleton_2_codon_pos_counts{UTR_3}++;
## 					}
## 				}
## 				
## 				
## 			} else { # non-singleton, multi-variant
## 			
## 			
## 				if($coding_start_site && $coding_stop_site) { # second and third argts
## 					if($site_id < $coding_start_site) {
## 						$non_singleton_multiV_codon_pos_counts{UTR_5}++;
## 					} elsif($site_id <= $coding_stop_site) {
## 						my $codon_position = (($site_id - $coding_start_site) % 3) + 1;
## 						$non_singleton_multiV_codon_pos_counts{$codon_position}++;
## 					} else {
## 						$non_singleton_multiV_codon_pos_counts{UTR_3}++;
## 					}
## 				}
## 				
## 				
## 				
## 				
## 				
## 			}
## 		}
## 		
## 		my $transi_transv_ratio;
## 		if($transversions > 0) {
## 			$transi_transv_ratio = $transitions / $transversions;
## 		} else {
## 			$transi_transv_ratio = '*';
## 		}
## 		
## #		print "\nThere are $transitions transitions and $transversions transversions\n".
## #			"The ratio is: $transi_transv_ratio\n\n";
## 		
## 		print OUT "$site_id\t$A\t$A_prop\t$C\t$C_prop\t$G\t$G_prop\t$T\t$T_prop\n";
## 		
## 	}
## 	
## }
## close OUT;
## 
## my $total_transi_transv_ratio;
## if($total_transversions > 0) {
## 	$total_transi_transv_ratio = $total_transitions / $total_transversions;
## } else {
## 	$total_transi_transv_ratio = '*';
## }
## 
## my $non_singleton_transi_transv_ratio;
## if($non_singleton_transversions > 0) {
## 	$non_singleton_transi_transv_ratio = $non_singleton_transitions / $non_singleton_transversions;
## } else {
## 	$non_singleton_transi_transv_ratio = '*';
## }
## 
## my $non_singleton_2_transi_transv_ratio;
## if($non_singleton_2_transversions > 0) {
## 	$non_singleton_2_transi_transv_ratio = $non_singleton_2_transitions / $non_singleton_2_transversions;
## } else {
## 	$non_singleton_2_transi_transv_ratio = '*';
## }
## 
## my $singleton_2_transi_transv_ratio;
## if($singleton_2_transversions > 0) {
## 	$singleton_2_transi_transv_ratio = $singleton_2_transitions / $singleton_2_transversions;
## } else {
## 	$singleton_2_transi_transv_ratio = '*';
## }
## 
## print "\nThere are $total_transitions transitions and $total_transversions transversions\n".
## 	"The ratio is: $total_transi_transv_ratio\nThe proportion is: " .
## 	$total_transitions / ($total_transitions + $total_transversions) . "\n\n";
## 			
## print "\nFor non-singletons, there are $non_singleton_transitions transitions and $non_singleton_transversions transversions\n".
## 	"The ratio is: $non_singleton_transi_transv_ratio\nThe proportion is: " .
## 	$non_singleton_transitions / ($non_singleton_transitions + $non_singleton_transversions) . "\n\n";
## 	
## print "\nFor non-singletons one-variant sites, there are $non_singleton_2_transitions transitions and $non_singleton_2_transversions transversions\n".
## 	"The ratio is: $non_singleton_2_transi_transv_ratio\nThe proportion is: " .
## 	$non_singleton_2_transitions / ($non_singleton_2_transitions + $non_singleton_2_transversions) . "\n".
## 	"These occur at the following sites:\n5\'NCR: ".$non_singleton_2_codon_pos_counts{UTR_5}.
## 	"\nCodon pos 1: ".$non_singleton_2_codon_pos_counts{1}.
## 	"\nCodon pos 2: ".$non_singleton_2_codon_pos_counts{2}.
## 	"\nCodon pos 3: ".$non_singleton_2_codon_pos_counts{3}.
## 	"\n3\'NCR: ".$non_singleton_2_codon_pos_counts{UTR_3}.
## 	"\nNon-singleton, MULTI-variant occur at the following sites:\n5\'NCR: ".$non_singleton_multiV_codon_pos_counts{UTR_5}.
## 	"\nCodon pos 1: ".$non_singleton_multiV_codon_pos_counts{1}.
## 	"\nCodon pos 2: ".$non_singleton_multiV_codon_pos_counts{2}.
## 	"\nCodon pos 3: ".$non_singleton_multiV_codon_pos_counts{3}.
## 	"\n3\'NCR: ".$non_singleton_multiV_codon_pos_counts{UTR_3}.
## 	"\n\n";
## 
## # We don't really need these data right now
## print "\nFor singletons one-variant sites, there are $singleton_2_transitions transitions and $singleton_2_transversions transversions\n".
## 	"The ratio is: $singleton_2_transi_transv_ratio\nThe proportion is: " .
## 	$singleton_2_transitions / ($singleton_2_transitions + $singleton_2_transversions) . "\n".
## 	"These occur at the following sites:\n5\'NCR: ".$singleton_2_codon_pos_counts{UTR_5}.
## 	"\nCodon pos 1: ".$singleton_2_codon_pos_counts{1}.
## 	"\nCodon pos 2: ".$singleton_2_codon_pos_counts{2}.
## 	"\nCodon pos 3: ".$singleton_2_codon_pos_counts{3}.
## 	"\n3\'NCR: ".$singleton_2_codon_pos_counts{UTR_3}.
## 	"\n\n";
## 
## print "num_var_seqs\tnum_distinct_nts\tcount\n";
## foreach my $num_var_seqs (sort {$a <=> $b} keys %num_vars2num_nts) {
## 	foreach my $num_distinct_nts (sort {$a <=> $b} keys %{$num_vars2num_nts{$num_var_seqs}}) {
## 		print "$num_var_seqs\t$num_distinct_nts\t" . $num_vars2num_nts{$num_var_seqs}->{$num_distinct_nts} . "\n";
## 	}
## }

print "\n\nCOMPLETED.\n\n";

#########################################################################################
# Obtains all file names in current directory ending in .fa and/or .fasta
sub get_fasta_file_names { 
	my @fasta_file_names = glob "*.fa";
	my @other_fasta_file_names;
	
	#if (scalar(@fasta_file_names) == 0) {
	#	@fasta_file_names = glob "*.fasta";
	#} else {
		@other_fasta_file_names = glob "*.fasta";
		push (@fasta_file_names,@other_fasta_file_names);
	#}
	
	#print "\n\n@fasta_file_names\n\n";
	
	if (scalar(@fasta_file_names) == 0) {
		#chdir('SNPGenie_Results');
		#open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		## FILE | PRODUCT | SITE | CODON | WARNING
		#print ERROR_FILE "N/A\tN/A\tN/A\t".
		#	"No FASTA (.fa or .fasta) files in directory. SNPGenie terminated.\n";
		#close ERROR_FILE;
		#chdir('..');
		
		die "\n\n## WARNING: There are no .fa or .fasta files. SNPGenie terminated.\n\n";
	}
	
	#print "\n\n@fasta_file_names\n\n";
	return 	@fasta_file_names;
}
