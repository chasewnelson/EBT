#! /usr/bin/perl

# INPUT: Takes in one FASTA file with multiple aligned NUCLEOTIDE seqs as an argument.
# OUTPUTS: one TAB-delimited file to the working directory, and brief summary 
#	statistics to the Terminal.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# aligned_fasta2site_nt_freqs.pl <aligned_seqs.fasta>
#########################################################################################

# Copyright (C) 2016 Chase W. Nelson
# DATE CREATED: November 2016
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural 
#	History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# 	the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use warnings;

my $filename = $ARGV[0];
my $coding_start_site = $ARGV[1];
my $coding_stop_site = $ARGV[2];

# Read in the sequence from the file
my $seq = '';
my @seqs_arr;
#my $header = '';
#my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN, "$filename") or die "Could not open FASTA file $filename\n";

print "\n\nGenerating summary nucleotide statistics...\n\n";

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

$seq_num = scalar(@seqs_arr);

my @A_counts;
my @C_counts;
my @G_counts;
my @T_counts;

# INITIALIZE HASH VALUES
for(my $site_id=1; $site_id<=length($seq); $site_id++) { # for each site
	my $curr_site_index = $site_id - 1;
	push(@A_counts,0);
	push(@C_counts,0);
	push(@G_counts,0);
	push(@T_counts,0);
}

for(my $seq_id=1; $seq_id <= scalar(@seqs_arr); $seq_id++) { # for each sequence
	my $curr_seq_index = $seq_id - 1;
	my $curr_seq = $seqs_arr[$curr_seq_index];
	#print "\ncurr_seq: $curr_seq\n";
	
	for(my $site_id=1; $site_id<=length($curr_seq); $site_id++) { # for each site
		my $curr_site_index = $site_id - 1;
		my $curr_nt = substr($curr_seq,$curr_site_index,1); # substr is 0-based indexing. VAR,OFFSET,LEN
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

# Do calculations and output to summary file
my $summary_file_name;
if($filename =~ '.fasta') {
	$summary_file_name = $` . "_site_summary.txt";
} elsif($filename =~ '.fa') {
	$summary_file_name = $` . "_site_summary.txt";
} else {
	$summary_file_name = "fasta_site_summary.txt";
}

# Also prepare a "SNP report"
my $snp_file_name;
if($filename =~ '.fasta') {
	$snp_file_name = $` . "_snp_summary.txt";
} elsif($filename =~ '.fa') {
	$snp_file_name = $` . "_snp_summary.txt";
} else {
	$snp_file_name = "fasta_snp_summary.txt";
}

# Do
my $total_transitions;
my $total_transversions;

my %num_vars2num_nts; # keep track of how many variant nucleotides

#my %non_singleton_sites;
my $non_singleton_transitions;
my $non_singleton_transversions;

# For non-singleton polymorphic sites with only 2 nucleotides present
my $non_singleton_2_transitions;
my $non_singleton_2_transversions;
my %non_singleton_2_codon_pos_counts;
my %non_singleton_multiV_codon_pos_counts;

# For singleton polymorphic sites with only 2 nucleotides present
my $singleton_2_transitions;
my $singleton_2_transversions;
my %singleton_2_codon_pos_counts;
my %singleton_multiV_codon_pos_counts;

open(OUT, ">>$summary_file_name");
print OUT "site\tA\tA_prop\tC\tC_prop\tG\tG_prop\tT\tT_prop\n";

open(OUT_SNP, ">>$snp_file_name");
print OUT_SNP "site\tA\tA_prop\tC\tC_prop\tG\tG_prop\tT\tT_prop\n";

for(my $site_id=1; $site_id<=scalar(@A_counts); $site_id++) { # for each site
	my $site_index = $site_id - 1;
	
	# The following method should automatically ignore gaps (-)
	my $A = $A_counts[$site_index];
	my $C = $C_counts[$site_index];
	my $G = $G_counts[$site_index];
	my $T = $T_counts[$site_index];
	
	my $total = ($A + $C + $G + $T);
	
#	print "My total is $total\n";
	
	if($total != $seq_num) {
		#die "\n\nNucleotide sum ($total) at site $site_id does not match number of sequences ($seq_num).\n".
		#	"A=$A\nC=$C\nG=$G\nT=$T\nTERMINATED\n\n";
		
#		print "\n\nNucleotide sum ($total) at site $site_id does not match number of sequences ($seq_num).\n".
#			"A=$A\nC=$C\nG=$G\nT=$T\nALIGNMENT GAP?\n\n";
	}
	
	if($total > 0) {
		
#		die "\nMAKE SURE SITES WITH ONLY N's WON'T OFFPUT THE ALIGNMENT ANALYSIS\n".
#			"(I think we're okay, but check)\n\n";
		
		# Which is the dominant nucleotide?
		my $maj_nt; # a variant nucleotide may have fixed
		my $maj_nt_count = 0;
		if($A > $maj_nt_count) {
			$maj_nt_count = $A;
			$maj_nt = 'A';
		}
		if($C > $maj_nt_count) {
			$maj_nt_count = $C;
			$maj_nt = 'C';
		} 
		if($G > $maj_nt_count) {
			$maj_nt_count = $G;
			$maj_nt = 'G';
		} 
		if($T > $maj_nt_count) {
			$maj_nt_count = $T;
			$maj_nt = 'T';
		}
		
		my $min_sum = $total - $maj_nt_count;
		
#		if($min_sum > 1) {
#			$non_singleton_sites{$site_id} = 1;
#		}

		# How many different variable nucleotides? Want to find only 2 (1 shared variant)
		my $num_alt_nts = 0;
		if($A > 0) {
			$num_alt_nts++;
		}
		if($C > 0) {
			$num_alt_nts++;
		}
		if($G > 0) {
			$num_alt_nts++;
		}
		if($T > 0) {
			$num_alt_nts++;
		}
		
		$num_vars2num_nts{$min_sum}->{$num_alt_nts}++;
		
		my $A_prop = ($A / $total);
		my $C_prop = ($C / $total);
		my $G_prop = ($G / $total);
		my $T_prop = ($T / $total);
		
		# Compute transitions and transversions
		my $num_pw_comps = (($total ** 2) - $total) / 2;
		my $AC = $A * $C;
		my $AG = $A * $G;
		my $AT = $A * $T;
		my $CG = $C * $G;
		my $CT = $C * $T;
		my $GT = $G * $T;
		
		my $transitions = $AG + $CT;
		my $transversions = $AC + $AT + $CG + $GT;
		
		$total_transitions += $transitions;
		$total_transversions += $transversions;
		
		
		
		
		if($min_sum == 1) { # SINGLETONS	
			if($num_alt_nts == 2) { # singleton AND one-variant
				$singleton_2_transitions += $transitions;
				$singleton_2_transversions += $transversions;
				
				if($coding_start_site && $coding_stop_site) { # second and third argts
					if($site_id < $coding_start_site) {
						$singleton_2_codon_pos_counts{UTR_5}++;
					} elsif($site_id <= $coding_stop_site) {
						my $codon_position = (($site_id - $coding_start_site) % 3) + 1;
						$singleton_2_codon_pos_counts{$codon_position}++;
					} else {
						$singleton_2_codon_pos_counts{UTR_3}++;
					}
				}
				
			} else { # singleton, multi-variant
				die "\n\nThis should never occur\n\n";
			}
		}
		
		
		
		
		
		if($min_sum > 1) { # non-singleton (more than one variable sequence)
			$non_singleton_transitions += $transitions;
			$non_singleton_transversions += $transversions;
			
			if($num_alt_nts == 2) { # non-singleton AND one-variant
				$non_singleton_2_transitions += $transitions;
				$non_singleton_2_transversions += $transversions;
				
				if($coding_start_site && $coding_stop_site) { # second and third argts
					if($site_id < $coding_start_site) {
						$non_singleton_2_codon_pos_counts{UTR_5}++;
					} elsif($site_id <= $coding_stop_site) { # CODING SITE
						my $codon_position = (($site_id - $coding_start_site) % 3) + 1;
						$non_singleton_2_codon_pos_counts{$codon_position}++;
						
						# DETERMINE CODON HERE FOR SYNONYMOUS COUNT?
						
					} else {
						$non_singleton_2_codon_pos_counts{UTR_3}++;
					}
				}
				
				
			} else { # non-singleton, multi-variant
			
			
				if($coding_start_site && $coding_stop_site) { # second and third argts
					if($site_id < $coding_start_site) {
						$non_singleton_multiV_codon_pos_counts{UTR_5}++;
					} elsif($site_id <= $coding_stop_site) {
						my $codon_position = (($site_id - $coding_start_site) % 3) + 1;
						$non_singleton_multiV_codon_pos_counts{$codon_position}++;
					} else {
						$non_singleton_multiV_codon_pos_counts{UTR_3}++;
					}
				}
				
				
				
				
				
			}
		}
		
		my $transi_transv_ratio;
		if($transversions > 0) {
			$transi_transv_ratio = $transitions / $transversions;
		} else {
			$transi_transv_ratio = '*';
		}
		
#		print "\nThere are $transitions transitions and $transversions transversions\n".
#			"The ratio is: $transi_transv_ratio\n\n";
		
		print OUT "$site_id\t$A\t$A_prop\t$C\t$C_prop\t$G\t$G_prop\t$T\t$T_prop\n";
		
	}
	
}
close OUT;

my $total_transi_transv_ratio;
if($total_transversions > 0) {
	$total_transi_transv_ratio = $total_transitions / $total_transversions;
} else {
	$total_transi_transv_ratio = '*';
}

my $non_singleton_transi_transv_ratio;
if($non_singleton_transversions > 0) {
	$non_singleton_transi_transv_ratio = $non_singleton_transitions / $non_singleton_transversions;
} else {
	$non_singleton_transi_transv_ratio = '*';
}

my $non_singleton_2_transi_transv_ratio;
if($non_singleton_2_transversions > 0) {
	$non_singleton_2_transi_transv_ratio = $non_singleton_2_transitions / $non_singleton_2_transversions;
} else {
	$non_singleton_2_transi_transv_ratio = '*';
}

my $singleton_2_transi_transv_ratio;
if($singleton_2_transversions > 0) {
	$singleton_2_transi_transv_ratio = $singleton_2_transitions / $singleton_2_transversions;
} else {
	$singleton_2_transi_transv_ratio = '*';
}

print "\nThere are $total_transitions transitions and $total_transversions transversions\n".
	"The ratio is: $total_transi_transv_ratio\nThe proportion is: " .
	$total_transitions / ($total_transitions + $total_transversions) . "\n\n";
			
print "\nFor non-singletons, there are $non_singleton_transitions transitions and $non_singleton_transversions transversions\n".
	"The ratio is: $non_singleton_transi_transv_ratio\nThe proportion is: " .
	$non_singleton_transitions / ($non_singleton_transitions + $non_singleton_transversions) . "\n\n";
	
print "\nFor non-singletons one-variant sites, there are $non_singleton_2_transitions transitions and $non_singleton_2_transversions transversions\n".
	"The ratio is: $non_singleton_2_transi_transv_ratio\nThe proportion is: " .
	$non_singleton_2_transitions / ($non_singleton_2_transitions + $non_singleton_2_transversions) . "\n".
	"These occur at the following sites:\n5\'NCR: ".$non_singleton_2_codon_pos_counts{UTR_5}.
	"\nCodon pos 1: ".$non_singleton_2_codon_pos_counts{1}.
	"\nCodon pos 2: ".$non_singleton_2_codon_pos_counts{2}.
	"\nCodon pos 3: ".$non_singleton_2_codon_pos_counts{3}.
	"\n3\'NCR: ".$non_singleton_2_codon_pos_counts{UTR_3}.
	"\nNon-singleton, MULTI-variant occur at the following sites:\n5\'NCR: ".$non_singleton_multiV_codon_pos_counts{UTR_5}.
	"\nCodon pos 1: ".$non_singleton_multiV_codon_pos_counts{1}.
	"\nCodon pos 2: ".$non_singleton_multiV_codon_pos_counts{2}.
	"\nCodon pos 3: ".$non_singleton_multiV_codon_pos_counts{3}.
	"\n3\'NCR: ".$non_singleton_multiV_codon_pos_counts{UTR_3}.
	"\n\n";

# We don't really need these data right now
print "\nFor singletons one-variant sites, there are $singleton_2_transitions transitions and $singleton_2_transversions transversions\n".
	"The ratio is: $singleton_2_transi_transv_ratio\nThe proportion is: " .
	$singleton_2_transitions / ($singleton_2_transitions + $singleton_2_transversions) . "\n".
	"These occur at the following sites:\n5\'NCR: ".$singleton_2_codon_pos_counts{UTR_5}.
	"\nCodon pos 1: ".$singleton_2_codon_pos_counts{1}.
	"\nCodon pos 2: ".$singleton_2_codon_pos_counts{2}.
	"\nCodon pos 3: ".$singleton_2_codon_pos_counts{3}.
	"\n3\'NCR: ".$singleton_2_codon_pos_counts{UTR_3}.
	"\n\n";

print "num_var_seqs\tnum_distinct_nts\tcount\n";
foreach my $num_var_seqs (sort {$a <=> $b} keys %num_vars2num_nts) {
	foreach my $num_distinct_nts (sort {$a <=> $b} keys %{$num_vars2num_nts{$num_var_seqs}}) {
		print "$num_var_seqs\t$num_distinct_nts\t" . $num_vars2num_nts{$num_var_seqs}->{$num_distinct_nts} . "\n";
	}
}

print "\n\nCOMPLETED.\n\n";
