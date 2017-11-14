#! /usr/bin/perl

# Extracts all variable codons from .ann.vcf produced from snpEff with -formatEff option,
# i.e., includes field "CODON: Codon change (e.g. 'ggT/ggG')"
# see "SnpEff 'EFF' fields" header at http://snpeff.sourceforge.net/SnpSift.html

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

# DATE CREATED: September 24, 2017
# AUTHOR: Chase W. Nelson
# CONTACT: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
#use warnings;
use Data::Dumper;

my @files_to_process = glob "*.ann.vcf";

print "\nAnalyzing " . scalar(@files_to_process) . " \*.ann.vcf files!\n";

my @codon_metadata;
my @ref_codons;
my @alt_codons;

# Loop through files
foreach my $file (@files_to_process) {

	# Extract Proband ID
	my $Proband;
	if($file =~/\.ann\.vcf/) { 
		$Proband = $`;
	} else {
		die "\nArgument must be a .vcf file\n\n";
	}
	
	# Open file
	open(IN_ANN_FILE, $file) or die "\n## Cannot open $file. TERMINATED.\n\n";
	while(<IN_ANN_FILE>) {
		chomp;
		
		my $line = $_;
		
		# Identify protein-coding mutations; store ref codon, alt codon, and metadata
		if($line =~ /\|([acgtACGT]{3})\/([acgtACGT]{3})\|/) {
			my $ref_codon = uc($1);
			my $alt_codon = uc($2);
			
			my $chromosome;
			my $position;
			if($line =~ /^(chr[XY\d]+)\t(\d+)/) {
				$chromosome = $1;
				$position = $2;
			}
			
			push(@codon_metadata,"$Proband\t$chromosome\t$position");
			push(@ref_codons,$ref_codon);
			push(@alt_codons,$alt_codon);
		}
		
	}
	
	close IN_ANN_FILE;

}

# Write all ref/alt codons and metadata to three files
open(REF_CODONS_OUTPUT,">>ref_codons_seq.txt");
open(ALT_CODONS_OUTPUT,">>alt_codons_seq.txt");
open(CODON_METADATA_OUTPUT,">>codon_metadata.txt");

print CODON_METADATA_OUTPUT "codon_num\tproband\tchr\tposition\tref\talt\n"; 

for(my $i=0; $i<scalar(@codon_metadata); $i++) {
	my $codon_num = $i+1;
	my $ref_codon = $ref_codons[$i];
	my $alt_codon = $alt_codons[$i];
	my $metadata = $codon_metadata[$i];
	
	print REF_CODONS_OUTPUT "$ref_codon";
	print ALT_CODONS_OUTPUT "$alt_codon";
	print CODON_METADATA_OUTPUT "$codon_num\t$metadata\t$ref_codon\t$alt_codon\n";
	
}

close REF_CODONS_OUTPUT;
close ALT_CODONS_OUTPUT;
close CODON_METADATA_OUTPUT;

print "\nDone extracting.\n###############################\n\n";

