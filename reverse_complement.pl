#! /usr/bin/perl

# Takes in one FASTA file as an argument and returns the reverse complement in a
# similarly named file

use strict;
use warnings;

# Read in the sequence from the command line
my $seq = $ARGV[0] or die "\n### PROVIDE A NUCLEOTIDE SEQUENCE\n\n";
chomp($seq);
$seq = uc($seq); # uc returns uppercase
$seq =~ tr/U/T/; # no RNA

my $rev_seq = reverse($seq);
my $rev_com_seq = $rev_seq;
$rev_com_seq =~ tr/ACGT/TGCA/;

print "$rev_com_seq\n";

exit;
