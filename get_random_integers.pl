#! /usr/bin/perl

# Takes in as arguments:
#	[0] number of random integer values wanted;
#	[1] bottom of range;
#	[2] top of range
# OUTPUTS: specified number of random integers in range.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# get_random_integers.pl 10 3 999
#########################################################################################

# Copyright (C) 2017 Chase W. Nelson
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# the Gerstner Family Foundation at the American Museum of Natural History, New York.

my $length = $ARGV[0] or die "\nProvide a first argument, the number of numbers\n\n";
my $start_range = $ARGV[1] or die "\nProvide a second argument, the beginning of range\n\n";
my $stop_range = $ARGV[2] or die "\nProvide a third argument, the end of range\n\n";

my @numbers;

for (my $i = $start_range; $i <= $stop_range; $i++) {
	push @numbers, $i;
}

#$number_to_add = $numbers[rand @numbers];

#print "\n\nMy numbers are: @numbers\n\n";

my @selection;

for (my $i = 0; $i < $length; $i++) {
	$number_to_add = $numbers[rand @numbers];
	push @selection, $number_to_add;
}

#print "\n\nMy random numbers are: @selection\n\n";

foreach (@selection) {
	print "$_\n";
}