#!/usr/bin/perl

# Read sstat output and write input data file for dominating set
# NOTE1: sstat repors similarities, not distances
# NOTE2: sstat reports each pair only once and not distance to itself

use strict;
use warnings;

my %adjt = ();

my $T = 0.000015;

my $n = 0;
my $n_minus = 0;

while (my $line = <STDIN>) {
    chomp $line;

    next if ($line =~ m/^#/);

    #print STDERR "$line\n";

    my ($query, $target, $smax, $ssum, $imax, $bmaxp, $alpha_a, $alpha_b) = split /\t/, $line;
    
    if ($ssum > $T) {
	$adjt{$query}->{$target} = 1; # close
	$adjt{$target}->{$query} = 1;
    } else {
	$adjt{$query}->{$target} = 0; # not close
	$adjt{$target}->{$query} = 0;	
    }

    $adjt{$query}->{$query} = 1;     # just put something as placeholders, 
    $adjt{$target}->{$target} = 1;
}



print "set NODES      :=";
for my $node (sort keys %adjt) {
    print "\t\"$node\"";
}
print ";\n";

# print adjacency tables with self edges
print "param AdjTable:";
for my $node (sort keys %adjt) {
    print "\t\"$node\"";
}
print ":=\n";
my $li = scalar(keys %adjt)-1;
for my $node1 (sort keys %adjt) {
    print "\"$node1\"";
    
    for my $node2 (sort keys %adjt) {
	print "\t", $adjt{$node1}->{$node2};

    }
    print "\n";
}
print ";\n";




