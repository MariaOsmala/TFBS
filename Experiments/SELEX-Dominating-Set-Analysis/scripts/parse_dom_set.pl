#!/usr/bin/perl
# Parse the dominating set from the gplsol output

# Usage
# $ perl parse_dom_set.pl <../solutions/SELEX_v3_1_TOMTOM_005_min_dom.sol

use strict;
use warnings;

my $obj = 0;

my $node_state = 0;
my $current_node = "";
my $sol_size = 0;
my $n = 0;

while (my $line = <STDIN>) {
    chomp $line; 

    if ($line =~ m/^Objective:\s+z = (\d+)/) {
	$obj = $1;
	print STDERR "Objective: $1\n";

    } elsif ($line =~ m/v\[([\w_\-]+)\]/) {
	$current_node = $1;
	#print STDERR "node: $current_node\n";
	$node_state = 1;

    } elsif ($node_state == 1) {
	# on two line, this one should be the line telling the value
	if ($line =~ m/\*\s+(\d)\s+\d/) {
	    my $val = $1;
	    #print STDERR "$current_node: $val\n";

	    if ($val == 1) {
		print "$current_node\n";
		$sol_size++;
	    }
	    $n++;
	    $node_state = 0;
	    die "Error: Weird node value" if ($val != 0 && $val != 1);
	} else {
	    die "Error: Could not parse value for node $current_node";
	}
    } else {
	#print STDERR "$line\n";
    }
}

print STDERR "Parsed nodes: $n\n";
print STDERR "Objective $obj, parsed solution size $sol_size\n";
print STDERR "NO MATCH!\n" if ($obj != $sol_size);
