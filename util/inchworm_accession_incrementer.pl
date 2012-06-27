#!/usr/bin/env perl

use strict;
use warnings;

my $x = 0;

while (<>) {
	if (/>/) {
		$x++;
		s/>/>s$x;/;
	}
	print;
}


exit(0);
