#!/usr/bin/env perl

use 5.014;
use utf8;
use warnings;
use strict;

my $input_file_len = $ARGV[0];
my $num_samples = $ARGV[1];

if(@ARGV == 3){
    srand(int($ARGV[2]));
}else{
    srand(20161201);
}

# prepare random index
my @rnd = ();
foreach(1..$num_samples){
    push(@rnd, int(rand($input_file_len)));
}
my @index = reverse sort @rnd;

# select lines
my $next = pop @index;
my $line = 0;

while(<STDIN>){
    chomp;
    if($line == $next){
	print $_, "\n";
	if(@index > 0){
	    $next = pop @index;
	}else{
	    last;
	}
    }
    $line++;
}
