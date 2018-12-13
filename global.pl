#!/usr/bin/env perl 

#
# dynamic programming global alignment
# usage = ./global.pl <seq> <seq> <scores_file>
# bcb 5300 - algorithms in comutational biology
#
# outside sources:
#   - teaching materials by dr. erin chambers
#   - teaching materials by dr. llana lareau
#

use strict;
use warnings;

my $s1 = $ARGV[0];
my $s2 = $ARGV[1];
my $matrix = $ARGV[2];

unless ($s1 and $s2 and $matrix) { 
	die "global.pl <seq> <seq> <scores_file>\n"
}

my @s1 = split('', $s1);
my @s2 = split('', $s2);

my $matrixref = get_matrix($matrix); 

(my $traceref, my $best) = get_tables(\@s1, \@s2, $matrixref);

(my $align1, my $align2) = trace($traceref, $best->[1], $best->[2], \@s1, \@s2); 

print "$align1\n";
print "$align2";

sub get_tables {
	my $seq1 = $_[0];
	my $seq2 = $_[1];
	my $matrix = $_[2];
	my @scores;
	my @traces;
	my @localbest = (0, 0, 0);

	$scores[0][0] = 0;
	$traces[0][0] = 0;
	for (my $i=1; $i <= ($#$seq1 + 1); $i++) { 
		$scores[$i][0] = -$i;
		$traces[$i][0] = '|'; 
	}

	for (my $j=1; $j <= ($#$seq2 + 1); $j++)   {
		$scores[0][$j] = -$j;   
		$traces[0][$j] = '-'; 
	}

	for (my $i=1; $i <= ($#$seq1 + 1); $i++) {
		for (my $j=1; $j <= ($#$seq2 + 1); $j++) {
			my $x = $scores[$i-1][$j-1] + $matrix->{$seq1->[$i-1]}{$seq2->[$j-1]};

			my $y = $scores[$i-1][$j] - 1; 

			my $z = $scores[$i][$j-1] - 1;

			my @tmp;
			@tmp = sort {$b->[0] <=> $a->[0]} ([$x,'\\'], [$y,'|'], [$z,'-']); 
			@tmp = sort {$b->[0] <=> $a->[0]} ([$x,'\\'], [$y,'|'], [$z,'-']); 

			if ( ($i == ($#$seq1 + 1) or $j == ($#$seq2 + 1)) and ($tmp[0][0] >= $localbest[0]) ) {
				@localbest = ($tmp[0][0], $i, $j);  
			}                                     

			@tmp = sort {$b->[0] <=> $a->[0]} ([$x,'\\'], [$y,'|'], [$z,'-'], [0,0]); 
			if ($tmp[0][0] >= $localbest[0]) {
				@localbest = ($tmp[0][0], $i, $j); 
			}

			$scores[$i][$j] = $tmp[0][0];     
			$traces[$i][$j] = $tmp[0][1];    
		}
	}

	@localbest = ( $scores[$#$seq1 + 1][$#$seq2 + 1], $#$seq1 + 1, $#$seq2 + 1); 
	return (\@traces, \@localbest);
}

sub get_matrix {
	my $file = $_[0];
	my %matrix;     
	my @letters;   
	my $len;      

	open ( MATRIX, "$file" ) or die "Cannot open score matrix.\n";

	while (<MATRIX>) {
		chomp;
		if (/^  \w/) { 
			@letters = /(\w)\s/g;            
			$len = $#letters;                
		} elsif (/^\w/) {
			(my $letter) = /^(\w)/;
			my @scores = /(-?\d+)/g;
			for (my $i=0; $i <= $len; $i++) {
				$matrix{$letter}{$letters[$i]} = $scores[$i];
			}
		}
	}    
	close MATRIX;
	return \%matrix;
}


sub trace {
	my $traces = $_[0];                      
	my $i = $_[1];
	my $j = $_[2];
	my $s1ref = $_[3];
	my $s2ref = $_[4];

	if ($traces->[$i][$j] eq '0') { 
		$align1 = "";
		$align2 = "";
	} elsif ($traces->[$i][$j] eq '\\') {
		my $a = $s1ref->[$i-1];
		my $b = $s2ref->[$j-1];

		($align1, $align2) = trace($traces, $i-1, $j-1, $s1ref, $s2ref);
		$align1 .= $a;
		$align2 .= $b;
	} elsif ($traces->[$i][$j] eq '-') {
		my $a = "-";
		my $b = $s2ref->[$j-1];

		($align1, $align2) = trace($traces, $i, $j-1, $s1ref, $s2ref);
		$align1 .= $a;
		$align2 .= $b;
	} elsif ($traces->[$i][$j] eq '|') {	
		my $a = $s1ref->[$i-1];
		my $b = "-";

		($align1, $align2) = trace($traces, $i-1, $j, $s1ref, $s2ref);
		$align1 .= $a;
		$align2 .= $b;
	}

	return ($align1, $align2);
}
