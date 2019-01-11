#! /usr/bin/perl

# Creates reverse complement versions of all 3 SNPGenie input files.

# Copyright (C) 2015 Chase W. Nelson

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

# AUTHOR: Chase W. Nelson
# CONTACT1: nelsoncw@email.sc.edu
# CONTACT2: cwnelson88@gmail.com
# AFFILIATION1: Austin L. Hughes lab, University of South Carolina (Columbia, SC, USA)
# AFFILIATION2: Wen-Hsiung Li lab, Academia Sinica (Taipei, Taiwan)

# ACKNOWLEDGMENTS: written by C.W.N. with support from a National Science Foundation
# Graduate Research Fellowship (DGE-0929297), a National Science Foundation East Asian 
# and Pacific Summer Institutes Fellowship, and a University of South Carolina 
# Presidential Fellowship.

use strict;
#use warnings;
use IO::Handle;

if(scalar @ARGV != 1) {
	die "\n\n## WARNING: This script requires exactly 1 ".
		"argument:\n## (1) A '+' strand FASTA file\n\n## For example: ".
		"fasta2revcom.pl my_sequence.fasta\n\n";
}

my $fasta_file = $ARGV[0];


#########################################################################################
# STORE / REVCOM / PRINT EACH SEQUENCE
my $seq = '';
my $header = '';
my $seq_num = 0;

my $output_fasta = 'revcom.fasta';
if($fasta_file =~/([\w\.]+)\.\w+/) { 
	$output_fasta = "$1\_revcom.fasta";
}

open(OUT_FASTA, ">>$output_fasta");
open(IN_FASTA, "$fasta_file") or die "Could not open FASTA file $fasta_file\n";

print "\nReading in sequence from FASTA...\n";

while(<IN_FASTA>) {
	chomp;
	
	if(/>/) {
		$seq_num++;
		
		if($seq_num == 1) { # just print the header
			$header = $_;
			print OUT_FASTA "$header\_revcom\n";
			
		} else { # not first sequence; header already printed
			my $next_header = $_;
			
			# Process and print the newly finished sequence
			my $this_seq_length = length($seq);
			print "Sequence $header is of length $this_seq_length\n";
			
			$seq = uc($seq);
			$seq = reverse($seq);
			$seq =~ tr/ACGT/TGCA/;
			
			# PRINT REVCOM SEQ
			while(length($seq) > 0) {
				print OUT_FASTA '' . substr($seq, 0, 50, '') . "\n";
			}
			
			# Ready new sequence
			$seq = '';
			$header = $next_header;
			
			# Print next header
			print OUT_FASTA "\n$header\_revcom\n";
		}
		
	} else {
		$seq .= $_;
	}
}

close IN_FASTA;

# Take care of last sequence
# Process and print the newly finished sequence
my $this_seq_length = length($seq);
print "Sequence $header is of length $this_seq_length\n";

$seq = uc($seq);
$seq = reverse($seq);
$seq =~ tr/ACGT/TGCA/;

# PRINT REVCOM SEQ
while(length($seq) > 0) {
	print OUT_FASTA '' . substr($seq, 0, 50, '') . "\n";
}

close OUT_FASTA;


print "\n## REVERSE COMPLEMENT FASTA has been written to:\n".
	"## FASTA: $output_fasta\n\n";


exit;

