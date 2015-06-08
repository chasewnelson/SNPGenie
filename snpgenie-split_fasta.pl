#! /usr/bin/perl

use strict;
use warnings;

# Title: splitFasta.pl
# Author: Chase W. Nelson
# Affiliation1: Austin L. Hughes lab, University of South Carolina (Columbia, SC 29208, USA)
# Affiliation2: Wen-Hsiung Li lab, Academia Sinica (Taipei, Taiwan)
# Contact1: cwnelson88@gmail.com
# Contact2: nelsoncw@email.sc.edu

# Date of creation: 09 August 2013
# Last update: 09 August 2013

# Explanation: This script takes a FASTA file as input and creates n FASTA files,
# one for each of the n individual sequences in the input FASTA. Thus, it is meant
# for separating a FASTA file with n sequences into n FASTA files, one for each sequence.

# Check that an argument is given
# User-given arguments are placed in the array @ARGV, indexed with 0

if(! $ARGV[0]) {
	die "\nAn argument must be supplied: seqFile.fasta";
}

my $fastaFile = $ARGV[0];
my $whichSeq = 1;

# Woei-Fu's request
my $newFilePrefix;
if($fastaFile =~/\.fasta/) { 
	$newFilePrefix = $`;
} elsif($fastaFile =~/\.txt/) {
	$newFilePrefix = $`;
}

open(REF_FASTA_FILE,"$fastaFile"); # Reference FASTA sequence

while (<REF_FASTA_FILE>) { # For each line of the FASTA
	#chomp;
	
	if($_ =~/^>(\w+)/) {
		if($whichSeq == 1) {
			open(CURR_OUTFILE,">>$newFilePrefix\_$1\.fasta");
			$whichSeq++;
			print CURR_OUTFILE "$_";
		} else {
			close CURR_OUTFILE;
			open(CURR_OUTFILE,">>$newFilePrefix\_$1\.fasta");
			print CURR_OUTFILE "$_";
		}
	} else {
		print CURR_OUTFILE "$_";
	}
}

close CURR_OUTFILE;
close REF_FASTA_FILE;