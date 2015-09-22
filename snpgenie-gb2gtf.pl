#! /usr/bin/perl

# Takes a Genbank file as its 1 ARGV argument, and creates a GTF file for its products.
# In this instantiation, we are concerned about CDS annotations only
# If reverse complements exists, those annotations are put in another file

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

use strict;
#use warnings;
use IO::Handle;

my $genbank_file_name = $ARGV[0];

# Generate new file name prefix
my $new_file_prefix;
if($genbank_file_name =~/\.gbk/) { 
	$new_file_prefix = $`;
} elsif($genbank_file_name =~/\.gb/) { 
	$new_file_prefix = $`;
} else {
	$new_file_prefix = "inFile";
}

my $new_file_name = $new_file_prefix . ".gtf";
my $new_file_compl_name = $new_file_prefix . "_rev_compl.gtf";
my $seq_length = 0;

# For the standard annotations
my $curr_ORF = '';
my $curr_start = 0;
my $curr_stop = 0;
my $curr_start2 = 0;
my $curr_stop2 = 0;

my %ORF_info_hh;
my @ORF_arr;

# For the COMPLEMENT annotations (i.e., antisense strand)
my $is_compl = 0;
my $curr_compl_ORF = '';
my $curr_compl_start = 0;
my $curr_compl_stop = 0;
my $curr_compl_start2 = 0;
my $curr_compl_stop2 = 0;

my %ORF_compl_info_hh;
my @ORF_compl_arr;

# Name
my $seen_label = 0;

open(GBK_FILE, $genbank_file_name) or die "Could not open file $genbank_file_name\n";

while(<GBK_FILE>) {
	if($_ =~ /^\s*LOCUS/) {
		if($' =~ /(\d+)\s*bp/) {
			$seq_length = $1;
		}
	} elsif($_ =~ /^\s*CDS\s+(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	} elsif($_ =~ /^\s*CDS\s+<(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	} elsif($_ =~ /^\s*CDS\s+join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$curr_start2 = $3;
		$curr_stop2 = $4;
		$seen_label = 0;
	} elsif($_ =~ /^\s*CDS\s+complement\((\d+)\.\.(\d+)\)/) { # for complement
		$curr_compl_start = $1;
		$curr_compl_stop = $2;
		$seen_label = 0;
		$is_compl = 1;
	} 

	# Now find the name of the feature, in order of preference, as:
	# [1] locus_tag; [2] label; and [3] product
	if($seen_label == 0 && $curr_start != 0 && $curr_compl_start == 0) {
		if($_ =~ /\s*\/locus_tag="([\w\s\.']+ [\w\s\.']+)"/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/locus_tag="([\w\s\.']+)"/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF ='';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/locus_tag=([\w\s\.']+)/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label="([\w\s\.']+ [\w\s\.']+)"/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label="([\w\s\.']+)"/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label=([\w\s\.']+)/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product="([\w\s\.']+ [\w\s\.']+)"/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product="([\w\s\.']+)"/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product=([\w\s\.']+)/) {
			if($curr_start2 == 0) {
				$curr_ORF = $1;
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = 0;
				$curr_stop = 0;
				$curr_start2 = 0;
				$curr_stop2 = 0;
				
				$seen_label = 1;
			}
		}
	} elsif($seen_label == 0 && $curr_start == 0 && $curr_compl_start != 0) { # COMPLEMENT
		if($_ =~ /\s*\/locus_tag="([\w\s\.']+ [\w\s\.']+)"/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/locus_tag="([\w\s\.']+)"/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF ='';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/locus_tag=([\w\s\.']+)/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				chomp($curr_compl_ORF);
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				chomp($curr_compl_ORF);
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label="([\w\s\.']+ [\w\s\.']+)"/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label="([\w\s\.']+)"/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label=([\w\s\.']+)/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				chomp($curr_compl_ORF);
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				chomp($curr_compl_ORF);
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product="([\w\s\.']+ [\w\s\.']+)"/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product="([\w\s\.']+)"/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product=([\w\s\.']+)/) {
			if($curr_compl_start2 == 0) {
				$curr_compl_ORF = $1;
				chomp($curr_compl_ORF);
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				
				$seen_label = 1;
			} else {
				$curr_compl_ORF = $1;
				chomp($curr_compl_ORF);
				push(@ORF_compl_arr,$curr_compl_ORF);
				$ORF_compl_info_hh{$curr_compl_ORF}->{start} = $curr_compl_start;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop} = $curr_compl_stop;
				$ORF_compl_info_hh{$curr_compl_ORF}->{start2} = $curr_compl_start2;
				$ORF_compl_info_hh{$curr_compl_ORF}->{stop2} = $curr_compl_stop2;
				
				$curr_compl_ORF = '';
				$curr_compl_start = 0;
				$curr_compl_stop = 0;
				$curr_compl_start2 = 0;
				$curr_compl_stop2 = 0;
				
				$seen_label = 1;
			}
		}
	} # end of COMPLEMENT
} # end of current GenBank file line

close GBK_FILE;

# SENSE STRAND to GTF file
my @ORF_arr_sorted = sort (@ORF_arr);

open(OUTFILE,">>$new_file_name");

foreach my $curr_product (@ORF_arr_sorted) {
	my $this_start = $ORF_info_hh{$curr_product}->{start};
	my $this_stop = $ORF_info_hh{$curr_product}->{stop};
	
	print OUTFILE "$genbank_file_name\tCLC\tCDS\t$this_start\t$this_stop\t\.\t\+\t0\tgene_id \"$curr_product\";\n";
	
	if (exists $ORF_info_hh{$curr_product}->{start2}) {
		my $this_start2 = $ORF_info_hh{$curr_product}->{start2};
		my $this_stop2 = $ORF_info_hh{$curr_product}->{stop2};
		
		print OUTFILE "$genbank_file_name\tCLC\tCDS\t$this_start2\t$this_stop2\t\.\t\+\t0\tgene_id \"$curr_product\";\n";
	}
}	

close OUTFILE;

# ANTISENSE STRAND to GTF file, to a separate file as the + strand
#if($is_compl == 1) {
#	my $genbank_compl_file_name = $new_file_prefix . "_revcompl.gb";
#	my @ORF_compl_arr_sorted = sort (@ORF_compl_arr);
#	
#	open(OUTFILE,">>$new_file_compl_name");
#	
#	foreach my $curr_product (@ORF_compl_arr_sorted) {
#		my $this_start = $ORF_compl_info_hh{$curr_product}->{start};
#		my $this_stop = $ORF_compl_info_hh{$curr_product}->{stop};
#		my $feature_length = ($this_stop - $this_start + 1);
#		
#		my $offset = ($seq_length - $this_stop);
#		my $rev_compl_start = ($offset + 1);
#		my $rev_compl_stop = ($rev_compl_start + $feature_length - 1);
#		
#		print OUTFILE "$genbank_compl_file_name\tCLC\tCDS\t$rev_compl_start\t$rev_compl_stop\t\.\t\+\t0\tgene_id \"$curr_product\";\n";
#		
#		if (exists $ORF_compl_info_hh{$curr_product}->{start2}) {
#			my $this_start2 = $ORF_compl_info_hh{$curr_product}->{start2};
#			my $this_stop2 = $ORF_compl_info_hh{$curr_product}->{stop2};
#			my $feature_length2 = ($this_stop2 - $this_start2 + 1);
#			
#			my $offset2 = ($seq_length - $this_stop2);
#			my $rev_compl_start2 = ($offset2 + 1);
#			my $rev_compl_stop2 = ($rev_compl_start2 + $feature_length2 - 1);
#		
#			print OUTFILE "$genbank_compl_file_name\tCLC\tCDS\t$rev_compl_start2\t$rev_compl_stop2\t\.\t\+\t0\tgene_id \"$curr_product\";\n";
#		}
#	}	
#	
#	close OUTFILE;
#}

# ANTISENSE STRAND to GTF file
if($is_compl == 1) {
	my $genbank_compl_file_name = $new_file_prefix . "_revcompl.gb";
	my @ORF_compl_arr_sorted = sort (@ORF_compl_arr);
	
	open(OUTFILE,">>$new_file_name");
	
	foreach my $curr_product (@ORF_compl_arr_sorted) {
		my $this_start = $ORF_compl_info_hh{$curr_product}->{start};
		my $this_stop = $ORF_compl_info_hh{$curr_product}->{stop};
		#my $feature_length = ($this_stop - $this_start + 1);
		
		#my $offset = ($seq_length - $this_stop);
		#my $rev_compl_start = ($offset + 1);
		#my $rev_compl_stop = ($rev_compl_start + $feature_length - 1);
		
		print OUTFILE "$genbank_file_name\tCLC\tCDS\t$this_start\t$this_stop\t\.\t\-\t0\tgene_id \"$curr_product\";\n";
		
		if (exists $ORF_compl_info_hh{$curr_product}->{start2}) {
			my $this_start2 = $ORF_compl_info_hh{$curr_product}->{start2};
			my $this_stop2 = $ORF_compl_info_hh{$curr_product}->{stop2};
			#my $feature_length2 = ($this_stop2 - $this_start2 + 1);
			
			#my $offset2 = ($seq_length - $this_stop2);
			#my $rev_compl_start2 = ($offset2 + 1);
			#my $rev_compl_stop2 = ($rev_compl_start2 + $feature_length2 - 1);
		
			print OUTFILE "$genbank_file_name\tCLC\tCDS\t$this_start2\t$this_stop2\t\.\t\-\t0\tgene_id \"$curr_product\";\n";
		}
	}	
	
	close OUTFILE;
}

# END