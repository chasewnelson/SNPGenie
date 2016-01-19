#! /usr/bin/perl

# Converts a GFF (argument 1) to a GTF file for use with SNPGenie.

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

if(scalar @ARGV != 1) {
	die "\n\n## WARNING: The SNPGenie script gff2gtf needs exactly 1 ".
		"argument:\n## A GFF file with '+' strand data relative to the ".
		"reference sequence against which SNPs were called.".
		"\n## Only the \"ID\" tag will be used to identify the gene name, e.g., \"ID=GENE_001\"".
		"\n\n## For example: ".
		"snpgenie-gff2gtf.pl my_cds_file.gff\n\n";
}

my $gff_file_nm = $ARGV[0];

# Generate new file name names
my $new_gtf_file_nm;
if($gff_file_nm =~/\.gff/) { 
	$new_gtf_file_nm = $` . "_gff_converted.gtf";
} elsif($gff_file_nm =~/gff\.txt/) { 
	$new_gtf_file_nm = $` . "_gff_converted.gtf";
} else {
	die "\nFirst and only argument must be a .gff file\n\n";
}

if(-e $new_gtf_file_nm) {
	die "\n## $new_gtf_file_nm already exists in this directory; delete before proceeding\n\n";
}

print "\n## Converting $gff_file_nm to $new_gtf_file_nm for SNPGenie...\n";

#my %all_ORF_info; # {name}->{all_attributes}

open(OUTFILE,">>$new_gtf_file_nm");
open (CURRINFILE, $gff_file_nm);
while (<CURRINFILE>) {
	if(/^##FASTA/ || /^\>/) {
		last;
	} else {
		unless(/^#/) {
			if(/CDS/) {
				chomp;
				
				my @line_arr = split(/\t/,$_);
				
				my $seqname = $line_arr[0];
				my $source = $line_arr[1]; 
				my $feature = $line_arr[2]; 
				my $start = $line_arr[3]; 
				my $end = $line_arr[4]; 
				my $score = $line_arr[5]; 
				my $strand = $line_arr[6]; 
				my $frame = $line_arr[7];
				my $group = $line_arr[8];
				
				my $gene_id;
				if($group =~ /ID=([\w\s\.']+)/) {
					$gene_id = $1;
				} elsif($group =~ /ID=([\w\s\.']+ [\w\s\.']+)/) {
					$gene_id = $1;
				} 
				
#				$all_ORF_info{$gene_id}->{$seqname} = $line_arr[0];
#				$all_ORF_info{$gene_id}->{$source} = $line_arr[1]; 
#				$all_ORF_info{$gene_id}->{$feature} = $line_arr[2]; 
#				$all_ORF_info{$gene_id}->{$start} = $line_arr[3]; 
#				$all_ORF_info{$gene_id}->{$end} = $line_arr[4]; 
#				$all_ORF_info{$gene_id}->{$score} = $line_arr[5]; 
#				$all_ORF_info{$gene_id}->{$strand} = $line_arr[6]; 
#				$all_ORF_info{$gene_id}->{$frame} = $line_arr[7];
				
				#	if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
				#		if($_ =~/gene_id \"gene\:([\w\s\.']+)\"/) {
				#			$products_hash{$1} = 1;
				#		} elsif($_ =~ /gene_id \"([\w\s\.']+ [\w\s\.']+)\"/) {
				#			$products_hash{$1} = 1;
				#		} elsif($_ =~/gene_id \"([\w\s\.']+)\"/) {
				#			$products_hash{$1} = 1;
				#		} 
				#	}
				
				my $this_line = "". $seqname . "\t" . $source . "\t" . $feature . "\t" . $start . 
					"\t" . $end . "\t" . $score . "\t" . $strand . "\t" . $frame . "\t" . 
					"gene_id \"$gene_id\"\;";
				
				print OUTFILE "$this_line\n";
			}
		}
	}
}
close CURRINFILE;
close OUTFILE;


print "\n## GTF for SNPGenie file has been written to $new_gtf_file_nm\n\n";