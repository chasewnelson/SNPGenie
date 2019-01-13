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

if(scalar @ARGV != 2) {
	die "\n\n## WARNING: This script requires exactly 2 ".
		"arguments, in this order:\n## (1) A '+' strand GTF file ".
		"containing both '+' and '–' strand products from the '+' strand point of view; ".
		"and\n## (2) The total sequence length.\n\n## For example: ".
		"gtf2revcom.pl my_cds_file.gtf 10000\n\n";
}

my $gtf_file_nm = $ARGV[0];
my $seq_length = $ARGV[1];

# Generate new file names
print "seq length is $seq_length\n";

my $new_gtf_file_name = &generate_reverse_complement_gtf($gtf_file_nm);

print "\n## Converting $new_gtf_file_name to reverse complement...\n";

#chromosome	lengths, GRCh38
#1	248956422
#2	242193529
#3	198295559
#4	190214555
#5	181538259
#6	170805979
#7	159345973
#8	145138636
#9	138394717
#10	133797422
#11	135086622
#12	133275309
#13	114364328
#14	107043718
#15	101991189
#16	90338345
#17	83257441
#18	80373285
#19	58617616
#20	64444167
#21	46709983
#22	50818468


#NC_000002.12	CCDS	CDS	1414409	1414502	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1423045	1423129	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1433438	1433607	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1436252	1436384	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1453694	1453823	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1456076	1456282	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1484596	1484854	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1487821	1487991	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1493802	1494039	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1495989	1496197	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1496595	1496765	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1503948	1504079	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1516883	1516982	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1540594	1540723	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	1542421	1542474	.	+	0	gene_id "7173";
#NC_000002.12	CCDS	CDS	3545785	3545871	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3547931	3548055	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3548640	3548724	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3549058	3549112	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3550373	3550472	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3552144	3552308	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3556789	3556904	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3558133	3558260	.	-	0	gene_id "246243";
#NC_000002.12	CCDS	CDS	3575610	3575684	.	+	0	gene_id "6201";
#NC_000002.12	CCDS	CDS	3575817	3575888	.	+	0	gene_id "6201";
#NC_000002.12	CCDS	CDS	3576487	3576630	.	+	0	gene_id "6201";
#NC_000002.12	CCDS	CDS	3577710	3577774	.	+	0	gene_id "6201";
#NC_000002.12	CCDS	CDS	3580110	3580260	.	+	0	gene_id "6201";
#NC_000002.12	CCDS	CDS	3580805	3580882	.	+	0	gene_id "6201";


print "\n## REVERSE COMPLEMENT FILES have been written to:\n".
	"## GTF: $new_gtf_file_name\n\n";

# Give NUM LINES CONVERTED tag/warning
#	my @stored_pos_sorted = sort {$a <=> $b} (keys %gtf_line_by_starting_position);
#	my $num_positions = scalar(@stored_pos_sorted);
#	print "\n\n## A total of $num_positions products were processed. Please verify that ".
#		"this is the correct number of products.\n\n";
#	
#	open(OUT_FILE_GTF, ">>$new_gtf_file_name");
#	foreach my $curr_pos (@stored_pos_sorted) {
#		print OUT_FILE_GTF $gtf_line_by_starting_position{$curr_pos};
#	}
#	close OUT_FILE_GTF;


#########################################################################################
sub detect_newline_char {
	my ($curr_filename) = @_;
	my $newline_char;
	
	open (CURRINFILE, $curr_filename);
	while (<CURRINFILE>) {	
		if($_ =~ /\r\n/) {
			$newline_char = "\r\n";
			last;
		} elsif($_ =~ /\r/) {
			$newline_char = "\r";
			last;
		} elsif($_ =~ /\n/) {
			$newline_char = "\n";
			last;
		}
	}
	#seek(CURRINFILE,0,0);
	close CURRINFILE;
	return $newline_char;
}

#########################################################################################
sub get_product_names_vcf {
	my ($gtf_file_nm) = @_;
	my %products_hash;
	open (CURRINFILE, $gtf_file_nm);
	while (<CURRINFILE>) {
		chomp;
		
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			if($_ =~/gene_id \"gene\:([\w\s\.']+)\"/) {
				my $product = $1;
				
				if ((! exists $products_hash{$product}) && ($product ne '')) {
					$products_hash{$product} = 1;
				}
			} elsif($_ =~ /gene_id \"([\w\s\.']+ [\w\s\.']+)\"/) {
				my $product = $1;
				
				if ((! exists $products_hash{$product}) && ($product ne '')) {
					$products_hash{$product} = 1;
				}
			} elsif($_ =~ /gene_id \"([\w\s\.']+)\"/) {
				my $product = $1;
				
				if ((! exists $products_hash{$product}) && ($product ne '')) {
					$products_hash{$product} = 1;
				}
			}
		}
	}
	close CURRINFILE;
	my @product_names = keys %products_hash;
	return @product_names;
}

#########################################################################################
#sub populate_product_information_hh { # Only + strand products, as-is
#	my ($filename) = @_;
#	my %hh_compl_position_info;
#	
#	open(GTF_FILE_AGAIN, "$filename") or die "\nCould not open the GTF file $filename - $!\n\n";
#	while(<GTF_FILE_AGAIN>) {
#		if($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\+\t\d+\tgene_id \"gene\:([\w\s\.\:']+)\"/) { # Line is - strand
#			my $start_pos = $1; # Where the gene itself actually STOPS
#			my $stop_pos = $2; # Where the gene itself actually STARTS
#			my $this_product = $3;
#			
#			if(exists $hh_compl_position_info{$this_product}->{start}) {
#				$hh_compl_position_info{$this_product}->{start_2} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop_2} = $stop_pos;
#			} else {
#				$hh_compl_position_info{$this_product}->{start} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop} = $stop_pos;
#			}
#		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\+\t\d+\tgene_id \"([\w\s\.\:']+ [\w\s\.\:']+)\"/) {
#			my $start_pos = $1; # Where the gene itself actually STOPS
#			my $stop_pos = $2; # Where the gene itself actually STARTS
#			my $this_product = $3;
#			
#			if(exists $hh_compl_position_info{$this_product}->{start}) {
#				$hh_compl_position_info{$this_product}->{start_2} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop_2} = $stop_pos;
#			} else {
#				$hh_compl_position_info{$this_product}->{start} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop} = $stop_pos;
#			}
#		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\+\t\d+\tgene_id \"([\w\s\.\:']+)\"/) {
#			my $start_pos = $1; # Where the gene itself actually STOPS
#			my $stop_pos = $2; # Where the gene itself actually STARTS
#			my $this_product = $3;
#			
#			if(exists $hh_compl_position_info{$this_product}->{start}) {
#				$hh_compl_position_info{$this_product}->{start_2} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop_2} = $stop_pos;
#			} else {
#				$hh_compl_position_info{$this_product}->{start} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop} = $stop_pos;
#			}
#		# NOW, IN CASE transcript_id comes first
#		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\+\t\d+\ttranscript_id \"[\w\s\.\:']+\"; gene_id \"gene\:([\w\s\.\:']+)\"/) {
#			my $start_pos = $1; # Where the gene itself actually STOPS
#			my $stop_pos = $2; # Where the gene itself actually STARTS
#			my $this_product = $3;
#			
#			if(exists $hh_compl_position_info{$this_product}->{start}) {
#				$hh_compl_position_info{$this_product}->{start_2} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop_2} = $stop_pos;
#			} else {
#				$hh_compl_position_info{$this_product}->{start} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop} = $stop_pos;
#			}
#		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\+\t\d+\ttranscript_id \"[\w\s\.\:']+\"; gene_id \"([\w\s\.\:']+ [\w\s\.\:']+)\"/) {
#			my $start_pos = $1; # Where the gene itself actually STOPS
#			my $stop_pos = $2; # Where the gene itself actually STARTS
#			my $this_product = $3;
#			
#			if(exists $hh_compl_position_info{$this_product}->{start}) {
#				$hh_compl_position_info{$this_product}->{start_2} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop_2} = $stop_pos;
#			} else {
#				$hh_compl_position_info{$this_product}->{start} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop} = $stop_pos;
#			}
#		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\+\t\d+\ttranscript_id \"[\w\s\.\:']+\"; gene_id \"([\w\s\.\:']+)\"/) {
#			my $start_pos = $1; # Where the gene itself actually STOPS
#			my $stop_pos = $2; # Where the gene itself actually STARTS
#			my $this_product = $3;
#			
#			if(exists $hh_compl_position_info{$this_product}->{start}) {
#				$hh_compl_position_info{$this_product}->{start_2} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop_2} = $stop_pos;
#			} else {
#				$hh_compl_position_info{$this_product}->{start} = $start_pos;
#				$hh_compl_position_info{$this_product}->{stop} = $stop_pos;
#			}
#		}
#	}
#	close GTF_FILE_AGAIN;
#	
#	return %hh_compl_position_info;
#}

#########################################################################################
sub reverse_complement_from_fasta {
	my ($filename) = @_;
	
	# Read in the sequence from the file
	my $seq = '';
	
	open(IN, "$filename") or die "\nCould not open FASTA file $filename\n\n";
	while(<IN>) {
		unless(/>/) {
			chomp;
			$seq .= $_;
		}
	}
	close IN;
	
	my $rev_seq = reverse($seq);
	my $rev_compl = $rev_seq;
	$rev_compl =~ tr/ACGT/TGCA/;
	
	return $rev_compl;
}

#########################################################################################
sub generate_reverse_complement_gtf {
	my ($filename) = @_;
	
	my $new_gtf_file_name;
	if($filename =~/\.gtf/) { 
		$new_gtf_file_name = $` . "_revcom.gtf";
	} else {
		#$new_gtf_file_name = "gtf_revcom.gtf";
		die "\nSecond argument must be a .gtf file\n\n";
	}
	
	if(-e $new_gtf_file_name) {
		die "\n## $new_gtf_file_name already exists in this directory; delete before proceeding\n\n";
	}
	
	my @gtf_lines;
	my $lines_total = 0;
	
	open(GTF_FILE, "$filename") or die "\nCould not open the GTF file $filename - $!\n\n";
	while(<GTF_FILE>) {
		if($_ =~ /CDS\t(\d+)\t(\d+)/) {
			my $old_start = $1; # Where the gene itself actually STOPS
			my $old_stop = $2; # Where the gene itself actually STARTS
			
			# Find the coordinates from the revcom point of view
			my $this_start = $seq_length - $old_stop + 1;
			my $this_stop = $seq_length - $old_start + 1;
			
			# Replace old sites and strand with new
			#$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\+/CDS\t$this_start\t$this_stop\t\.\t\-/;
			#$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\-/CDS\t$this_start\t$this_stop\t\.\t\+/;
			$_ =~ s/CDS\t$old_start\t$old_stop\t([\d\.]+)\t\+/CDS\t$this_start\t$this_stop\t\1\t\-/;
			$_ =~ s/CDS\t$old_start\t$old_stop\t([\d\.]+)\t\-/CDS\t$this_start\t$this_stop\t\1\t\+/;
			
			push(@gtf_lines,$_);
			$lines_total++;
		} 
	}
	close GTF_FILE;
	
	# Give NUM LINES CONVERTED tag/warning
	print "\n## A total of $lines_total lines were processed. Please verify that ".
		"this is the correct number of lines.\n";
	
	open(OUT_FILE_GTF, ">>$new_gtf_file_name");
	foreach my $gtf_line (@gtf_lines) {
		print OUT_FILE_GTF $gtf_line;
	}
	close OUT_FILE_GTF;
	
	return $new_gtf_file_name;
}

