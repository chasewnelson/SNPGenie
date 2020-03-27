#!/usr/bin/env perl

# Creates reverse complement versions of all 3 SNPGenie input files.

# Copyright (C) 2018 Chase W. Nelson

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

# DATE: December, 2018
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

STDOUT->autoflush(1);

if(scalar @ARGV != 2) {
	die "\n\n## WARNING: The SNPGenie script vcfformat1_to_revcom needs exactly 2 ".
		"arguments, in this order:\n## ". 
		"## (1) A '+' strand SNP report in VCF format 1 (SNPGenie descriptions).\n" . 
		"## (2) The exact length of the sequence in the FASTA file.\n". 
		"\n## For example: ".
		"vcfformat1_to_revcom.pl my_snp_report.vcf 248956422\n\n";
}

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

my $curr_snp_report_name = $ARGV[0]; # what we're reading from: the original file
my $seq_length = $ARGV[1];

# Generate new file name names
my $new_vcf_file_name;
if($curr_snp_report_name =~/\.vcf/) { 
	$new_vcf_file_name = $` . "_revcom.vcf";
} else {
	#$new_vcf_file_name = "vcf_revcom.txt";
	die "\nInput must be a .vcf file\n\n";
}

if(-e $new_vcf_file_name) {
	die "\n## $new_vcf_file_name already exists in this directory; delete before proceeding\n\n";
}
print "Provided seq length is $seq_length\n";

print "\n## Converting $curr_snp_report_name to reverse complement format...\n";

# Populate a hash with product information for + strand
#my %hh_compl_position_info = &populate_product_information_hh($gtf_file_nm);
#my @curr_compl_products_ordered_by_start = sort { $hh_compl_position_info{$a}->{start} <=> $hh_compl_position_info{$b}->{start} } keys %hh_compl_position_info;


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
my $header_line;
open (ORIGINAL_SNP_REPORT, $curr_snp_report_name);
FIRST_LOOP: while (<ORIGINAL_SNP_REPORT>) {	
	chomp;
	if($_ =~ /^#\w+/) {
		$header_line = $_;
		if(!($_ =~/\t/)) {
			die "\n\n## WARNING:\n## The SNP Report $curr_snp_report_name is ".
				"not TAB-delimited (\\t), there is only one column, or it does not begin with '#<seqname>'.\n\n";
		}
		last FIRST_LOOP;
	}
}
close ORIGINAL_SNP_REPORT;

$header_line =~ s/^#//;
#print "\nHEADER LINE IS: $header_line\n\n";
my @header_arr = split("\t", $header_line);
#print "\nHEADER ARRAY IS: @header_arr\n\n";

# Now we want to cycle through the file and convert to revcom
open(NEW_VCF_FILE, ">>$new_vcf_file_name");
open(ORIGINAL_SNP_REPORT, $curr_snp_report_name);
while(<ORIGINAL_SNP_REPORT>) {
	chomp;
	my $line = $_;
	
	if(/^#/) { # lines that begins with "##" or "#" are metadata headers
		print NEW_VCF_FILE "$line\n";
	} else {
		my @line_arr = split(/\t/, $line);
		
#		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
#		22	10874444	rs9617549	C	T	100	PASS	AC=981;AF=0.195887;AN=5008;NS=2504;DP=15958;EAS_AF=0.1419;AMR_AF=0.111;AFR_AF=0.3359;EUR_AF=0.165;SAS_AF=0.1544;AA=.|||;VT=SNP;GRCH37_POS=16074105;GRCH37_REF=G;STRAND_FLIP
#		22	11122151	rs5747224	G	A	100	PASS	AC=1358;AF=0.271166;AN=5008;NS=2504;DP=10844;EAS_AF=0.4395;AMR_AF=0.2651;AFR_AF=0.2867;EUR_AF=0.173;SAS_AF=0.182;AA=G|||;VT=SNP;DEPRECATED_RSID=rs572534874;GRCH37_POS=16965224;GRCH37_REF=G;GRCH37_38_REF_STRING_MATCH
#		22	11122417	rs2456393	G	T	100	PASS	AC=1358;AF=0.271166;AN=5008;NS=2504;DP=9564;EAS_AF=0.4395;AMR_AF=0.2651;AFR_AF=0.2867;EUR_AF=0.173;SAS_AF=0.182;AA=G|||;VT=SNP;DEPRECATED_RSID=rs546325452;GRCH37_POS=16965490;GRCH37_REF=G;GRCH37_38_REF_STRING_MATCH
#		22	11123361	rs191556358	C	A	100	PASS	AC=1;AF=0.000199681;AN=5008;NS=2504;DP=7798;EAS_AF=0.001;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0;AA=C|||;VT=SNP;DEPRECATED_RSID=rs531729453;GRCH37_POS=16966434;GRCH37_REF=C;GRCH37_38_REF_STRING_MATCH
		
		my $CHROM = $line_arr[0];
		my $POS = $line_arr[1];
		my $ID = $line_arr[2];
		my $REF = $line_arr[3];
		my $ALT = $line_arr[4];
		my $QUAL = $line_arr[5];
		my $FILTER = $line_arr[6];
		my $INFO = $line_arr[7];
		
		my $line_revcom = "$CHROM\t";
		
		# New REF and ALT
		my $REF_new = $REF;
#		$REF_new = uc($REF_new); # just in case
		$REF_new = reverse($REF_new); # just in case
		$REF_new =~ tr/acgtuACGTU/tgcaaTGCAA/;
		my $REF_length = length($REF_new);
		
		my $ALT_new = $ALT;
#		$ALT_new = uc($ALT_new); # just in case
		$ALT_new = reverse($ALT_new); # just in case
		$ALT_new =~ tr/acgtuACGTU/tgcaaTGCAA/;
		
		# New POS
		my $POS_new = $seq_length - $POS + 1 - ($REF_length - 1);
		
		# Same ID
		$line_revcom .= "$POS_new\t$ID\t$REF_new\t$ALT_new\t$QUAL\t$FILTER\t";
		
		# Finally, new INFO
		my @AA_matches = ($INFO =~ /AA=([\.\w\-]+)/g);
		
		foreach my $AA (@AA_matches) {
			#print "$_ ";
			my $AA_new = $AA;
#			my $AA_new = uc($AA);
			$AA_new = reverse($AA_new); # just in case
			$AA_new =~ tr/acgtuACGTU/tgcaaTGCAA/;
			$INFO =~ s/AA=$AA/AA=$AA_new/;
			#print "POS=$POS FIND=$AA REPLACE=$AA_new\n";
		}
		
#		while($INFO =~ /AA=([\.\w\-]+)/g) {
#			my $AA_new = uc($1);
#			$AA_new = reverse($AA_new); # just in case
#			$AA_new =~ tr/ACGTU/TGCAA/;
#			$INFO =~ s/AA=$1/AA=$AA_new/;
#			print "stuck here $POS ";
#		}
		
		$line_revcom .= "$INFO\t";
		
		# Any other columns after 8 (index 7)?
		for(my $i = 8; $i < scalar(@line_arr); $i++) {
			$line_revcom .= "$line_arr[$i]\t";
		}
		
		chop $line_revcom;
		
		print NEW_VCF_FILE "$line_revcom\n";
	}
}

close ORIGINAL_SNP_REPORT;
close NEW_VCF_FILE;

print "\n## REVERSE COMPLEMENT VCF FILE has been written to:\n".
	"## VCF: $new_vcf_file_name\n\n";

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
#sub generate_reverse_complement_gtf {
#	my ($filename) = @_;
#	
#	my $new_gtf_file_name;
#	if($filename =~/\.gtf/) { 
#		$new_gtf_file_name = $` . "_revcom.gtf";
#	} else {
#		#$new_gtf_file_name = "gtf_revcom.gtf";
#		die "\nSecond argument must be a .gtf file\n\n";
#	}
#	
#	if(-e $new_gtf_file_name) {
#		die "\n## $new_gtf_file_name already exists in this directory; delete before proceeding\n\n";
#	}
#	
#	my %gtf_line_by_starting_position;
#	
#	open(GTF_FILE, "$filename") or die "\nCould not open the GTF file $filename - $!\n\n";
#	while(<GTF_FILE>) {
#		if($_ =~ /CDS\t(\d+)\t(\d+)/) {
#			my $old_start = $1; # Where the gene itself actually STOPS
#			my $old_stop = $2; # Where the gene itself actually STARTS
#			
#			# Find the coordinates from the revcom point of view
#			my $this_start = $seq_length - $old_stop + 1;
#			my $this_stop = $seq_length - $old_start + 1;
#			
#			if(exists $gtf_line_by_starting_position{$this_start}) {
#				die "\n\nTwo products have the same starting position, causing an error.\n".
#					"Please contact script author for a revision.\n\n";
#			} else {
#				$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\+/CDS\t$this_start\t$this_stop\t\.\t\-/;
#				$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\-/CDS\t$this_start\t$this_stop\t\.\t\+/;
#				$gtf_line_by_starting_position{$this_start} = $_;
#			}
#		} 
#	}
#	close GTF_FILE;
#	
#	# Give NUM LINES CONVERTED tag/warning
#	my @stored_pos_sorted = sort {$a <=> $b} (keys %gtf_line_by_starting_position);
#	my $num_positions = scalar(@stored_pos_sorted);
#	print "\n## A total of $num_positions products were processed. Please verify that ".
#		"this is the correct number of products.\n";
#	
#	open(OUT_FILE_GTF, ">>$new_gtf_file_name");
#	foreach my $curr_pos (@stored_pos_sorted) {
#		print OUT_FILE_GTF $gtf_line_by_starting_position{$curr_pos};
#	}
#	close OUT_FILE_GTF;
#	
#	return $new_gtf_file_name;
#}

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
			$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\+/CDS\t$this_start\t$this_stop\t\.\t\-/;
			$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\-/CDS\t$this_start\t$this_stop\t\.\t\+/;
			
			push(@gtf_lines,$_);
			$lines_total++;
		} 
	}
	close GTF_FILE;
	
	# Give NUM LINES CONVERTED tag/warning
	print "\n## A total of $lines_total products were processed. Please verify that ".
		"this is the correct number of products.\n";
	
	open(OUT_FILE_GTF, ">>$new_gtf_file_name");
	foreach my $gtf_line (@gtf_lines) {
		print OUT_FILE_GTF $gtf_line;
	}
	close OUT_FILE_GTF;
	
	return $new_gtf_file_name;
}

