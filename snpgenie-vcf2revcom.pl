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

use strict;
#use warnings;
use IO::Handle;

if(scalar @ARGV != 3) {
	die "\n\n## WARNING: The SNPGenie script vcf2revcom needs exactly 3 ".
		"arguments, in this order:\n## (1) A '+' strand FASTA (.fa or .fasta) file containing the ".
		"reference sequence against which SNPs were called\n## (2) A '+' strand GTF file ".
		"containing both '+' and '–' strand products from the '+' strand point of view; ".
		"and\n## (3) A '+' strand SNP report in VCF format.\n\n## For example: ".
		"snpgenie-vcf2revcom.pl my_snp_report.vcf my_reference_sequence.fasta my_cds_file.gtf\n\n";
}

my $fasta_file_nm = $ARGV[0];
my $gtf_file_nm = $ARGV[1];
my $curr_snp_report_name = $ARGV[2]; # what we're reading from: the original file

# Generate new file name names
my $new_vcf_file_name;
if($curr_snp_report_name =~/\.vcf/) { 
	$new_vcf_file_name = $` . "_revcom.txt";
} else {
	#$new_vcf_file_name = "vcf_revcom.txt";
	die "\nThird argument must be a .vcf file\n\n";
}

if(-e $new_vcf_file_name) {
	die "\n## $new_vcf_file_name already exists in this directory; delete before proceeding\n\n";
}

my $revcom_seq = &reverse_complement_from_fasta($fasta_file_nm);
my $seq_length = length($revcom_seq);
print "seq length is $seq_length\n";

my $new_gtf_file_name = &generate_reverse_complement_gtf($gtf_file_nm);

my $new_fasta_file_name = &generate_reverse_complement_fasta($fasta_file_nm);

print "\n## Converting $curr_snp_report_name to reverse complement SNPGenie format...\n";

# Populate a hash with product information for + strand
#my %hh_compl_position_info = &populate_product_information_hh($gtf_file_nm);
#my @curr_compl_products_ordered_by_start = sort { $hh_compl_position_info{$a}->{start} <=> $hh_compl_position_info{$b}->{start} } keys %hh_compl_position_info;

my $newline_char = &detect_newline_char($curr_snp_report_name);
$/ = $newline_char;
my $newline_type;

if($newline_char eq "\r\n") {
	$newline_type = "Windows (CRLF, \\r\\n\)";
} elsif($newline_char eq "\r") {
	$newline_type = "Mac (CR, \\r\)";
} elsif($newline_char eq "\n") {
	$newline_type = "Unix (LF, \\n\)";
}

print "\n## In file $curr_snp_report_name, the newline type is: $newline_type\n";

my $index_chrom;
my $index_pos;
my $index_id;
my $index_ref;
my $index_alt;
my $index_qual;
my $index_filter;
my $index_info;
my $index_format;
my $index_sample;

my $seen_index_chrom;
my $seen_index_pos;
my $seen_index_id;
my $seen_index_ref;
my $seen_index_alt;
my $seen_index_qual;
my $seen_index_filter;
my $seen_index_info;
my $seen_index_format;
my $seen_index_sample;

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
my $header_line;
open (ORIGINAL_SNP_REPORT, $curr_snp_report_name);
while (<ORIGINAL_SNP_REPORT>) {	
	chomp;
	if($_ =~ /^#(\w+)/) {
		$header_line = $_;
		if(!($_ =~/\t/)) {
			die "\n\n## WARNING:\n## The SNP Report $curr_snp_report_name is ".
				"not TAB-delimited (\\t), or there is only one column.\n\n";
		}
		last;
	}
}
close ORIGINAL_SNP_REPORT;

$header_line =~ s/^#//;
#print "\nHEADER LINE IS: $header_line\n\n";
my @header_arr = split("\t",$header_line);
#print "\nHEADER ARRAY IS: @header_arr\n\n";

# Determine the index of each column
for (my $i=0; $i<scalar(@header_arr); $i++) {
	if ($header_arr[$i] =~ /CHROM/) {
		$index_chrom = $i;
		$seen_index_chrom = 1;
	} elsif ($header_arr[$i] =~ /POS/) { 
		$index_pos = $i;
		$seen_index_pos = 1;
	} elsif ($header_arr[$i] =~ /ID/) {
		$index_id = $i;
		$seen_index_id = 1;
	} elsif ($header_arr[$i] =~ /REF/) {
		$index_ref = $i;
		$seen_index_ref = 1;
	} elsif ($header_arr[$i] =~ /ALT/) { # Will check the value BEGINS with 'SNP'
		$index_alt = $i;
		$seen_index_alt = 1;
	} elsif ($header_arr[$i] =~ /QUAL/) { # Will check the value BEGINS with 'SNP'
		$index_qual = $i;
		$seen_index_qual = 1;
	} elsif ($header_arr[$i] =~ /FILTER/) {
		$index_filter = $i;
		$seen_index_filter = 1;
	} elsif ($header_arr[$i] =~ /INFO/) {
		$index_info = $i;
		$seen_index_info = 1;
	} elsif ($header_arr[$i] =~ /FORMAT/) {
		$index_format = $i;
		$seen_index_format = 1;
		$index_sample = $i+1; # ASSUME the SAMPLE column follows the FORMAT column
	}
}

# DIE if one of the headers have not been seen
if ($seen_index_chrom == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"CHROM\". SNPGenie terminated\n\n";	
} elsif ($seen_index_pos == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"POS\". SNPGenie terminated\n\n";	
} elsif ($seen_index_id == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"ID\". SNPGenie terminated\n\n";	
} elsif ($seen_index_ref == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"REF\". SNPGenie terminated\n\n";	
} elsif ($seen_index_alt == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"ALT\". SNPGenie terminated\n\n";	
} elsif ($seen_index_qual == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"QUAL\". SNPGenie terminated\n\n";	
} elsif ($seen_index_filter == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FILTER\". SNPGenie terminated\n\n";	
} elsif ($seen_index_info == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"INFO\". SNPGenie terminated\n\n";	
} elsif ($seen_index_format == 0) {
	die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FORMAT\". SNPGenie terminated\n\n";	
}

# Product names are not included directly in the VCF. Thus we must use the GTF file
my @product_names_arr = &get_product_names_vcf($gtf_file_nm);
@product_names_arr = sort(@product_names_arr);

# NEED TO BUILD A HASH WITH keys as ALL PRODUCT POSITIONS IN THE GENOME,
# which we will later convert to minus '-' strand coordinates

# Loop through GTF and store product positions with product names in array
my %positions_with_product_hash;
my %products_hash;

# Must be the NEW GTF FILE, with the new + strand products
#print "\nNew GTF file name is: $new_gtf_file_name\n";
open (GTF_INFILE, $new_gtf_file_name);
while (<GTF_INFILE>) {
	chomp;
	
	if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Make sure it's on the + strand
		my $product;
		if($_ =~ /gene_id \"gene\:([\w\s\.']+)\"/) {
			$product = $1;
			
			if ((! exists $products_hash{$product}) && ($product ne '')) {
				$products_hash{$product} = 1;
			}
		} elsif($_ =~ /gene_id \"([\w\s\.']+ [\w\s\.']+)\"/) {
			$product = $1;
			
			if ((! exists $products_hash{$product}) && ($product ne '')) {
				$products_hash{$product} = 1;
			}
		} elsif($_ =~ /gene_id \"([\w\s\.']+)\"/) {
			$product = $1;
			
			if ((! exists $products_hash{$product}) && ($product ne '')) {
				$products_hash{$product} = 1;
			}
		}
		
		if($product) {
			if($_ =~ /CDS\t(\d+)\t(\d+)/) {
				my $start = $1;
				my $stop = $2;
				
				for(my $i = $start; $i <= $stop; $i++) {
					push(@{$positions_with_product_hash{$i}},$product);
				}
			}
		}
	}
}
close GTF_INFILE;

# NEEDED?
my @product_names_arr = sort(keys %products_hash);

# Now we want to cycle through the file again in order to build a SNPGenie (CLC style) 
# version file for output
open(TEMP_FILE_VCF,">>$new_vcf_file_name\_TEMP");
open (ORIGINAL_SNP_REPORT, $curr_snp_report_name);
while (<ORIGINAL_SNP_REPORT>) {
	unless(/^#/) { # lines that begins with "##" or "#" are metadata
		chomp;
		
		my @line_arr = split(/\t/,$_);
		
		# ONLY DO THE THING IF THE GENEIOUS TYPE IS "Polymorphism"; not CDS.
		my $id = $line_arr[$index_id];
		
		# Save the pure records; the GENEIOUS ORDER is:
		# chrom => NC_002516.2
		# pos => 154
		# id => .
		# ref => T
		# alt => C
		# qual => 222
		# filter => PASS
		# info => DP=262;VDB=0.266664;RPB=0.201403;AF1=1;AC1=2;DP4=1,0,219,38;MQ=60;FQ=-282;PV4=1,1,1,0.29
		# format => GT:PL:GQ
		# sample1 => 1/1:255,255,0:99
		
		my $chrom_value = $line_arr[$index_chrom];
		my $pos_value = $line_arr[$index_pos];
		my $id_value = $line_arr[$index_id];
		my $ref_value = $line_arr[$index_ref];
		my $alt_value = $line_arr[$index_alt];
		my $qual_value = $line_arr[$index_qual];
		my $filter_value = $line_arr[$index_filter];
		my $info_value = $line_arr[$index_info];
		my $format_value = $line_arr[$index_format];
		my $sample_value = $line_arr[$index_sample];
		#my $sample1_value = $line_arr[$index_sample1];
		
		# Extract the needed values; the CLC ORDER will be:
		# file_nm
		# ref_pos => 350
		# THERE WILL BE NO CDS HERE, because not Geneious
		# type => SNV OR MNV OR Insertion OR etc.
		# ref => G
		# allele => A
		# count => 270
		# cov => 1313
		# freq => 20.56
		# over_annot => CDS: gag
	
		# REF_POS
		my $ref_pos;
		# In case there's a less-than sign
		if($line_arr[$index_pos] =~ /\<(\d+)/) {
			$ref_pos = $1;
		}
		$ref_pos = $pos_value;
		
		# TYPE
		# Find length, nucleotides involved, etc.
		my $reference_nts = $ref_value;
		my $variant_nts = $alt_value;
		my $reference_length = length($reference_nts);
		#my $variant_length = length($variant_nts);
		my $is_change = 0;

		if($reference_nts ne $variant_nts) {
			$is_change = 1;
		}
		
		# Change the site numbers to the reverse complement strand
		# UPDATE $ref_pos for revcom
		my $ref_pos = ($seq_length - $ref_pos + 1);
		
		# DETERMINE CLC type and WHICH PRODUCTS OVERLAP THIS SITE ON THIS STRAND
		my $clc_type;
		my @this_site_products;
		if($reference_length == 1) {
			$clc_type = 'SNV';
			if(exists $positions_with_product_hash{$ref_pos}) {
				@this_site_products = @{$positions_with_product_hash{$ref_pos}};
			}
		} elsif($reference_length > 1) {
			$clc_type = 'MNV';
			#print "\nThere is a multi-nucleotide variant at $ref_pos;\n".
			#	"this is not fully supported for VCF\n\n";
			# ELABORATE HERE for sites i+1, i+2, etc.
			# Actually, we can add lines to those sites, IF VCF even does MNVs
			if(exists $positions_with_product_hash{$ref_pos}) {
				@this_site_products = @{$positions_with_product_hash{$ref_pos}};
			}
		}
		
		my $variant1;
		my $variant2;
		my $variant3;
		
		# CHANGE nucleotides to complements
		$reference_nts =~ tr/ACGT/TGCA/;
		$variant_nts =~ tr/ACGT/TGCA/;
		
		# EXTRACT variants 
		if($variant_nts =~ /(\w+),(\w+),(\w+)/) {
			$variant1 = $1;
			$variant2 = $2;
			$variant3 = $3;
		} elsif($variant_nts =~ /(\w+),(\w+)/) {
			$variant1 = $1;
			$variant2 = $2;
		} elsif($variant_nts =~ /\w+/) {
			$variant1 = $&;
		}
		
		# Check all variants are equal lengths, if exist
		my $equal_lengths_variants = 0;
		if($variant3) {
			if($reference_length == length($variant1) && length($variant1) == length($variant2) &&
				length($variant2) == length($variant3)) {
				$equal_lengths_variants = 1;
			}
		} elsif($variant2) {
			if($reference_length == length($variant1) && length($variant1) == length($variant2)) {
				$equal_lengths_variants = 1;
			}
		} elsif($reference_length == length($variant1)) {
			$equal_lengths_variants = 1;
		}
		
		# IF the variant(s) are equal lengths, then we can proceed to procees and
		# print to file
		my $product_entry = '';
		if($equal_lengths_variants) {
			# Double check reference nucleotide identity
			if(($reference_length == 1) && (substr($revcom_seq,($ref_pos-1),1) ne $reference_nts)) {
				print "\nSubstring: ".substr($revcom_seq,($ref_pos-1),1)."\nRef: $reference_nts\n";
				die "\n## WARNING: New nucleotide at site $ref_pos does not match reverse ".
					"complement. SNPGenie terminatred\n\n";
			}
			
			if($variant3) { # THERE ARE THREE VARIANTS -- ADD A FLAG!
				if($info_value =~ /NS=(\d+)/) { # We've got a VCF SUMMARIZING INDIVIDUALS
					print "\n### FILE TYPE NOT FULLY SUPPORTED###\n";
					my $num_samples;
					$num_samples = $1;
					
					my $variant_freq1;
					my $variant_freq2;
					my $variant_freq3;
					if($info_value =~ /AF=([\d\.]+),([\d\.]+),([\d\.]+)/) {
						$variant_freq1 = $1;
						$variant_freq2 = $2;
						$variant_freq3 = $3;
					}
					
					# COUNTS and PERCENTS
					my $ref_freq = (1 - ($variant_freq1 + $variant_freq2 + $variant_freq3));
					my $ref_count = ($ref_freq * $num_samples);
					
					my $variant_count1 = ($variant_freq1 * $num_samples);
					my $variant_count2 = ($variant_freq2 * $num_samples);
					my $variant_count3 = ($variant_freq3 * $num_samples);
					
					my $variant_pct1 = (100 * $variant_freq1);
					my $variant_pct2 = (100 * $variant_freq2);
					my $variant_pct3 = (100 * $variant_freq3);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							# PRINT 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\tCDS: $product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\tCDS: $product_entry\n";
								
								print TEMP_FILE_VCF "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$num_samples\t$variant_pct3\tCDS: $product_entry\n";
								
								print TEMP_FILE_VCF "$this_line3";
							}
						}
					} else {
						# PRINT 3 LINES TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\tCDS: $product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
						
						if($variant_freq2 > 0) {
							my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\tCDS: $product_entry\n";
							
							print TEMP_FILE_VCF "$this_line2";
						}
						
						if($variant_freq3 > 0) {
							my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant3\t$variant_count3\t$num_samples\t$variant_pct3\tCDS: $product_entry\n";
							
							print TEMP_FILE_VCF "$this_line3";
						}
					}
					
				} elsif($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) { # We've got a VCF of POOL
					# Warn
					print "\n## File $curr_snp_report_name site $ref_pos has three variants ".
					"in a pooled (DP4-tag) VCF file.\n## Variant frequencies have been approximated\n";
					
					# These are high-quality reads, so may be less that the actual coverage
					my $fwd_ref_reads = $1;
					my $rev_ref_reads = $2;
					my $fwd_alt_reads = $3;
					my $rev_alt_reads = $4;
					
					# COUNTS and FREQS
					my $ref_count = ($fwd_ref_reads + $rev_ref_reads);
					my $alt_count = ($fwd_alt_reads + $rev_alt_reads);
					my $coverage = ($ref_count + $alt_count);
					
					my $variant_count1 = ($alt_count / 3);
					my $variant_count2 = ($alt_count / 3);
					my $variant_count3 = ($alt_count / 3);
					
					my $variant_freq1 = ($variant_count1 / $coverage);
					my $variant_freq2 = ($variant_count2 / $coverage);
					my $variant_freq3 = ($variant_count3 / $coverage);
					
					my $variant_pct1 = (100 * $variant_freq1);
					my $variant_pct2 = (100 * $variant_freq2);
					my $variant_pct3 = (100 * $variant_freq3);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line3";
							}
						}
					} else {
						# PRINT 3 LINES TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
						
						if($variant_freq2 > 0) {
							my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line2";
						}
						
						if($variant_freq3 > 0) {
							my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line3";
						}
					}
					
				} elsif($format_value =~ /AD/) { # We've got a VCF of POOL; we need AD and DP
					# Find out how many ":" appear before AD, if any
					my $prior_to_AD = $`;
					my @colons_prior_to_AD = $prior_to_AD =~ /\:/g; 
					my $colon_count_before_AD = @colons_prior_to_AD;
					
					my $prior_to_DP;
					my @colons_prior_to_DP;
					my $colon_count_before_DP;
						
					if($format_value =~ /DP/) {
						$prior_to_DP = $`;
						@colons_prior_to_DP = $prior_to_DP =~ /\:/g;
						$colon_count_before_DP = @colons_prior_to_DP;
					} else {
						
						die "\n\n## WARNING: $curr_snp_report_name contains AD but not DP data. SNPGenie terminated\n\n";	 
					}
					
					# EXTRACT the VALUE of AD
					my @sample_value_arr = split(/\:/,$sample_value);
					my $AD1;
					my $AD2;
					my $AD3;
					my $AD4;
					
					if($sample_value_arr[$colon_count_before_AD] =~ /(\d+)\,(\d+)\,(\d+)\,(\d+)/) { # REF and 2 ALTS
						$AD1 = $1;
						$AD2 = $2;
						$AD3 = $3;
						$AD4 = $4;
					}
					
					# EXTRACT the VALUE of DP (coverage)
					my $DP = $sample_value_arr[$colon_count_before_DP];
			
					# COUNTS and FREQS
					my $ref_count = $AD1;
					my $variant_count1 = $AD2;
					my $variant_count2 = $AD3;
					my $variant_count3 = $AD4;
					my $coverage = ($AD1+$AD2+$AD3+$AD4);
					
					# Warn if the total reads don't equal the coverage
					if ($DP != $coverage) {
						print "\n## WARNING: In $curr_snp_report_name site $ref_pos".
								",\n## the reads total should ".
								"equal the coverage ($DP) but is instead: $coverage.".
								"\n## The reads total has been used. Please verify your data.\n";
					}
							
					my $variant_freq1 = ($variant_count1 / $coverage);
					my $variant_freq2 = ($variant_count2 / $coverage);
					my $variant_freq3 = ($variant_count3 / $coverage);
					
					my $variant_pct1 = (100 * $variant_freq1);
					my $variant_pct2 = (100 * $variant_freq2);
					my $variant_pct3 = (100 * $variant_freq3);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line3";
							}
						}
					} else {
						# PRINT 3 LINES TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
						
						if($variant_freq2 > 0) {
							my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line2";
						}
						
						if($variant_freq3 > 0) {
							my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line3";
						}
					}
				}
				
			} elsif($variant2) { # THERE ARE TWO VARIANTS -- ADD A FLAG!
				if($info_value =~ /NS=(\d+)/) { # We've got a VCF summarizing INDIVIDUALS
					print "\n### FILE TYPE NOT FULLY SUPPORTED###\n";
					my $num_samples;
					$num_samples = $1;
					
					my $variant_freq1;
					my $variant_freq2;
					if($info_value =~ /AF=([\d\.]+),([\d\.]+)/) {
						$variant_freq1 = $1;
						$variant_freq2 = $2;
					}
					
					# COUNTS and PERCENTS
					my $ref_freq = (1 - ($variant_freq1 + $variant_freq2));
					my $ref_count = ($ref_freq * $num_samples);
					
					my $variant_count1 = ($variant_freq1 * $num_samples);
					my $variant_count2 = ($variant_freq2 * $num_samples);
					
					my $variant_pct1 = (100 * $variant_freq1);
					my $variant_pct2 = (100 * $variant_freq2);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\tCDS: $product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\tCDS: $product_entry\n";
								
								print TEMP_FILE_VCF "$this_line2";
							}
						}
					} else {
						# PRINT 2 LINES TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\tCDS: $product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
						
						if($variant_freq2 > 0) {
							my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\tCDS: $product_entry\n";
							
							print TEMP_FILE_VCF "$this_line2";
						}
					}
					
				} elsif($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) { # We've got a VCF of POOL
					# Warn 
					print "\n## File $curr_snp_report_name site $ref_pos has two variants ".
					"in a pooled (DP4-tag) VCF file.\n## Variant frequencies have been approximated\n";
					
					# These are high-quality reads, so may be less that the actual coverage
					my $fwd_ref_reads = $1;
					my $rev_ref_reads = $2;
					my $fwd_alt_reads = $3;
					my $rev_alt_reads = $4;
					
					# COUNTS and FREQS
					my $ref_count = ($fwd_ref_reads + $rev_ref_reads);
					my $alt_count = ($fwd_alt_reads + $rev_alt_reads);
					my $coverage = ($ref_count + $alt_count);
					
					my $variant_count1 = ($alt_count / 2);
					my $variant_count2 = ($alt_count / 2);
					
					my $variant_freq1 = ($variant_count1 / $coverage);
					my $variant_freq2 = ($variant_count2 / $coverage);
					
					my $variant_pct1 = (100 * $variant_freq1);
					my $variant_pct2 = (100 * $variant_freq2);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line2";
							}
						}
					} else {
						# PRINT 2 LINES TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
							
						if($variant_freq2 > 0) {
							my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line2";
						}
					}
					
				} elsif($format_value =~ /AD/) { # We've got a VCF of POOL; we need AD and DP
					# Find out how many ":" appear before AD, if any
					my $prior_to_AD = $`;
					my @colons_prior_to_AD = $prior_to_AD =~ /\:/g; 
					my $colon_count_before_AD = @colons_prior_to_AD;
					#print "\n\nColon count before AD: $colon_count_before_AD\n";
					
					my $prior_to_DP;
					my @colons_prior_to_DP;
					my $colon_count_before_DP;
						
					if($format_value =~ /DP/) {
						$prior_to_DP = $`;
						@colons_prior_to_DP = $prior_to_DP =~ /\:/g;
						$colon_count_before_DP = @colons_prior_to_DP;
						#print "\n\nColon count before DP: $colon_count_before_DP\n";
					} else {
						
						die "\n\n## WARNING: $curr_snp_report_name contains AD but not DP data. SNPGenie terminated\n\n";	 
					}
					
					# GENERATE a regex for what comes before the VALUE of AD
					#my $prior_to_AD_regex;
					#for (my $i=0; $i<$colon_count_before_AD; $i++) {
					#	$prior_to_AD_regex .= "[\d\w\/\,]+\:";
					#}
					
					# EXTRACT the VALUE of AD
					my @sample_value_arr = split(/\:/,$sample_value);
					my $AD1;
					my $AD2;
					my $AD3;
					
					if($sample_value_arr[$colon_count_before_AD] =~ /(\d+)\,(\d+)\,(\d+)/) { # REF and 2 ALTS
						$AD1 = $1;
						$AD2 = $2;
						$AD3 = $3;
					}
					
					# EXTRACT the VALUE of DP (coverage)
					my $DP = $sample_value_arr[$colon_count_before_DP];
			
					# COUNTS and FREQS
					my $ref_count = $AD1;
					my $variant_count1 = $AD2;
					my $variant_count2 = $AD3;
					my $coverage = ($AD1+$AD2+$AD3);
					
					# Warn if the total reads don't equal the coverage
					if ($DP != $coverage) {
						print "\n## WARNING: In $curr_snp_report_name site $ref_pos".
								",\n## the reads total should ".
								"equal the coverage ($DP) but is instead: $coverage.".
								"\n## The reads total has been used. Please verify your data.\n";
					}
							
					my $variant_freq1 = ($variant_count1 / $coverage);
					my $variant_freq2 = ($variant_count2 / $coverage);
					
					my $variant_pct1 = (100 * $variant_freq1);
					my $variant_pct2 = (100 * $variant_freq2);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line2";
							}
						}
					} else {
						# PRINT 2 LINES TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
						
						if($variant_freq2 > 0) {
							my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line2";
						}
					}
				}
				
			} elsif($variant1) { # THERE IS ONE VARIANT -- no flag needed
				if($info_value =~ /NS=(\d+)/) { # We've got a VCF summarizing INDIVIDUALS
					print "\n### FILE TYPE NOT FULLY SUPPORTED###\n";
					my $num_samples;
					$num_samples = $1;
					
					my $variant_freq1;
					if($info_value =~ /AF=([\d\.]+)/) {
						$variant_freq1 = $1;
					}
					
					# COUNTS and PERCENTS
					my $ref_freq = (1 - $variant_freq1);
					my $ref_count = ($ref_freq * $num_samples);
					
					my $variant_count1 = ($variant_freq1 * $num_samples);
					
					my $variant_pct1 = (100 * $variant_freq1);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 1 LINE TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\tCDS: $product_entry\n";
						
								print TEMP_FILE_VCF "$this_line1";
							}
						}
					} else {
						# PRINT 1 LINE TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\tCDS: $product_entry\n";
						
							print TEMP_FILE_VCF "$this_line1";
						}
					}
					
				} elsif($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) { # We've got a VCF of POOL
					# These are high-quality reads, so may be less that the actual coverage
					my $fwd_ref_reads = $1;
					my $rev_ref_reads = $2;
					my $fwd_alt_reads = $3;
					my $rev_alt_reads = $4;
					
					# COUNTS and FREQS
					my $ref_count = ($fwd_ref_reads + $rev_ref_reads);
					my $alt_count = ($fwd_alt_reads + $rev_alt_reads);
					my $coverage = ($ref_count + $alt_count);
					
					my $variant_count1 = $alt_count;
					
					my $variant_freq1 = ($variant_count1 / $coverage);
					
					my $variant_pct1 = (100 * $variant_freq1);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 1 LINE TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
								print TEMP_FILE_VCF "$this_line1";
							}
						}
					} else {
						# PRINT 1 LINE TO FILE; $product_entry remains BLANK
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
					}
				} elsif($format_value =~ /AD/) { # We've got a VCF of POOL; we need AD and DP
					# Find out how many ":" appear before AD, if any
					my $prior_to_AD = $`;
					my @colons_prior_to_AD = $prior_to_AD =~ /\:/g; 
					my $colon_count_before_AD = @colons_prior_to_AD;
					#print "\n\nColon count before AD: $colon_count_before_AD\n";
					
					my $prior_to_DP;
					my @colons_prior_to_DP;
					my $colon_count_before_DP;
						
					if($format_value =~ /DP/) {
						$prior_to_DP = $`;
						@colons_prior_to_DP = $prior_to_DP =~ /\:/g;
						$colon_count_before_DP = @colons_prior_to_DP;
						#print "\n\nColon count before DP: $colon_count_before_DP\n";
					} else {
						
						die "\n\n## WARNING: $curr_snp_report_name contains AD but not DP data. SNPGenie terminated\n\n";	 
					}
					
					# EXTRACT the VALUE of AD
					my @sample_value_arr = split(/\:/,$sample_value);
					my $AD1;
					my $AD2;
					
					if($sample_value_arr[$colon_count_before_AD] =~ /(\d+)\,(\d+)/) { # REF and 1 ALT
						$AD1 = $1;
						$AD2 = $2;
					}
					
					# EXTRACT the VALUE of DP (coverage)
					my $DP = $sample_value_arr[$colon_count_before_DP];
			
					# COUNTS and FREQS
					my $ref_count = $AD1;
					my $variant_count1 = $AD2;
					my $coverage = ($AD1+$AD2);
					
					# Warn if the total reads don't equal the coverage
					if ($DP != $coverage) {
						print "\n## WARNING: In $curr_snp_report_name site $ref_pos".
								",\n## the reads total should ".
								"equal the coverage ($DP) but is instead: $coverage.".
								"\n## The reads total has been used. Please verify your data.\n";
					}
							
					my $variant_freq1 = ($variant_count1 / $coverage);
					
					my $variant_pct1 = (100 * $variant_freq1);
					
					$product_entry = '';
					if(@this_site_products) {
						foreach my $product (@this_site_products) {
							$product_entry = 'CDS: ' . $product;
							
							# PRINT 1 LINE TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE_VCF "$this_line1";
							}
						}
					} else {
						# PRINT 1 LINE TO FILE
						if($variant_freq1 > 0) {
							my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
								"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
							
							print TEMP_FILE_VCF "$this_line1";
						}
					}
				} # End of AD format
			} # End of 1 variant
		}
	}
}
close ORIGINAL_SNP_REPORT;
close TEMP_FILE_VCF;

open(VCF_TEMP, "$new_vcf_file_name\_TEMP");
my @new_lines = <VCF_TEMP>;
close VCF_TEMP;
unlink "$new_vcf_file_name\_TEMP";

open(OUT_FILE_VCF,">>$new_vcf_file_name");
my $clc_format_header = "File\tReference Position\tType\tReference\tAllele\tCount\t".
	"Coverage\tFrequency\tOverlapping annotations\n";
print OUT_FILE_VCF "$clc_format_header";
while (my $line = pop @new_lines) {
	if($line =~ /\n$/) {
		print OUT_FILE_VCF $line;
	} else {
		print OUT_FILE_VCF "$line\n";
	}
}
close OUT_FILE_VCF;

print "\n## REVERSE COMPLEMENT FILES have been written to:\n".
	"## VCF: $new_vcf_file_name\n".
	"## GTF: $new_gtf_file_name\n".
	"## FASTA: $new_fasta_file_name\n\n";

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
sub generate_reverse_complement_fasta {
	my ($filename) = @_;
	
	# Read in the sequence from the file
	my $seq = '';
	my $header = '';
	my $num_headers = 0;
	
	open(IN, "$filename") or die "\nCould not open file $filename\n";
	while(<IN>) {
		if(/>/) {
			if($num_headers == 0) {
				$header = $_;
				$num_headers ++;
			} else {
				die "\nThere is more than one FASTA header in the file $filename\n";
			}
		} else {
			chomp;
			$seq .= $_;
		}
	}
	
	chomp($header);
	my $rev_seq = reverse($seq);
	my $rev_compl = $rev_seq;
	$rev_compl =~ tr/ACGT/TGCA/;
	
	my $rev_compl_filename;
	if($fasta_file_nm =~/\.fasta/) { 
		$rev_compl_filename = $` . "_revcom.fasta";
	} elsif($fasta_file_nm =~/\.fa/) { 
		$rev_compl_filename = $` . "_revcom.fa";
	} else {
		die "\nFirst argument must be a .fa or .fasta file\n\n";
		#$rev_compl_filename = "fasta_revcom.fasta";
	}
	
	if(-e $rev_compl_filename) {
		die "\n## $new_fasta_file_name already exists in this directory; delete before proceeding\n\n";
	}
	
	my $seq_length = length($rev_compl);
	
	open(OUT_FILE_FASTA, ">>$rev_compl_filename");
	
	my $new_header = $header . " REVERSE COMPLEMENT\n";
	print OUT_FILE_FASTA $new_header;
	
	for(my $i=0; $i<($seq_length); $i+=60) {
		if($i>($seq_length - 60)) {
			my $line = substr($rev_compl, $i);
			print OUT_FILE_FASTA "$line\n";
			last;
		} else {
			my $line = substr($rev_compl, $i, 60);
			print OUT_FILE_FASTA "$line\n";
		}
	}
	
	close OUT_FILE_FASTA;
	
	return $rev_compl_filename;
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
	
	my %gtf_line_by_starting_position;
	
	open(GTF_FILE, "$filename") or die "\nCould not open the GTF file $filename - $!\n\n";
	while(<GTF_FILE>) {
		if($_ =~ /CDS\t(\d+)\t(\d+)/) {
			my $old_start = $1; # Where the gene itself actually STOPS
			my $old_stop = $2; # Where the gene itself actually STARTS
			
			# Find the coordinates from the revcom point of view
			my $this_start = $seq_length - $old_stop + 1;
			my $this_stop = $seq_length - $old_start + 1;
			
			if(exists $gtf_line_by_starting_position{$this_start}) {
				die "\n\nTwo products have the same starting position, causing an error.\n".
					"Please contact script author for a revision.\n\n";
			} else {
				$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\+/CDS\t$this_start\t$this_stop\t\.\t\-/;
				$_ =~ s/CDS\t$old_start\t$old_stop\t\.\t\-/CDS\t$this_start\t$this_stop\t\.\t\+/;
				$gtf_line_by_starting_position{$this_start} = $_;
			}
		} 
	}
	close GTF_FILE;
	
	# Give NUM LINES CONVERTED tag/warning
	my @stored_pos_sorted = sort {$a <=> $b} (keys %gtf_line_by_starting_position);
	my $num_positions = scalar(@stored_pos_sorted);
	print "\n## A total of $num_positions products were processed. Please verify that ".
		"this is the correct number of products.\n";
	
	open(OUT_FILE_GTF, ">>$new_gtf_file_name");
	foreach my $curr_pos (@stored_pos_sorted) {
		print OUT_FILE_GTF $gtf_line_by_starting_position{$curr_pos};
	}
	close OUT_FILE_GTF;
	
	return $new_gtf_file_name;
}