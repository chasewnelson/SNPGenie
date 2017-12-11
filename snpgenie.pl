#! /usr/bin/perl

# PROGRAM: SNPGenie is a Perl program that calculates dN/dS, piN/piS, and gene diversity
# from NGS SNP Reports generated from pooled DNA samples.

# Copyright (C) 2015, 2016, 2017 Chase W. Nelson

#########################################################################################
## LICENSE
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################################

# DATE CREATED: April 10, 2015

# AUTHOR: Chase W. Nelson

# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com

# AFFILIATION1: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA
# AFFILIATION2: Special Volunteer, Division of Cancer Epidemiology & Genetics, National
#     Cancer Institute, National Institutes of Health, Rockville, MD 20850, USA
# AFFILIATION3: BigPlant Consortium, Center for Genomics and Systems Biology, New York 
#     University, New York, NY 10003, USA

# CITATION1: SNPGenie, https://github.com/chasewnelson/snpgenie
# CITATION2: Nelson CW, Moncla LH, Hughes AL (2015) SNPGenie: estimating evolutionary 
#	parameters to detect natural selection using pooled next-generation sequencing data. 
#	Bioinformatics 31(22):3709-11, doi: 10.1093/bioinformatics/btv449.

use strict;
#use warnings;
use IO::Handle;
use Data::Dumper;
use File::Temp qw(tempfile);
use Getopt::Long;
use List::Util qw(max sum);

my $time1 = time;
my $local_time1 = localtime;

# Die if no arguments, explaining what they should be
#unless (@ARGV) {
#	die "\nThis program accepts arguments:\n[1] Blah;\n[2] Blah2;\n[3] Blah3.\n\n";
#}

# GET THE TIME

# Initialize variables
my $minfreq; # the min. allele freq. to be considered (prop., not percentage)
my $snpreport;
my $vcfformat; ##SAMVCF
my $fastafile;
my $gtffile;
my $sepfiles;
my $slidingwindow;
my $ratiomode;
my $sitebasedmode; # not supported or recommended; site (not codon) contexts
my $complementmode;
#my $prop_diff_threshold = 0.01;

# Flag variables
my $clc_mode = 0;
my $geneious_mode = 0;
my $vcf_mode = 0;
my $generate_random_fastas = 0;
#my $progress_period_count = 0;
my @snp_report_file_names_arr;
my $the_fasta_file = '';
my @fasta_file_names_arr;
my $fasta_arr_size;
my $cds_file;
my $multi_fasta_mode;

my $param_file_contents = "SNPGenie version 1.2 parameter log.\n\n";

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions(	"minfreq:f" => \$minfreq, # optional floating point parameter
			"snpreport:s" => \$snpreport, # optional string parameter
			"vcfformat:i" => \$vcfformat, # optional integer parameter ##SAMVCF
			"fastafile:s" => \$fastafile, # optional string parameter
			"multi_fasta_mode" => \$multi_fasta_mode, # optional Boolean parameter
			"gtffile:s" => \$gtffile, # optional string parameter
			"sepfiles" => \$sepfiles, # optional Boolean; set to false if not given
			"slidingwindow:i" => \$slidingwindow, # optional integer parameter
			"ratiomode" => \$ratiomode, # optional Boolean; set to false if not given
			"sitebasedmode" => \$sitebasedmode) # optional Boolean; set to false if not given
#			"complementmode" => \$complementmode) # optional Boolean; set to false if not given
#			"prop_diff_threshold:f" => \$prop_diff_threshold)
			
			or die "\n## WARNING: Error in command line arguments. SNPGenie terminated.\n\n"; 

# N.B.: When an argument, e.g., slidingwindow, is called only as a flag, its value is 0
# When it is not called at all, it is null

# Create a directory for the results
if (-d "SNPGenie_Results") { # Can also use "./SNPGenie_Results"; use "-e" to check file
	die "\n\n## WARNING:\n## The directory SNPGenie_Results already exists in the ".
			"working directory.\n## Please rename or move this directory so a new one ".
			"can be created.\n\n";
} else {
	mkdir('SNPGenie_Results');
	
	## Set OPTIONS given the user's INPUT and RECORD PARAMETERS	
	if(! $minfreq) {
		$minfreq = 0;
		$param_file_contents .= "MINIMUM ALLELE FREQUENCY: Default used; all SNPs included\n";
	} elsif(($minfreq >= 1) || ($minfreq < 0)) {
		die "\n## WARNING: The --minfreq option must be a decimal between 0 and 1\n".
			"## SNPGenie terminated.\n\n";
	} else {
		$param_file_contents .= "MINIMUM ALLELE FREQUENCY: $minfreq\n";
	}
	
	if(! $sepfiles) {
		$sepfiles = 0; # default behavior: no separate codon files for each SNP Report
		$param_file_contents .= "SEPARATE FILES OUTPUT: No\n";
	} else {
		$sepfiles = 1;
		$param_file_contents .= "SEPARATE FILES OUTPUT: Yes\n";
	}
	
	if($slidingwindow == 0) { # Called as a flag, but given no value
		$slidingwindow = 9; # default behavior: nonamer peptide
		$param_file_contents .= "SLIDING WINDOW LENGTH: Default used; 9 codons\n";
	} elsif($slidingwindow > 0) {
		$param_file_contents .= "SLIDING WINDOW LENGTH: $slidingwindow\n";
	} else {
		$param_file_contents .= "SLIDING WINDOW LENGTH: None\n";
	}
	
	if(! $ratiomode) {
		$ratiomode = 0; # default behavior: no separate codon files for each SNP Report
		$param_file_contents .= "RATIO MODE: Default used; No\n";
	} else {
		$ratiomode = 1;
		$param_file_contents .= "RATIO MODE: Yes\n";
	}
	
	if(! $sitebasedmode) {
		$sitebasedmode = 0; # default behavior: no separate codon files for each SNP Report
		$param_file_contents .= "SITE-BASED MODE: Default used; No\n";
	} else {
		$sitebasedmode = 1;
		$param_file_contents .= "SITE-BASED MODE: Yes\n";
	}
	
	if(! $multi_fasta_mode) {
		$multi_fasta_mode = 0; # default behavior: no separate codon files for each SNP Report
		$param_file_contents .= "MULTI-FASTA MODE: Default used; No\n";
	} else {
		$multi_fasta_mode = 1;
		$param_file_contents .= "MULTI-FASTA MODE: Yes\n";
	}
	
	# Get SNP Report name(s)
	if(! $snpreport) {
		@snp_report_file_names_arr = &get_txt_file_names;
		my @snp_report_file_names_ADD_arr = &get_csv_file_names;
		my @snp_report_file_names_ADD_VCF_arr = &get_vcf_file_names;
		push(@snp_report_file_names_arr,@snp_report_file_names_ADD_arr);
		push(@snp_report_file_names_arr,@snp_report_file_names_ADD_VCF_arr);
		$param_file_contents .= "SNP REPORTS: Default auto-detected file(s): @snp_report_file_names_arr\n";
	} else {
		@snp_report_file_names_arr = ($snpreport);
		$param_file_contents .= "SNP REPORTS: User submitted file: $snpreport\n";
	}
	#print "\n@snp_report_file_names_arr\n\n";
	
	if (scalar (@snp_report_file_names_arr) == 0) {
		die "\n\n## WARNING: There are no SNP Reports. SNPGenie terminated.\n\n";
	} else {
		foreach(@snp_report_file_names_arr) {
			if($_ =~ /\.vcf/ && ! $vcfformat) {
				rmdir('SNPGenie_Results');
				die "\n\n### WARNING: User must specify the specific SNP report format when using\n".
					"### VCF files. Use the --vcfformat option. SNPGENIE TERMINATED.\n\n";
			}
		}
	}
	
	if(! $fastafile) {
		@fasta_file_names_arr = &get_fasta_file_names;
		#print "\nWorking directory fasta files are: @fasta_file_names_arr\n\n";
		$param_file_contents .= "REFERENCE FASTA FILE: Default auto-detected file(s): @fasta_file_names_arr\n";
	} else {
		$fasta_file_names_arr[0] = $fastafile;
		$param_file_contents .= "REFERENCE FASTA FILE: User submitted file: $fastafile\n";
	}
	
	$fasta_arr_size = scalar(@fasta_file_names_arr);
	#print "\nThe size of the fasta array is $fasta_arr_size\n";
	
	
	
	if($fasta_arr_size > 1) {
		
		if(! $multi_fasta_mode) {
			die "\n\n## WARNING: There are multiple FASTA (.fa or .fasta) files in the working directory.\n".
			"## There must be only one reference genome in single-FASTA mode. SNPGenie terminated.\n\n";
		}
		
	} elsif($fasta_arr_size == 1) {
		$the_fasta_file = $fasta_file_names_arr[0];
		
		# Go through FASTA to make sure there's only one sequence
		my $seen_seq = 0;
		open (INFILE, $the_fasta_file);
		while (<INFILE>) {
			if (/>/) {
				if($seen_seq == 0) {
					$seen_seq = 1;
				} else {
					die "\n\n### WARNING: There are multiple sequences in the FASTA file.\n".
						"### There must be only one reference genome. SNPGenie terminated.\n\n";
				}
			}
		}
		close INFILE;
		
		
	} else {
		die "\n\n## WARNING: There are no FASTA (.fa or .fasta) files in the working directory. ".
		"SNPGenie terminated.\n\n";
	}
	
	if(! $gtffile) { # default behavior
		$cds_file = &get_cds_file_name;
		$param_file_contents .= "GTF FILE: Default auto-detected file: $cds_file\n";
	} else {
		$cds_file = $gtffile;
		$param_file_contents .= "GTF FILE: User submitted file: $cds_file\n";
	}
	
	#print "\n@fasta_file_names_arr\n\n";
	
	STDOUT->autoflush(1);
	
	chdir('SNPGenie_Results');
	
	open(PARAM_FILE,">>SNPGenie\_parameters\.txt");
	print PARAM_FILE "$param_file_contents";
	close PARAM_FILE;
	
	### NUCLEOTIDE DIVERSITY FILE
	open(OUTFILE_NT_DIV,">>codon\_results\.txt");
	my $ntd_headers_to_print = "file\tproduct\tsite\tcodon\tnum_overlap_ORF_nts\t".
			"N_diffs\tS_diffs\t";
	if($sitebasedmode == 1) {
		$ntd_headers_to_print .= "N_diffs_site_based\tS_diffs_site_based\t";
	}
	$ntd_headers_to_print .= "N_sites\tS_sites\t";
	if($sitebasedmode == 1) {
		$ntd_headers_to_print .= "N_sites_site_based\tS_sites_site_based\t";
	}
	$ntd_headers_to_print .= "N_sites_ref\tS_sites_ref\t";
	#$ntd_headers_to_print .= "piN\tpiS\t";
	if($ratiomode == 1) {
		$ntd_headers_to_print .= "piN/piS\t";
	}
	$ntd_headers_to_print .= "N_diffs_vs_ref\tS_diffs_vs_ref\t".
		"gdiv\tN_gdiv\tS_gdiv\n";
	#$ntd_headers_to_print .= "mean_dN_vs_ref\tmean_dS_vs_ref\n"; # \tAverage_cov
	
	print OUTFILE_NT_DIV "$ntd_headers_to_print";
	close OUTFILE_NT_DIV;
	
	### GENE DIVERSITY FILE
	open(OUTFILE_GENE_DIV,">>site\_results\.txt");
	print OUTFILE_GENE_DIV "file\tproduct\tsite\tref_nt\tmaj_nt\t".
		"position_in_codon\t".
		"overlapping_ORFs\tcodon_start_site\tcodon\tpi\t".
		#"Polymorphic (Y=1; N=0)\t".
		"gdiv\t".
		"class_vs_ref\tclass\t".
		"coverage\t".
		"A\tC\tG\tT\n";
	close OUTFILE_GENE_DIV;
	
	### PRODUCT SUMMARY FILE
	open(PRODUCT_SUMMARY,">>product\_results\.txt");
	print PRODUCT_SUMMARY "file\tproduct\tN_diffs\tS_diffs\t".
		"N_diffs_vs_ref\tS_diffs_vs_ref\t".
		"N_sites\tS_sites\t".
		"piN\tpiS\tmean_dN_vs_ref\tmean_dS_vs_ref\t".
		"mean_gdiv_polymorphic\tmean_N_gdiv\tmean_S_gdiv\n";
	close PRODUCT_SUMMARY;
	
	### POPULATION SUMMARY NONCODING RESULTS
	open(POP_SUMMARY,">>population\_summary\.txt");
	print POP_SUMMARY "file\tsites\tsites_coding\tsites_noncoding\t".
		"pi\tpi_coding\tpi_noncoding\t".
		#"mean_nonsyn_diffs\tmean_syn_diffs\t".
		#"mean_nonsyn_diffs_vs_ref\tmean_syn_diffs_vs_ref\t".
		"N_sites\tS_sites\t".
		"piN\tpiS\tmean_dN_vs_ref\tmean_dS_vs_ref\t".
		"mean_gdiv_polymorphic\tmean_N_gdiv\tmean_S_gdiv\t".
		"mean_gdiv\t".
		"sites_polymorphic\t".
		"mean_gdiv_coding_poly\t".
		"sites_coding_poly\t".
		"mean_gdiv_noncoding_poly\t".
		"sites_noncoding_poly\n";
	close POP_SUMMARY;
	
	
	### A LOG file
	open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
	print ERROR_FILE "file\tproduct\tsite\t".
			"LOG\n";
	close ERROR_FILE;
	
	chdir('..');
}

# Hash for storing which product we've seen, just for error-reporting purposes
my %seen_product_early_stop_hash;
my %seen_no_start_hash;
my %seen_no_stop_hash;

# Executive error switches
my $exec_errors = 0;
my $warn_5nt = 0;
my $warn_frequencies = 0;
my $warn_file_type_not_supported = 0;
my $seen_sense_strand_products = 0;

my $SNP_report_counter = 0;

print "\n\n################################################################################".
	"\n##                                                                            ##".
	"\n##                             SNPGenie Initiated!                            ##".
	"\n##                                                                            ##".
	"\n################################################################################\n";

# Print LICENSE
print "\n  ###############################  LICENSE:  #################################\n";
print "  ##              SNPGenie Copyright (C) 2015 Chase W. Nelson               ##\n".
	"  ##            This program comes with ABSOLUTELY NO WARRANTY;             ##\n".
	"  ##     This is free software, and you are welcome to redistribute it      ##\n".
	"  ##               under certain conditions; see LICENSE.txt.               ##";
print "\n  ############################################################################\n";

# GET THE TIME
chdir('SNPGenie_Results');
open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
print ERROR_FILE "NA\tNA\tNA\t".
		"SNPGenie initiated at local time $local_time1\n";
close ERROR_FILE;
chdir('..');

if($minfreq > 0) {
	print "\nYour MIN. MINOR ALLELE FREQ. is $minfreq. All variants falling below this frequency will be ignored.\n";
} else {
	print "\nYou have not selected a MIN. MINOR ALLELE FREQ. All variants in the SNP report(s) will be included.\n";
}

# We will 
# [1] determine whether COMPLEMENT ENTRIES are to be considered, and 
# [2] if so, construct some way to DO EVERYTHING BELOW, but do it for the
# rev. complement SNP reports against the rev. complement FASTA with respect to the 
# "-" strand records in the GTF file.

# Complement mode?
$complementmode = &determine_complement_mode($cds_file);

# Announce and initialize REVERSE COMPLEMENT MODE
my %hh_compl_position_info; # REGARDLESS OF SNP REPORT. saved with respect to the + strand
#my @curr_compl_products;
my @curr_compl_products_ordered_by_start;

if($complementmode) {
	print "\nThere are antisense ('-') strand records in the GTF file. COMPLEMENT MODE activated...\n";
	
	chdir('SNPGenie_Results');
	open(PARAM_FILE,">>SNPGenie\_parameters\.txt");
	print PARAM_FILE "COMPLEMENT MODE: Yes\n";
	close PARAM_FILE;
	chdir('..');
	
	# Look through the GTF file for the - strand entries
	# translate the start and stop sites to + strand sites using $fasta_length
	#my $rev_complement_seq = &reverse_complement_from_fasta($fasta_to_open);
	#my $rev_compl_seq = &reverse_complement_from_fasta($the_fasta_file);
	#my $seq_length = length($rev_compl_seq);
	open(GTF_FILE_AGAIN, "$cds_file") or die "\nCould not open the GTF file $cds_file - $!\n\n";
	while(<GTF_FILE_AGAIN>) { # each record in the GTF file
		my $this_product;
		
		# Reverse complement - strand
		my $rev_compl_start; # Where the gene itself actually STOPS
		my $rev_compl_stop; # Where the gene itself actually STARTS

		if($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\-\t\d+\t\s*gene_id\s*\"gene\:([\w\s\.\-\:']+)\"/) { # Line is - strand
			$rev_compl_start = $1; # Where the gene itself actually STOPS
			$rev_compl_stop = $2; # Where the gene itself actually STARTS
			$this_product = $3;
		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\-\t\d+\t\s*gene_id\s*\"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
			$rev_compl_start = $1; # Where the gene itself actually STOPS
			$rev_compl_stop = $2; # Where the gene itself actually STARTS
			$this_product = $3;
		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\-\t\d+\t\s*gene_id\s*\"([\w\s\.\-\:']+)\"/) {
			$rev_compl_start = $1; # Where the gene itself actually STOPS
			$rev_compl_stop = $2; # Where the gene itself actually STARTS
			$this_product = $3;
		# NOW, IN CASE transcript_id comes first
		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\-\t\d+\ttranscript_id \"[\w\s\.\-\:']+\"\s*;\s*gene_id\s*\"gene\:([\w\s\.\-\:']+)\"/) {
			$rev_compl_start = $1; # Where the gene itself actually STOPS
			$rev_compl_stop = $2; # Where the gene itself actually STARTS
			$this_product = $3;
		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\-\t\d+\ttranscript_id \"[\w\s\.\-\:']+\"\s*;\s*gene_id\s*\"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
			$rev_compl_start = $1; # Where the gene itself actually STOPS
			$rev_compl_stop = $2; # Where the gene itself actually STARTS
			$this_product = $3;
		} elsif($_ =~ /CDS\t(\d+)\t(\d+)\t\.\t\-\t\d+\ttranscript_id \"[\w\s\.\-\:']+\"\s*;\s*gene_id\s*\"([\w\s\.\-\:']+)\"/) {
			$rev_compl_start = $1; # Where the gene itself actually STOPS
			$rev_compl_stop = $2; # Where the gene itself actually STARTS
			$this_product = $3;
		}
		
		# New segments approach for reverse complement - strand
		if($rev_compl_start) { # it's a revcom record
			#print "\nDo we ever get here 1? Yes, for $this_product\n";
			# so we have this record's product name, start, and stop
			my $curr_start_key = 'start_1';
			my $curr_stop_key = 'stop_1';
			my $segment_number = 1;
			
			if(! $hh_compl_position_info{$this_product}->{$curr_start_key}) { # no first segment yet
				#print "\nDo we ever get here 2? Yes, for $this_product\n";
				if($this_product ne '') {
					$hh_compl_position_info{$this_product}->{$curr_start_key} = $rev_compl_start;
					$hh_compl_position_info{$this_product}->{$curr_stop_key} = $rev_compl_stop;
				}
			} else { # already a start_1
				while($hh_compl_position_info{$this_product}->{$curr_start_key}) {
					#print "\nDo we ever get here 3? Yes, for $this_product\n";
					my $curr_segment_number = $segment_number;
					$segment_number++;
					$curr_start_key =~ s/\_$curr_segment_number/\_$segment_number/;
					$curr_stop_key =~ s/\_$curr_segment_number/\_$segment_number/;
				}
				
				if($this_product ne '') {
					$hh_compl_position_info{$this_product}->{$curr_start_key} = $rev_compl_start;
					$hh_compl_position_info{$this_product}->{$curr_stop_key} = $rev_compl_stop;
				}
			}
			
			# Store (initiate OR update) number of segments
			$hh_compl_position_info{$this_product}->{num_segments} = $segment_number;
		} 
	}
	close GTF_FILE_AGAIN;
	
	#@curr_compl_products = sort(keys %hh_compl_position_info);
	@curr_compl_products_ordered_by_start = sort { $hh_compl_position_info{$a}->{start_1} <=> $hh_compl_position_info{$b}->{start_1} } keys %hh_compl_position_info;
	
	#print "\nproduct\tstart\tstop\n";
	#foreach (@curr_compl_products_ordered_by_start) {
	#	print "$_\t" . $hh_compl_position_info{$_}->{start_1} . "\t" . $hh_compl_position_info{$_}->{stop_1} . "\n";
	#}
} else {
	#print "\nThere are NO - strand records in the GTF file. COMPLEMENT MODE NOT activated...\n";
	chdir('SNPGenie_Results');
	open(PARAM_FILE,">>SNPGenie\_parameters\.txt");
	print PARAM_FILE "COMPLEMENT MODE: No\n";
	close PARAM_FILE;
	chdir('..');
}

#foreach my $this_product(keys %hh_compl_position_info) {
#	foreach(sort keys $hh_compl_position_info{$this_product}) {
#		print "\n product $this_product key position $_\n";
#	}
#}

# Streamline attainment of product coordinates
# All products should be in GTF file, so don't need to get them from SNP reports as before
my @product_names_arr = &get_product_names_from_gtf($cds_file);
@product_names_arr = sort {$a <=> $b} @product_names_arr;

# DIE if no sense + strand products seen
if(! $seen_sense_strand_products) {
	chdir('SNPGenie_Results');
	open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
	print ERROR_FILE "$cds_file\tNA\tNA\t".
		"Does not contain any sense (+) strand products. SNPGenie terminated.\n";
	close ERROR_FILE;
	chdir('..');
	
	die "\n\n## WARNING: $cds_file does not contain any sense (+) strand products. SNPGenie terminated.\n\n";
}

# Determine all product coordinates from the start
my %product_coordinates_harr;
foreach my $product_name (@product_names_arr) {
	#print "\nproduct name is: $product_name\n";
	my @product_coord_arr = &get_product_coordinates($product_name);
	#print "\n$product_name product_coord arr: @product_coord_arr\n";
	
	# Save to a hash of arrays
	push(@{$product_coordinates_harr{$product_name}->{product_coord_arr}},@product_coord_arr);
	#print "\n$product_name product_coord arr in harr: @{$product_coordinates_harr{$product_name}->{product_coord_arr}} \n\n";
}

# Streamline bulding of sequence
# Using the fasta file, record the sequence in a variable

my $seq;
my @seq_by_index_arr;

if(! $multi_fasta_mode) {

	print "\nReading in FASTA file... ";
	open (INFILE, $the_fasta_file);
	while (<INFILE>) {
		unless (/>/) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			$seq .= $_;
		}
	}
	close INFILE;
	print "COMPLETED.\n";
	
	# In case it's lowercase
	#$seq =~ tr/acgt/ACGT/;
	$seq = uc($seq);
	$seq =~ tr/U/T/;
	
	# Record in an array by index (old %seq_by_pos_hash)
	print "\nIndexing sequence... ";
	
	for (my $i = 0; $i < length($seq); $i++) {
		$seq_by_index_arr[$i] = substr($seq,$i,1); # This is $position - 1
	}
	print "COMPLETED.\n";
	
}

#########################################################################################
# GENERATE TEMP VCF SNP REPORTS IFF FORMAT 4

my @temp_vcf4_file_names;

# If it's VCF file format 4, it's possible that there are multiple samples in the same
# VCF file. Unfortunately, SNPGenie is modular in the sense that I designed it to be run
# on ONE SNP REPORT. Thus, to make a quick fix and avoid rewriting the whole concept,
# we'll just make each sample into its own temporary SNP report.
if($vcfformat == 4) { # generate as many SNP reports as there are sample columns
	
	# Find out if there is more than one sample by counting columns after FORMAT
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
	my $header_line;
	open (ORIGINAL_SNP_REPORT, $snp_report_file_names_arr[0]); # always just one
	while (<ORIGINAL_SNP_REPORT>) {	
		chomp;
			
		if($_ =~ /^#(\w+)/) {
			$header_line = $_;
			if(!($_ =~/\t/)) {
				die "\n\n## WARNING:\n## The SNP Report $snp_report_file_names_arr[0] is ".
					"not TAB-delimited (\\t), or there is only one column.\n\n";
			}
			last;
		}
	}
	close ORIGINAL_SNP_REPORT;
	
	$header_line =~ s/^#//;
	#print "\nHEADER LINE IS: $header_line\n\n";
	my @header_arr = split("\t",$header_line,-1);
	#print "\nHEADER ARRAY IS: @header_arr\n\n";
	
	# Determine the INDEX of FORMAT
	my $FORMAT_index;
	for (my $vcf_i = 0; $vcf_i < scalar(@header_arr); $vcf_i++) {
		if($header_arr[$vcf_i] eq 'FORMAT') {
			$FORMAT_index = $vcf_i;
			last;
		}
	}
	
	# Count the number of samples
	my $num_samples = 0;
	my @sample_names;
	for(my $samp_index = $FORMAT_index+1; $samp_index < scalar(@header_arr); $samp_index++) {
		if($header_arr[$samp_index] =~ /[\w\d\_\.]+/) {
			$num_samples++;
			push(@sample_names, $&);
		}
	}
	
	#print "\n### VCF FORMAT 4 with $num_samples samples.\n\n";
	
	# Generate TEMP VCF SNP reports, ranging from index ($FORMAT_index+1) to ($FORMAT_index+$num_samples)
	for(my $sample_index = ($FORMAT_index+1); $sample_index <= ($FORMAT_index+$num_samples); $sample_index++) {
		my $new_temp_vcf4_file_name = 'temp_vcf4_' . $header_arr[$sample_index] . '.vcf';
		
		push(@temp_vcf4_file_names, $new_temp_vcf4_file_name);
		
		print "\nCreating $new_temp_vcf4_file_name\...\n\n";
		
		if(-f "$new_temp_vcf4_file_name") {
			warn "\n\n### WARNING: $new_temp_vcf4_file_name already present. Replacing...\n\n";
			unlink $new_temp_vcf4_file_name; # just in case
		}
		
		open(THIS_NEW_VCF,">>$new_temp_vcf4_file_name");
		
		# Find out if there is more than one sample by counting columns after FORMAT
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
		my $header_line;
		
		open (ORIGINAL_SNP_REPORT, $snp_report_file_names_arr[0]); # always just one
		my $seen_vcf4_header = 0;
		while (<ORIGINAL_SNP_REPORT>) {	
			chomp;
			
			if($_ =~ /^##(\w+)/) { # STILL METADATA; print
				print THIS_NEW_VCF "$_\n";
			} elsif($_ =~ /^#(\w+)/) { # header
				my $this_new_header_line = '#';
				
				for(my $vcf_i = 0; $vcf_i <= $FORMAT_index; $vcf_i++) {
					$this_new_header_line .= "$header_arr[$vcf_i]\t";
				}
				
				$this_new_header_line .= "$header_arr[$sample_index]";
				
				print THIS_NEW_VCF "$this_new_header_line\n";
				
				$seen_vcf4_header = 1;
				
			} elsif($seen_vcf4_header == 1) {
				my @this_line_arr = split(/\t/,$_,-1);
				
				my $this_new_variant_line;
				unless($this_line_arr[$sample_index] =~ /\.\:\./) { # ZERO COVERAGE
					for(my $vcf_i = 0; $vcf_i <= $FORMAT_index; $vcf_i++) {
						$this_new_variant_line .= "$this_line_arr[$vcf_i]\t";
					}
					
					$this_new_variant_line .= "$this_line_arr[$sample_index]";
					
					print THIS_NEW_VCF "$this_new_variant_line\n";
				}
			}
		}
		close ORIGINAL_SNP_REPORT;
		
		close THIS_NEW_VCF;
		
	}
	
	@snp_report_file_names_arr = @temp_vcf4_file_names; # THESE are the SNP REPORTS now
	
}

#########################################################################################
# FOR METAPOPULATION STORAGE -- provided all sites are with respect to the same reference
my %master_frequencies_hh; # $master_frequencies_hh -> {site_num} -> {@A_props_arr}/{@C_props_arr}/{@G_props_arr}/{@T_props_arr}

#########################################################################################
# PROCESS THE SNP REPORTS
foreach my $curr_snp_report_name (@snp_report_file_names_arr) {
	my $file_nm = $curr_snp_report_name;
	#my $curr_newline = &detect_newline_char($curr_snp_report_name);
	#$/ = $curr_newline;
	
	my %h_nc_results;
	my $seen_percent_warning = 0;
	
	print "\n\n###########################  CURRENTLY PROCESSING:   ###########################\n".
	"$file_nm... ";
	
	#print "\n\n$_\n\n";
	
	if($multi_fasta_mode) {
		my $did_we_detect = 0;
		
		my $snp_report_prefix;
		if($file_nm =~ /^([\w\.]+?)_/) {
			$snp_report_prefix = $1;
		}
		
		foreach my $potential_fasta (@fasta_file_names_arr) {
			my $potential_fasta_prefix;
			
			if($potential_fasta =~ /^$snp_report_prefix\_/) {
				$the_fasta_file = $potential_fasta;
				$did_we_detect = 1;
				last;
			}
		}
		
		unless($did_we_detect == 1) {
			die "\n\nWe did not detect a FASTA for $file_nm in multi-fasta mode\n\n";
		}
		
		print "\nReading in FASTA file $the_fasta_file\... ";
		
		$seq = '';
		
		open (INFILE, $the_fasta_file);
		while (<INFILE>) {
			unless (/>/) {
				chomp;
				# CHOMP for 3 operating systems
				if($_ =~ /\r\n$/) {
					$_ =~ s/\r\n//;
				} elsif($_ =~ /\r$/) {
					$_ =~ s/\r//;
				} elsif($_ =~ /\n$/) {
					$_ =~ s/\n//;
				}
				
				$seq .= $_;
			}
		}
		close INFILE;
		print "COMPLETED.\n";
		
		# In case it's lowercase
		#$seq =~ tr/acgt/ACGT/;
		$seq = uc($seq);
		$seq =~ tr/U/T/;
		
		# Record in an array by index (old %seq_by_pos_hash)
		print "\nIndexing sequence... ";
		
		# THIS WILL OVERWRITE
		for (my $i = 0; $i < length($seq); $i++) {
			$seq_by_index_arr[$i] = substr($seq,$i,1); # This is $position - 1
		}
		print "COMPLETED.\n";
		
	} # end recording of individual FASTA in multi-FASTA mode
	
	
	# Generate new file name prefix
	my $new_file_prefix;
	if($curr_snp_report_name =~/\.txt/) { 
		$new_file_prefix = $`;
	} elsif($curr_snp_report_name =~/\.csv/) { 
		$new_file_prefix = $`; 
	} else {
		$new_file_prefix = "inFile";
	}
	
	#print "\nPrefix for $curr_snp_report_name: $new_file_prefix\n\n";
	
	# At first, created the tempfile before getting its headers; we've reversed the 
	# order so we can determine the format of the SNP Report (e.g., Geneious or CLC)
	# before creating the tempfile. Should not create problems.
	
	my @header_names_arr = &get_header_names($curr_snp_report_name,$curr_snp_report_name);
	#print "@header_names_arr";
	
	# CHECK FOR NON-CLC SNP REPORT, and reset the MODE if so
	foreach (@header_names_arr) {
		if($_ eq 'Minimum' || $_ eq 'Polymorphism Type' || $_ eq 'Change') {
			# Switch to GENEIOUS mode
			$geneious_mode = 1;
			$clc_mode = 0;
			$vcf_mode = 0;
			#print "\n\nWe are in GENEIOUS MODE!\n\n";
		} elsif($_ eq 'Reference Position' || $_ eq 'Overlapping annotations') {
			# We remain in CLC mode
			$geneious_mode = 0;
			$clc_mode = 1;
			$vcf_mode = 0;
			#print "\n\nWe are in CLC MODE!\n\n";
		} elsif($curr_snp_report_name =~ /\.vcf/) {
			# Switch to VCF mode
			$geneious_mode = 0;
			$clc_mode = 0;
			$vcf_mode = 1;
		}
	}
	
	##SAMVCF: changed to simple Boolean
	if($geneious_mode && ! $clc_mode && ! $vcf_mode) {
		print "GENEIOUS format detected\n";
	} elsif(! $geneious_mode && $clc_mode && ! $vcf_mode) {
		print "CLC GENOMICS WORKBENCH format detected\n";
	} elsif(! $geneious_mode && ! $clc_mode && $vcf_mode) {
		print "VCF format detected\n";
	} else {
		die "\n## WARNING: Conflicting SNP Report formats detected. Please contact author. ".
			"## SNPGenie TERMINATED.\n\n";
	}
	
	# Create tempfiles
	my $temp_snp_report_TEMPLATE = "temp_snp_report_XXXX";
	# Remember, the following automatically OPENS the file, too.
	my ($TEMP_SNP_REPORT_HANDLE,$temp_snp_report_name) = 
		tempfile($temp_snp_report_TEMPLATE, SUFFIX => ".txt", UNLINK => 1);
	#print "\n\n\nTEMP FILE: $temp_snp_report_name\n\n\n";
	
	##SAMVCF: changed to simple Boolean
	# POPULATE the temporary file based on FORMAT MODE
	if($clc_mode) {
		&populate_tempfile_clc($curr_snp_report_name,$temp_snp_report_name);
		$curr_snp_report_name = $temp_snp_report_name;
	} elsif($geneious_mode) {
		# (1) snpgenie_prep_geneious;
		# (2) snpgnie_geneious_to_clc;
		&populate_tempfile_geneious($curr_snp_report_name,$temp_snp_report_name);
		$curr_snp_report_name = $temp_snp_report_name;
	} elsif($vcf_mode) {
		if(! $vcfformat) {
			die "\n## WARNING: User must specify VCF format, e.g., --vcfformat=4. ".
				"## SNPGenie TERMINATED.\n\n";
		}
		
		&populate_tempfile_vcf($curr_snp_report_name,$temp_snp_report_name,$cds_file);
		$curr_snp_report_name = $temp_snp_report_name;
	}
	
	# Includes the automatic deletion of the tempfile afterwards. 
	my @new_header_names_arr = &get_header_names($curr_snp_report_name,$file_nm);
	@header_names_arr = @new_header_names_arr;
	#print "\n\nHEADER:\n@header_names_arr\n\n";
	
	my $index_ref_pos;
	my $index_type;
	my $index_ref;
	my $index_allele;
	my $index_count;
	my $index_cov;
	my $index_freq;
	my $index_over_annot;
	#my $index_cod_reg_chg;
	#my $index_ami_aci_chg;
	
	my $seen_index_ref_pos = 0;
	my $seen_index_type = 0;
	my $seen_index_ref = 0;
	my $seen_index_allele = 0;
	my $seen_index_count = 0;
	my $seen_index_cov = 0;
	my $seen_index_freq = 0;
	my $seen_index_over_annot = 0;
	#my $seen_index_cod_reg_chg = 0;
	#my $seen_index_ami_aci_chg = 0;
	
	#print "\nref_pos: $index_ref_pos length: $index_len ref: $index_ref allele: $index_allele ".
	#	"count: $index_count cov: $index_cov freq: $index_freq over_annot: $index_over_annot ".
	#	"cod_reg_chg: $index_cod_reg_chg amino acid chg: $index_ami_aci_chg\n";
	
	#print "\n\n$_\n\n";
	
	# Determine the index of each column -- this is CLC
	for (my $i=0; $i<scalar(@header_names_arr); $i++) {
		if ($header_names_arr[$i] =~ /Reference Position/) {
			$index_ref_pos = $i;
			$seen_index_ref_pos = 1;
		} elsif ($header_names_arr[$i] =~ /Type/) {
			$index_type = $i;
			$seen_index_type = 1;
		} elsif ($header_names_arr[$i] =~ /Reference/) { # Since this comes AFTER 
													# "Reference Position, we're fine
			$index_ref = $i;
			$seen_index_ref = 1;
		} elsif ($header_names_arr[$i] =~ /Allele/) {
			$index_allele = $i;
			$seen_index_allele = 1;
		} elsif ($header_names_arr[$i] =~ /Count/) {
			$index_count = $i;
			$seen_index_count = 1;
		} elsif ($header_names_arr[$i] =~ /Coverage/) {
			$index_cov = $i;
			$seen_index_cov = 1;
		} elsif ($header_names_arr[$i] =~ /Frequency/) {
			$index_freq = $i;
			$seen_index_freq = 1;
		} elsif ($header_names_arr[$i] =~ /Overlapping annotations/) {
			$index_over_annot = $i;
			$seen_index_over_annot = 1;
#		} elsif ($header_names_arr[$i] =~ /Coding region change/) {
#			$index_cod_reg_chg = $i;
#			$seen_index_cod_reg_chg = 1;
#		} elsif ($header_names_arr[$i] =~ /Amino acid change/) {
#			$index_ami_aci_chg = $i;
#			$seen_index_ami_aci_chg = 1;
		}
	}
	
	if($seen_index_ref_pos == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Reference Position\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Reference Position\". SNPGenie terminated.\n\n";
	} elsif($seen_index_type == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Type\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Type\". SNPGenie terminated.\n\n";
	} elsif($seen_index_ref == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Reference\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Reference\". SNPGenie terminated.\n\n";
	} elsif($seen_index_allele == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Allele\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Allele\". SNPGenie terminated.\n\n";
	} elsif($seen_index_count == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Count\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Count\". SNPGenie terminated.\n\n";
	} elsif($seen_index_cov == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Coverage\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Coverage\". SNPGenie terminated.\n\n";
	} elsif($seen_index_freq == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Frequency\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Frequency\". SNPGenie terminated.\n\n";
	} elsif($seen_index_over_annot == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tNA\t".
			"Does not contain the column header \"Overlapping annotations\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $file_nm does not contain the column header \"Overlapping annotations\". SNPGenie terminated.\n\n";
	}
	
	#print "\n\n$_\n\n";
	
	#print "\n\nAll my product names are: @product_names_arr\n\n";
	# Now we have a specific file we're looking at ($curr_snp_report_name), and we 
	# know all the names of products within that file (stored in @product_names_arr),
	# and we know all the names of fasta files in this directory (stored in 
	# @fasta_file_names_arr). So we have to identify the products in their names,
	# taking special care with the HA's -- HA, HA1, HA2
	
	# Now we have, for the current SNP Report we're examining, the fasta files to which
	# each of the products refer, in the hash %products_to_fasta_file
	
	# Next, let's store the contents of the SNP Report in a hash of hashes, with the
	# outer key as the product, followed by and making sure to store the name of the 
	# associated fasta file
	
	my %hh_product_position_info; # For THIS SNP REPORT only. SEE EXAMPLE BELOW.
	
	# Early example:
	#my %hh_example = (
	#	'HA' => {
	#		'start' => 666,
	#		'stop' =>,
	#		158 => {
	#			'fasta' => $product_to_fasta_file{'HA'},
	#			'ref' => 'C',
	#			'cov' => 3002,
	#			'A' => 0,
	#			'C' => 0,
	#			'G' => 0,
	#			'T' => 0
	#		},
	#		160 => {
	#			'A' => 0,
	#			'C' => 0,
	#			'G' => 0,
	#			'T' => 0
	#		}
	#	},
	#	'NA' => {
	#		158 => {
	#			'A' => 0,
	#			'C' => 0,
	#			'G' => 0,
	#			'T' => 0
	#		},
	#		160 => {
	#			'A' => 0,
	#			'C' => 0,
	#			'G' => 0,
	#			'T' => 0
	#		}
	#	}
	#);
	
	# For site-by-site plain nucleotide diversity of all sites in the FASTA.
	# No partitioning into synonymous and nonsynonymous.
	my %hh_nc_position_info;
	
	##### DATA STORAGE ##################################################################
	# Open current SNP Report to store data for each line
	print "\nCalculating and storing protein-coding genome and variant data (that dogma stuff)... ";
	my $line = 0;
	#open (INFILE, $curr_snp_report_name);
	while (<$TEMP_SNP_REPORT_HANDLE>) {
		if($line == 0) {
			$line++;
		} else {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			#if ($_ =~/(.+)\r$/) {
			#	$_ = $1;
			#	print "\nYes, in $curr_snp_report_name we have an ending \\r\n";
			#}
			
			my @line_arr = split(/\t/,$_,-1);
			
			#$line++;
			#if ($line_arr[$index_over_annot] eq '') {
			#	print "\nIn $curr_snp_report_name, line $line, line_arr[index_over_annot] ".
			#		"is ne ''";
			#}
			# LESSON: even when the first line has no entry for Overlapping annotations,
			# that datum is still (ne '').
			
			# Check that we're dealing with a single nucleotide variant, or
			# SNV in the "Type" column, and that we have a specific product, meaning
			# the "Overlapping annotations" column is not blank
			if(($line_arr[$index_type] eq 'SNV') && ($line_arr[$index_count] > 0) && ($line_arr[$index_freq] > 0)) {
				# Store data in variables to save room on screen and verify
				my $position = $line_arr[$index_ref_pos];
				my $coverage = $line_arr[$index_cov];
				my $count = $line_arr[$index_count];
				my $reference_nt = $line_arr[$index_ref];
				$reference_nt =~ tr/acgt/ACGT/;
				my $variant_nt = $line_arr[$index_allele];
				$variant_nt =~ tr/acgt/ACGT/;
				my $frequency = $line_arr[$index_freq];
				my $var_prop = ($count/$coverage);
				my $ref_prop = (1-$var_prop);
				
				my $ref_prop_key = "$reference_nt\_prop";
				my $var_prop_key = "$variant_nt\_prop";
				
				#print "\nSNV site $position\. ref_prop_key: $ref_prop_key | ref_prop: $ref_prop | var_prop_key: $var_prop_key | var_prop: $var_prop";
				
				if(length($reference_nt) == length($variant_nt)) {
					# Now, nonsyn and syn nucleotide diversity for coding variants
					if($line_arr[$index_over_annot] =~/CDS/) {
						
						# Get the product name(s) for this line
						my $product_name;
						my $mature_peptide_name;
						my $over_annot = $line_arr[$index_over_annot];
						#print "\n\nover_annot is: $over_annot\n\n";
						
						my @product_coord_arr;
						my @peptide_coord_arr;
						
						# EXPERIMENTAL
						if ($over_annot =~/Mature peptide: ([\w\s\.\-\:\'\"]+)/) {
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							print ERROR_FILE "$file_nm\tNA\tNA\t".
								"\"Mature peptide\" annotation must take into account ".
								"MNV records for CLC SNP Reports\n";
							close ERROR_FILE;
							chdir('..');
							
							print "\n## WARNING: \"Mature peptide\" annotation must take into account\n".
							"### MNV records for CLC SNP Reports.\n";
							
							$mature_peptide_name = $1;
							
							# Process product_name for leading or trailing quote marks
							$mature_peptide_name =~ s/^"//;
							$mature_peptide_name =~ s/^'//;
							$mature_peptide_name =~ s/"$//;
							$mature_peptide_name =~ s/'$//;
							
							@peptide_coord_arr = @{$product_coordinates_harr{$mature_peptide_name}->{product_coord_arr}};
						} 
						
						if ($over_annot =~/CDS: ([\w\s\.\-\:\'\"]+)/) {
							$product_name = $1;
							
							# Process product_name for leading or trailing quote marks
							$product_name =~ s/^"//;
							$product_name =~ s/^'//;
							$product_name =~ s/"$//;
							$product_name =~ s/'$//;
							
							#print "\n\nproduct name: $product_name\n\n";
							@product_coord_arr = @{$product_coordinates_harr{$product_name}->{product_coord_arr}};
						} 
						
						# In some files, we have such atrocities as:
						# ORF: gag
						# ORF: gag, ORF: pol
						# ORF: Nef, ORF: Env, ORF: Rev
						#if ($over_annot =~/ORF: (\w+)/) {
						#	$product_name = $1;
						#	@product_coord_arr = &get_product_coordinates($product_name);
						#}
						#
						#print "The array: @product_coord_arr";
						
						if($seen_percent_warning == 0 && $frequency < 0.01) {
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							print ERROR_FILE "$file_nm\t$product_name\t$position\t".
								"There is a variant frequency <0.01%. If this ".
								"was unexpected, make sure that the values in the ".
								"\"Frequency\" column are percentages, not ".
								"decimals.\n";
							close ERROR_FILE;
							chdir('..');
							
							warn "\n## WARNING There is a variant frequency <0.01%. If this ".
								"was unexpected, make sure that the \n## values in the ".
								"\"Frequency\" column are percentages, not ".
								"decimals.\n";
								
							$seen_percent_warning = 1;
						}
						
						# New segments approach
						my %product_starts;
						my %product_stops;
						
						my $num_segments = (@product_coord_arr / 2);
						for(my $i=1; $i<=$num_segments; $i++) { # $i<=scalar(@product_coord_arr)
							$product_starts{$i} = $product_coord_arr[2*$i-2];
							$product_stops{$i} = $product_coord_arr[2*$i-1];
						}
						
						#my $A_count = 0;
						#my $C_count = 0;
						#my $G_count = 0;
						#my $T_count = 0;
						#print "\n$fasta_file\n";
						#print "\n$coverage\n";
						
						if (! exists $hh_product_position_info{$product_name}->{$position}) { # Product/position HAVEN'T been seen
							# New segments approach
							foreach(sort {$a <=> $b} (keys %product_starts)) {
								my $this_start_key = 'start_' . $_;
								my $this_stop_key = 'stop_' . $_;
								$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
								$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
							}
							$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
							
							$hh_product_position_info{$product_name}->{$position}->{A} = 0;
							$hh_product_position_info{$product_name}->{$position}->{C} = 0;
							$hh_product_position_info{$product_name}->{$position}->{G} = 0;
							$hh_product_position_info{$product_name}->{$position}->{T} = 0;
							$hh_product_position_info{$product_name}->{$position}->{A_prop} = 0;
							$hh_product_position_info{$product_name}->{$position}->{C_prop} = 0;
							$hh_product_position_info{$product_name}->{$position}->{G_prop} = 0;
							$hh_product_position_info{$product_name}->{$position}->{T_prop} = 0;
							$hh_product_position_info{$product_name}->{$position}->{reference} = $reference_nt;
							$hh_product_position_info{$product_name}->{$position}->{cov} = $coverage;
							push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
							
							$hh_nc_position_info{$position}->{A} = 0;
							$hh_nc_position_info{$position}->{C} = 0;
							$hh_nc_position_info{$position}->{G} = 0;
							$hh_nc_position_info{$position}->{T} = 0;
							$hh_nc_position_info{$position}->{A_prop} = 0;
							$hh_nc_position_info{$position}->{C_prop} = 0;
							$hh_nc_position_info{$position}->{G_prop} = 0;
							$hh_nc_position_info{$position}->{T_prop} = 0;
							$hh_nc_position_info{$position}->{reference} = $reference_nt;
							$hh_nc_position_info{$position}->{cov} = $coverage;
							push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
							$hh_nc_position_info{$position}->{polymorphic} = 1;
							
							
							if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
								$hh_product_position_info{$product_name}->{$position}->{$variant_nt} = ($count);
								$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage-($count);
								
								$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} = $var_prop;
								$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = $ref_prop;
								
								$hh_nc_position_info{$position}->{$variant_nt} = ($count);
								$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
								$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
								$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
							} else { # Ref and var nts are identical
								$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage;
								$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = 1;
								
								$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
								$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
							}
							
						} else { # Product/position HAVE been seen before
							# If the variant and reference are the same, nothing required
							if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
								#$hh_product_position_info{$product_name}->{$position}->{$variant_nt} += ($coverage*$frequency/100);
								$hh_product_position_info{$product_name}->{$position}->{$variant_nt} += ($count);
								#$hh_product_position_info{$product_name}->{$position}->{$reference_nt} -= ($coverage*$frequency/100);
								$hh_product_position_info{$product_name}->{$position}->{$reference_nt} -= ($count);
								
								$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} += $var_prop;
								$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} -= $var_prop;
	
								$hh_nc_position_info{$position}->{$variant_nt} += ($count);
								$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
								$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
								$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
							}
							
							push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
							push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
							
							if ($hh_product_position_info{$product_name}->{$position}->{cov} != $coverage) {
								chdir('SNPGenie_Results');
								open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
								print ERROR_FILE "$file_nm\t$product_name\t$position\t".
									"There are conflicting coverages reported. ".
										"An averaging has taken place\n";
								close ERROR_FILE;
								chdir('..');
								
								warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position\n".
										"## An averaging has taken place.\n";
							}
						}
						
						# NOW for the "Mature peptide," if it exists: DEPRECATED
						my %peptide_starts;
						my %peptide_stops;
						
						my $num_segments = (@peptide_coord_arr / 2);
						for(my $i=1; $i<=scalar(@peptide_coord_arr); $i++) {
							$peptide_starts{$i} = $peptide_coord_arr[2*$i-2];
							$peptide_stops{$i} = $peptide_coord_arr[2*$i-1];
						}
						
						# NOTE: Mature peptide is never the same as CDS in examples.
						if ($mature_peptide_name) {
							if (! exists $hh_product_position_info{$mature_peptide_name}->{$position}) { # Product/position HAVEN'T been seen
								chdir('SNPGenie_Results');
								open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
								print ERROR_FILE "$file_nm\t$mature_peptide_name\t$position\t".
									"A \'mature peptide\' annotation is being used. Contact ".
										"the author to update this deprecated function\n";
								close ERROR_FILE;
								chdir('..');
								
								warn "\n## WARNING: A \'mature peptide\' annotation is being used at ".
										"$file_nm|$product_name|$position\n".
										"## Contact the author to update this deprecated function.\n";
								
								# New segments approach
								foreach(sort {$a <=> $b} (keys %peptide_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$mature_peptide_name}->{$this_start_key} = $peptide_starts{$_};
									$hh_product_position_info{$mature_peptide_name}->{$this_stop_key} = $peptide_stops{$_};
								}
								$hh_product_position_info{$mature_peptide_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$mature_peptide_name}->{$position}->{A} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{C} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{G} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{T} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{A_prop} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{C_prop} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{G_prop} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{T_prop} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{reference} = $reference_nt;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{cov} = $coverage;
								
								push(@{$hh_product_position_info{$mature_peptide_name}->{$position}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 1;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									#$hh_product_position_info{$mature_peptide_name}->{$position}->{$variant_nt} = ($coverage*$frequency/100);
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$variant_nt} = ($count);
									#$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} = $coverage-($coverage*$frequency/100);
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} = $coverage-($count);
									
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$var_prop_key} = $var_prop;
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$ref_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} = $coverage;
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$ref_prop_key} = 1;
									
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								# Made this first one += because some SNPs are showing up multiple
								# times in the reports, e.g., 3 entries in Virus_20... for PB1
								# position 982, all G->T
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									#$hh_product_position_info{$mature_peptide_name}->{$position}->{$variant_nt} = ($coverage*$frequency/100);
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$variant_nt} += ($count);
									#$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} = $coverage-($coverage*$frequency/100);
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} -= ($count);
									
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$var_prop_key} += $var_prop;
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$ref_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								}
								
								push(@{$hh_product_position_info{$mature_peptide_name}->{$position}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$mature_peptide_name}->{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$mature_peptide_name\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$mature_peptide_name|$position\n".
										"## An averaging has taken place.\n";
								}
							}
						}
					} elsif($line_arr[$index_over_annot] eq '' || $line_arr[$index_over_annot] =~ "\' UTR") { # NON-CODING SNV VARIANT
						#print "\n\nWe have a non-coding variant at site $position\n";
	
						if(! exists $hh_nc_position_info{$position}) { # First time seeing this 
							# site. First, we populate %hh_position_info for plain old 
							# nucleotide diversity.
							$hh_nc_position_info{$position}->{A} = 0;
							$hh_nc_position_info{$position}->{C} = 0;
							$hh_nc_position_info{$position}->{G} = 0;
							$hh_nc_position_info{$position}->{T} = 0;
							$hh_nc_position_info{$position}->{A_prop} = 0;
							$hh_nc_position_info{$position}->{C_prop} = 0;
							$hh_nc_position_info{$position}->{G_prop} = 0;
							$hh_nc_position_info{$position}->{T_prop} = 0;
							$hh_nc_position_info{$position}->{reference} = $reference_nt;
							$hh_nc_position_info{$position}->{cov} = $coverage;
							push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
							#$hh_nc_position_info{$position}->{coding} = 0;
							$hh_nc_position_info{$position}->{polymorphic} = 1;
							
							if($variant_nt ne $reference_nt) { # Make sure the reference and
								# variant nts aren't identical.
								$hh_nc_position_info{$position}->{$variant_nt} = ($count);
								$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
	
								$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
								$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
							} else { # If the ref and var happen to be the same
								$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
								$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
							}
							#print "\n\nTEST\nhh_nc_position_info for $position:\n".
							#	"A_prop=".$hh_nc_position_info{$position}->{A_prop}."\n".
							#	"C_prop=".$hh_nc_position_info{$position}->{C_prop}."\n".
							#	"G_prop=".$hh_nc_position_info{$position}->{G_prop}."\n".
							#	"T_prop=".$hh_nc_position_info{$position}->{T_prop}."\n";
						} else { # the site has been seen before
							if($variant_nt ne $reference_nt) { # Make sure the reference and 
							# variant nts aren't identical.
								$hh_nc_position_info{$position}->{$variant_nt} += ($count);
								$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
	
								$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
								$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
							} # NO ACTION REQUIRED if the ref and var are the same
							
							push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
							
							if ($hh_nc_position_info{$position}->{cov} != $coverage) {
								chdir('SNPGenie_Results');
								open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
								print ERROR_FILE "$file_nm\tN\/A\t$position\t".
									"There are conflicting coverages reported. ".
										"An averaging has taken place\n";
								close ERROR_FILE;
								chdir('..');
								
								warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N\/A|$position\n".
										"## An averaging has taken place.\n";
							}
							#print "\n\nTEST\nhh_nc_position_info for $position:\n".
							#	"A_prop=".$hh_nc_position_info{$position}->{A_prop}."\n".
							#	"C_prop=".$hh_nc_position_info{$position}->{C_prop}."\n".
							#	"G_prop=".$hh_nc_position_info{$position}->{G_prop}."\n".
							#	"T_prop=".$hh_nc_position_info{$position}->{T_prop}."\n";
						}
					}
				} # for SNV, ref length == variant length
			} elsif(($line_arr[$index_type] eq 'MNV') && ($line_arr[$index_count] > 0) && 
				($line_arr[$index_freq] > 0)) {  # THERE ARE MNV's both 2-5 nt in length
				
				# Store data in variables to save room on screen and verify
				my $position = $line_arr[$index_ref_pos];
				my $coverage = $line_arr[$index_cov];
				my $count = $line_arr[$index_count];
				my $reference_nt = $line_arr[$index_ref];
				$reference_nt =~ tr/acgt/ACGT/;
				my $variant_nt = $line_arr[$index_allele];
				$variant_nt =~ tr/acgt/ACGT/;
				my $frequency = $line_arr[$index_freq];
				my $var_prop = ($count/$coverage);
				#my $var_prop = ($line_arr[$index_freq])/100;
				my $ref_prop = (1-$var_prop);
								
				if(length($reference_nt) == length($variant_nt)) {			
					# Now, syn. and nonsyn. nucleotide diversity for coding variants
					if($line_arr[$index_over_annot] =~/CDS/) { # CODING VARIANT
						# was if(($line_arr[$index_type] eq 'SNV') && ($line_arr[$index_over_annot] ne ''))
						
						if($warn_5nt == 0) {
							warn "\n## WARNING: SNPGenie only considers MNV's up to 5nt in length.\n";
							$warn_5nt = 1;
						}
						
						# Get the product name(s) for this line
						my $product_name;
						my $mature_peptide_name;
						my $over_annot = $line_arr[$index_over_annot];
						
						my @peptide_coord_arr;
						my @product_coord_arr;
									
						# EXPERIMENTAL
						if ($over_annot =~/Mature peptide: ([\w\s\.']+)/) {
							$mature_peptide_name = $1;
							#print "product name is: $product_name\n";
							@peptide_coord_arr = @{$product_coordinates_harr{$mature_peptide_name}->{product_coord_arr}};
							warn "\n## WARNING: \"Mature peptide\" annotation must take into account\n".
							"### MNV records for CLC SNP Reports.\n";
						} 
						
						if ($over_annot =~/CDS: ([\w\s\.\-\:\'\"]+)/) {
							$product_name = $1;
							#print "product name is: $product_name\n";
							
							# Process product_name for leading or trailing quote marks
							$product_name =~ s/^"//;
							$product_name =~ s/^'//;
							$product_name =~ s/"$//;
							$product_name =~ s/'$//;
							
							@product_coord_arr = @{$product_coordinates_harr{$product_name}->{product_coord_arr}};
						}
						
						# New segments approach
						my %product_starts;
						my %product_stops;
						
						my $num_segments = (@product_coord_arr / 2);
						for(my $i=1; $i<=scalar(@product_coord_arr); $i++) {
							$product_starts{$i} = $product_coord_arr[2*$i-2];
							$product_stops{$i} = $product_coord_arr[2*$i-1];
						}
						
						#my $A_count = 0;
						#my $C_count = 0;
						#my $G_count = 0;
						#my $T_count = 0;
						#print "\n$fasta_file\n";
						#print "\n$coverage\n";
						
						if(length($reference_nt) == 2) {
							#print "\nSaw a 2-nt MNV! $reference_nt at $product_name site $position\n\n";
							# Add extra variables with values for the 2-nt MNV
							my $position2 = $position+1;
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt = substr($reference_nt,0,1);
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt = substr($variant_nt,0,1);
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{reference} = $reference_nt;
								$hh_product_position_info{$product_name}->{$position}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 1;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} = ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else { # Ref and var nts are identical
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = 1;
									
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
								
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} += ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} -= ($count);
									
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position\n".
										"## An averaging has taken place.\n";
								}
								
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position2}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position2}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{reference} = $reference_nt2;
								$hh_product_position_info{$product_name}->{$position2}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 1;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} = ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = 1;
									
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
								
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} += ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} -= ($count);
									
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								}  # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
						} elsif(length($reference_nt) == 3) {
							#print "\nSaw a 3-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 3-nt MNV
							my $position2 = $position+1;
							my $position3 = $position+2;
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt3 = substr($reference_nt,2,1);
							my $reference_nt = substr($reference_nt,0,1);
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt3 = substr($variant_nt,2,1);
							my $variant_nt = substr($variant_nt,0,1);
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $ref3_prop_key = "$reference_nt3\_prop";
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							my $var3_prop_key = "$variant_nt3\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{reference} = $reference_nt;
								$hh_product_position_info{$product_name}->{$position}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 1;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} = ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = 1;
									
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} += ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} -= ($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position}->{cov} != $coverage) {					
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position2}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position2}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{reference} = $reference_nt2;
								$hh_product_position_info{$product_name}->{$position2}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 1;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} = ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = 1;
									
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								# Made this first one += because some SNPs are showing up multiple
								# times in the reports, e.g., 3 entries in Virus_20... for PB1
								# position 982, all G->T
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} += ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} -= ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 3 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position3}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position3}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{reference} = $reference_nt3;
								$hh_product_position_info{$product_name}->{$position3}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position3}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position3}->{A} = 0;
								$hh_nc_position_info{$position3}->{C} = 0;
								$hh_nc_position_info{$position3}->{G} = 0;
								$hh_nc_position_info{$position3}->{T} = 0;
								$hh_nc_position_info{$position3}->{A_prop} = 0;
								$hh_nc_position_info{$position3}->{C_prop} = 0;
								$hh_nc_position_info{$position3}->{G_prop} = 0;
								$hh_nc_position_info{$position3}->{T_prop} = 0;
								$hh_nc_position_info{$position3}->{reference} = $reference_nt3;
								$hh_nc_position_info{$position3}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position3}->{coding} = 1;
								$hh_nc_position_info{$position3}->{polymorphic} = 1;
								
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position3}->{$variant_nt3} = ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position3}->{$var3_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position3}->{$variant_nt3} = ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} = $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} = $coverage;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} = 1;
									
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position3}->{$variant_nt3} += ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} -= ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$var3_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position3}->{$variant_nt3} += ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} -= ($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} += $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position3}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position3}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position3\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position3\n".
										"## An averaging has taken place.\n";
								}
							}
							
						} elsif(length($reference_nt) == 4) {
							#print "\nSaw a 3-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 3-nt MNV
							my $position2 = $position+1;
							my $position3 = $position+2;
							my $position4 = $position+3;
							
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt3 = substr($reference_nt,2,1);
							my $reference_nt4 = substr($reference_nt,3,1);
							my $reference_nt = substr($reference_nt,0,1);
							
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt3 = substr($variant_nt,2,1);
							my $variant_nt4 = substr($variant_nt,3,1);
							my $variant_nt = substr($variant_nt,0,1);
							
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $ref3_prop_key = "$reference_nt3\_prop";
							my $ref4_prop_key = "$reference_nt4\_prop";
							
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							my $var3_prop_key = "$variant_nt3\_prop";
							my $var4_prop_key = "$variant_nt4\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{reference} = $reference_nt;
								$hh_product_position_info{$product_name}->{$position}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 1;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} = ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = 1;
									
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} += ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} -= ($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position2}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position2}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{reference} = $reference_nt2;
								$hh_product_position_info{$product_name}->{$position2}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 1;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} = ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = 1;
									
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} += ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} -= ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 3 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position3}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position3}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{reference} = $reference_nt3;
								$hh_product_position_info{$product_name}->{$position3}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position3}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position3}->{A} = 0;
								$hh_nc_position_info{$position3}->{C} = 0;
								$hh_nc_position_info{$position3}->{G} = 0;
								$hh_nc_position_info{$position3}->{T} = 0;
								$hh_nc_position_info{$position3}->{A_prop} = 0;
								$hh_nc_position_info{$position3}->{C_prop} = 0;
								$hh_nc_position_info{$position3}->{G_prop} = 0;
								$hh_nc_position_info{$position3}->{T_prop} = 0;
								$hh_nc_position_info{$position3}->{reference} = $reference_nt3;
								$hh_nc_position_info{$position3}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position3}->{coding} = 1;
								$hh_nc_position_info{$position3}->{polymorphic} = 1;
								
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position3}->{$variant_nt3} = ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position3}->{$var3_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position3}->{$variant_nt3} = ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} = $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} = $coverage;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} = 1;
									
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position3}->{$variant_nt3} += ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} -= ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$var3_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position3}->{$variant_nt3} += ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} -= ($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} += $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position3}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position3}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position3\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position3\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 4 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position4}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position4}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{reference} = $reference_nt4;
								$hh_product_position_info{$product_name}->{$position4}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position4}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position4}->{A} = 0;
								$hh_nc_position_info{$position4}->{C} = 0;
								$hh_nc_position_info{$position4}->{G} = 0;
								$hh_nc_position_info{$position4}->{T} = 0;
								$hh_nc_position_info{$position4}->{A_prop} = 0;
								$hh_nc_position_info{$position4}->{C_prop} = 0;
								$hh_nc_position_info{$position4}->{G_prop} = 0;
								$hh_nc_position_info{$position4}->{T_prop} = 0;
								$hh_nc_position_info{$position4}->{reference} = $reference_nt4;
								$hh_nc_position_info{$position4}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position4}->{coding} = 1;
								$hh_nc_position_info{$position4}->{polymorphic} = 1;
								
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position4}->{$variant_nt4} = ($count);
									$hh_product_position_info{$product_name}->{$position4}->{$reference_nt4} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position4}->{$var4_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position4}->{$ref4_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position4}->{$variant_nt4} = ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage-($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} = $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position4}->{$reference_nt4} = $coverage;
									$hh_product_position_info{$product_name}->{$position4}->{$ref4_prop_key} = 1;
									
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position4}->{$variant_nt4} += ($count);
									$hh_product_position_info{$product_name}->{$position4}->{$reference_nt4} -= ($count);
									$hh_product_position_info{$product_name}->{$position4}->{$var4_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position4}->{$ref4_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position4}->{$variant_nt4} += ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} -= ($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} += $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position4}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position4}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position4\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position4\n".
										"## An averaging has taken place.\n";
								}
							}
							
						} elsif(length($reference_nt) == 5) {
							#print "\nSaw a 3-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 3-nt MNV
							my $position2 = $position+1;
							my $position3 = $position+2;
							my $position4 = $position+3;
							my $position5 = $position+4;
							
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt3 = substr($reference_nt,2,1);
							my $reference_nt4 = substr($reference_nt,3,1);
							my $reference_nt5 = substr($reference_nt,4,1);
							my $reference_nt = substr($reference_nt,0,1);
							
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt3 = substr($variant_nt,2,1);
							my $variant_nt4 = substr($variant_nt,3,1);
							my $variant_nt5 = substr($variant_nt,4,1);
							my $variant_nt = substr($variant_nt,0,1);
							
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $ref3_prop_key = "$reference_nt3\_prop";
							my $ref4_prop_key = "$reference_nt4\_prop";
							my $ref5_prop_key = "$reference_nt5\_prop";
							
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							my $var3_prop_key = "$variant_nt3\_prop";
							my $var4_prop_key = "$variant_nt4\_prop";
							my $var5_prop_key = "$variant_nt5\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position}->{reference} = $reference_nt;
								$hh_product_position_info{$product_name}->{$position}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 1;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} = ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} = 1;
									
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position}->{$variant_nt} += ($count);
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} -= $coverage-($count);
									$hh_product_position_info{$product_name}->{$position}->{$var_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position}->{$ref_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position2}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position2}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position2}->{reference} = $reference_nt2;
								$hh_product_position_info{$product_name}->{$position2}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 1;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} = ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} = $coverage;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} = 1;
									
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position2}->{$variant_nt2} += ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$reference_nt2} -= ($count);
									$hh_product_position_info{$product_name}->{$position2}->{$var2_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position2}->{$ref2_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position2}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 3 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position3}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position3}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position3}->{reference} = $reference_nt3;
								$hh_product_position_info{$product_name}->{$position3}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position3}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position3}->{A} = 0;
								$hh_nc_position_info{$position3}->{C} = 0;
								$hh_nc_position_info{$position3}->{G} = 0;
								$hh_nc_position_info{$position3}->{T} = 0;
								$hh_nc_position_info{$position3}->{A_prop} = 0;
								$hh_nc_position_info{$position3}->{C_prop} = 0;
								$hh_nc_position_info{$position3}->{G_prop} = 0;
								$hh_nc_position_info{$position3}->{T_prop} = 0;
								$hh_nc_position_info{$position3}->{reference} = $reference_nt3;
								$hh_nc_position_info{$position3}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position3}->{coding} = 1;
								$hh_nc_position_info{$position3}->{polymorphic} = 1;
								
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position3}->{$variant_nt3} = ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position3}->{$var3_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position3}->{$variant_nt3} = ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} = $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} = $coverage;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} = 1;
									
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position3}->{$variant_nt3} += ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$reference_nt3} -= ($count);
									$hh_product_position_info{$product_name}->{$position3}->{$var3_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position3}->{$ref3_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position3}->{$variant_nt3} += ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} -= ($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} += $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position3}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position3}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position3\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position3\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 4 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position4}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position4}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position4}->{reference} = $reference_nt4;
								$hh_product_position_info{$product_name}->{$position4}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position4}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position4}->{A} = 0;
								$hh_nc_position_info{$position4}->{C} = 0;
								$hh_nc_position_info{$position4}->{G} = 0;
								$hh_nc_position_info{$position4}->{T} = 0;
								$hh_nc_position_info{$position4}->{A_prop} = 0;
								$hh_nc_position_info{$position4}->{C_prop} = 0;
								$hh_nc_position_info{$position4}->{G_prop} = 0;
								$hh_nc_position_info{$position4}->{T_prop} = 0;
								$hh_nc_position_info{$position4}->{reference} = $reference_nt4;
								$hh_nc_position_info{$position4}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position4}->{coding} = 1;
								$hh_nc_position_info{$position4}->{polymorphic} = 1;
								
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position4}->{$variant_nt4} = ($count);
									$hh_product_position_info{$product_name}->{$position4}->{$reference_nt4} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position4}->{$var4_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position4}->{$ref4_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position4}->{$variant_nt4} = ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage-($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} = $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position4}->{$reference_nt4} = $coverage;
									$hh_product_position_info{$product_name}->{$position4}->{$ref4_prop_key} = 1;
									
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position4}->{$variant_nt4} += ($count);
									$hh_product_position_info{$product_name}->{$position4}->{$reference_nt4} -= ($count);
									$hh_product_position_info{$product_name}->{$position4}->{$var4_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position4}->{$ref4_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position4}->{$variant_nt4} += ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} -= ($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} += $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_product_position_info{$product_name}->{$position4}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position4}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position4\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position4\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 5 of the variant
							if (! exists $hh_product_position_info{$product_name}->{$position5}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %product_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
									$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
								}
								$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$product_name}->{$position5}->{A} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{C} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{G} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{T} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{A_prop} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{C_prop} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{G_prop} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{T_prop} = 0;
								$hh_product_position_info{$product_name}->{$position5}->{reference} = $reference_nt5;
								$hh_product_position_info{$product_name}->{$position5}->{cov} = $coverage;
								push(@{$hh_product_position_info{$product_name}->{$position5}->{cov_arr}},$coverage);
								
								$hh_nc_position_info{$position5}->{A} = 0;
								$hh_nc_position_info{$position5}->{C} = 0;
								$hh_nc_position_info{$position5}->{G} = 0;
								$hh_nc_position_info{$position5}->{T} = 0;
								$hh_nc_position_info{$position5}->{A_prop} = 0;
								$hh_nc_position_info{$position5}->{C_prop} = 0;
								$hh_nc_position_info{$position5}->{G_prop} = 0;
								$hh_nc_position_info{$position5}->{T_prop} = 0;
								$hh_nc_position_info{$position5}->{reference} = $reference_nt5;
								$hh_nc_position_info{$position5}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position5}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position5}->{coding} = 1;
								$hh_nc_position_info{$position5}->{polymorphic} = 1;
								
								if($variant_nt5 ne $reference_nt5) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position5}->{$variant_nt5} = ($count);
									$hh_product_position_info{$product_name}->{$position5}->{$reference_nt5} = $coverage-($count);
									$hh_product_position_info{$product_name}->{$position5}->{$var5_prop_key} = $var_prop;
									$hh_product_position_info{$product_name}->{$position5}->{$ref5_prop_key} = $ref_prop;
									
									$hh_nc_position_info{$position5}->{$variant_nt5} = ($count);
									$hh_nc_position_info{$position5}->{$reference_nt5} = $coverage-($count);
									$hh_nc_position_info{$position5}->{$var5_prop_key} = $var_prop;
									$hh_nc_position_info{$position5}->{$ref5_prop_key} = $ref_prop;
								} else {
									$hh_product_position_info{$product_name}->{$position5}->{$reference_nt5} = $coverage;
									$hh_product_position_info{$product_name}->{$position5}->{$ref5_prop_key} = 1;
									
									$hh_nc_position_info{$position5}->{$reference_nt5} = $coverage;
									$hh_nc_position_info{$position5}->{$ref5_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt5 ne $reference_nt5) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$product_name}->{$position5}->{$variant_nt5} += ($count);
									$hh_product_position_info{$product_name}->{$position5}->{$reference_nt5} -= ($count);
									$hh_product_position_info{$product_name}->{$position5}->{$var5_prop_key} += $var_prop;
									$hh_product_position_info{$product_name}->{$position5}->{$ref5_prop_key} -= $var_prop;
									
									$hh_nc_position_info{$position5}->{$variant_nt5} += ($count);
									$hh_nc_position_info{$position5}->{$reference_nt5} -= ($count);
									$hh_nc_position_info{$position5}->{$var5_prop_key} += $var_prop;
									$hh_nc_position_info{$position5}->{$ref5_prop_key} -= $var_prop;
								}
								
								push(@{$hh_product_position_info{$product_name}->{$position5}->{cov_arr}},$coverage);
								push(@{$hh_nc_position_info{$position5}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$product_name}->{$position5}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$product_name\t$position5\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$product_name|$position5\n".
										"## An averaging has taken place.\n";
								}
							}
						}
						
						#if($curr_snp_report_name eq 'SIVkrc_RC08_mapping_0.5_0.7 variants_new.txt' && $position == 6940) {
						#	print "\n6940 T = ".$hh_product_position_info{$product_name}->{$position+2}->{T}."\n\n";
						#}
						
						# New segments approach for mature peptide
						my %peptide_starts;
						my %peptide_stops;
						
						my $num_segments = (@peptide_coord_arr / 2);
						for(my $i=1; $i<=scalar(@peptide_coord_arr); $i++) {
							$peptide_starts{$i} = $peptide_coord_arr[2*$i-2];
							$peptide_stops{$i} = $peptide_coord_arr[2*$i-1];
						}
						
						# NOTE: Mature peptide is never the same as CDS. If Mature peptide is
						# present, then it is HA1, and CDS is HA
						if ($mature_peptide_name) {
							warn "/n## WARNING: \"Mature peptide\" annotation must take into account\n".
							"### MNV records for CLC SNP Reports.\n";
							if (! exists $hh_product_position_info{$mature_peptide_name}->{$position}) { # Product/position HAVEN'T been seen
								# New segments approach
								foreach(sort {$a <=> $b} (keys %peptide_starts)) {
									my $this_start_key = 'start_' . $_;
									my $this_stop_key = 'stop_' . $_;
									$hh_product_position_info{$mature_peptide_name}->{$this_start_key} = $peptide_starts{$_};
									$hh_product_position_info{$mature_peptide_name}->{$this_stop_key} = $peptide_stops{$_};
								}
								$hh_product_position_info{$mature_peptide_name}->{num_segments} = $num_segments;
								
								$hh_product_position_info{$mature_peptide_name}->{$position}->{A} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{C} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{G} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{T} = 0;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{reference} = $reference_nt;
								$hh_product_position_info{$mature_peptide_name}->{$position}->{cov} = $coverage;
								
								# add to coding?
								
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$variant_nt} = ($count);
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} = $coverage-($count);
								} else {
									$hh_product_position_info{$product_name}->{$position}->{$reference_nt} = $coverage;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$variant_nt} += ($count);
									$hh_product_position_info{$mature_peptide_name}->{$position}->{$reference_nt} -= ($count);
								}
								
								push(@{$hh_product_position_info{$product_name}->{$position}->{cov_arr}},$coverage);
								
								if ($hh_product_position_info{$mature_peptide_name}->{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\t$mature_peptide_name\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|$mature_peptide_name|$position\n".
										"## An averaging has taken place.\n";
								}
							}
						}	
					
					# NON-CODING MULTI-NUCLEOTIDE VARIANTS (nc MNVs)	
					} elsif($line_arr[$index_over_annot] eq '' || $line_arr[$index_over_annot] =~ "\' UTR") { # NON-CODING VARIANT
						#print "\n\nWe have a non-coding variant at site $position\n";
						
						if($warn_5nt == 0) {
							warn "\n## WARNING: SNPGenie only considers MNV's up to 5nt in length.\n";
							$warn_5nt = 1;
						}
						
						#my $fasta_file = $product_to_fasta_file{$product_name}; # SHOULDN'T MATTER
						#my $fasta_file = $the_fasta_file; # COMEBACK [why?]
						
						#my $A_count = 0;
						#my $C_count = 0;
						#my $G_count = 0;
						#my $T_count = 0;
						#print "\n$fasta_file\n";
						#print "\n$coverage\n";
						
						if(length($reference_nt) == 2) {
							#print "\nSaw a 2-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 2-nt MNV
							my $position2 = $position+1;
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt = substr($reference_nt,0,1);
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt = substr($variant_nt,0,1);
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_nc_position_info{$position}) { # Product/position HAVEN'T been seen		#########
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 0;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else { # Ref and var nts are identical
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
								
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_nc_position_info{$position2}) { # Product/position HAVEN'T been seen									
								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 0;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
								
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								}  # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
						} elsif(length($reference_nt) == 3) {  # END MNV length 2 #####
							#print "\nSaw a 3-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 3-nt MNV
							my $position2 = $position+1;
							my $position3 = $position+2;
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt3 = substr($reference_nt,2,1);
							my $reference_nt = substr($reference_nt,0,1);
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt3 = substr($variant_nt,2,1);
							my $variant_nt = substr($variant_nt,0,1);
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $ref3_prop_key = "$reference_nt3\_prop";
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							my $var3_prop_key = "$variant_nt3\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_nc_position_info{$position}) { # Product/position HAVEN'T been seen
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 0;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position}->{cov} != $coverage) {					
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_nc_position_info{$position2}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 0;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 3 of the variant
							if (! exists $hh_nc_position_info{$position3}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position3}->{A} = 0;
								$hh_nc_position_info{$position3}->{C} = 0;
								$hh_nc_position_info{$position3}->{G} = 0;
								$hh_nc_position_info{$position3}->{T} = 0;
								$hh_nc_position_info{$position3}->{A_prop} = 0;
								$hh_nc_position_info{$position3}->{C_prop} = 0;
								$hh_nc_position_info{$position3}->{G_prop} = 0;
								$hh_nc_position_info{$position3}->{T_prop} = 0;
								$hh_nc_position_info{$position3}->{reference} = $reference_nt3;
								$hh_nc_position_info{$position3}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position3}->{coding} = 0;
								$hh_nc_position_info{$position3}->{polymorphic} = 1;
								
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position3}->{$variant_nt3} = ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} = $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position3}->{$variant_nt3} += ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} -= ($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} += $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position3}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position3\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position3\n".
										"## An averaging has taken place.\n";
								}
							}
							
						} elsif(length($reference_nt) == 4) {
							#print "\nSaw a 3-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 3-nt MNV
							my $position2 = $position+1;
							my $position3 = $position+2;
							my $position4 = $position+3;
							
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt3 = substr($reference_nt,2,1);
							my $reference_nt4 = substr($reference_nt,3,1);
							my $reference_nt = substr($reference_nt,0,1);
							
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt3 = substr($variant_nt,2,1);
							my $variant_nt4 = substr($variant_nt,3,1);
							my $variant_nt = substr($variant_nt,0,1);
							
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $ref3_prop_key = "$reference_nt3\_prop";
							my $ref4_prop_key = "$reference_nt4\_prop";
							
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							my $var3_prop_key = "$variant_nt3\_prop";
							my $var4_prop_key = "$variant_nt4\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_nc_position_info{$position}) { # Product/position HAVEN'T been seen
								
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 0;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= ($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_nc_position_info{$position2}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 0;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 3 of the variant
							if (! exists $hh_nc_position_info{$position3}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position3}->{A} = 0;
								$hh_nc_position_info{$position3}->{C} = 0;
								$hh_nc_position_info{$position3}->{G} = 0;
								$hh_nc_position_info{$position3}->{T} = 0;
								$hh_nc_position_info{$position3}->{A_prop} = 0;
								$hh_nc_position_info{$position3}->{C_prop} = 0;
								$hh_nc_position_info{$position3}->{G_prop} = 0;
								$hh_nc_position_info{$position3}->{T_prop} = 0;
								$hh_nc_position_info{$position3}->{reference} = $reference_nt3;
								$hh_nc_position_info{$position3}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position3}->{coding} = 0;
								$hh_nc_position_info{$position3}->{polymorphic} = 1;
								
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position3}->{$variant_nt3} = ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} = $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position3}->{$variant_nt3} += ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} -= ($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} += $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position3}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position3\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position3\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 4 of the variant
							if (! exists $hh_nc_position_info{$position4}) { # Product/position HAVEN'T been seen
								
								$hh_nc_position_info{$position4}->{A} = 0;
								$hh_nc_position_info{$position4}->{C} = 0;
								$hh_nc_position_info{$position4}->{G} = 0;
								$hh_nc_position_info{$position4}->{T} = 0;
								$hh_nc_position_info{$position4}->{A_prop} = 0;
								$hh_nc_position_info{$position4}->{C_prop} = 0;
								$hh_nc_position_info{$position4}->{G_prop} = 0;
								$hh_nc_position_info{$position4}->{T_prop} = 0;
								$hh_nc_position_info{$position4}->{reference} = $reference_nt4;
								$hh_nc_position_info{$position4}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position4}->{coding} = 0;
								$hh_nc_position_info{$position4}->{polymorphic} = 1;
								
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position4}->{$variant_nt4} = ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage-($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} = $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = $ref_prop;
									
								} else {
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position4}->{$variant_nt4} += ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} -= ($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} += $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position4}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position4\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position4\n".
										"## An averaging has taken place.\n";
								}
							}

						} elsif(length($reference_nt) == 5) {
							#print "\nSaw a 3-nt MNV! $reference_nt at $product_name site $position\n\n";
							
							# Add extra variables with values for the 3-nt MNV
							my $position2 = $position+1;
							my $position3 = $position+2;
							my $position4 = $position+3;
							my $position5 = $position+4;
							
							my $reference_nt2 = substr($reference_nt,1,1);
							my $reference_nt3 = substr($reference_nt,2,1);
							my $reference_nt4 = substr($reference_nt,3,1);
							my $reference_nt5 = substr($reference_nt,4,1);
							my $reference_nt = substr($reference_nt,0,1);
							
							my $variant_nt2 = substr($variant_nt,1,1);
							my $variant_nt3 = substr($variant_nt,2,1);
							my $variant_nt4 = substr($variant_nt,3,1);
							my $variant_nt5 = substr($variant_nt,4,1);
							my $variant_nt = substr($variant_nt,0,1);
							
							my $ref_prop_key = "$reference_nt\_prop";
							my $ref2_prop_key = "$reference_nt2\_prop";
							my $ref3_prop_key = "$reference_nt3\_prop";
							my $ref4_prop_key = "$reference_nt4\_prop";
							my $ref5_prop_key = "$reference_nt5\_prop";
							
							my $var_prop_key = "$variant_nt\_prop";
							my $var2_prop_key = "$variant_nt2\_prop";
							my $var3_prop_key = "$variant_nt3\_prop";
							my $var4_prop_key = "$variant_nt4\_prop";
							my $var5_prop_key = "$variant_nt5\_prop";
							
							# Do for nucleotide 1 of the variant
							if (! exists $hh_nc_position_info{$position}) { # Product/position HAVEN'T been seen
							
								$hh_nc_position_info{$position}->{A} = 0;
								$hh_nc_position_info{$position}->{C} = 0;
								$hh_nc_position_info{$position}->{G} = 0;
								$hh_nc_position_info{$position}->{T} = 0;
								$hh_nc_position_info{$position}->{A_prop} = 0;
								$hh_nc_position_info{$position}->{C_prop} = 0;
								$hh_nc_position_info{$position}->{G_prop} = 0;
								$hh_nc_position_info{$position}->{T_prop} = 0;
								$hh_nc_position_info{$position}->{reference} = $reference_nt;
								$hh_nc_position_info{$position}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position}->{coding} = 0;
								$hh_nc_position_info{$position}->{polymorphic} = 1;
								
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} = ($count);
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} = $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position}->{$reference_nt} = $coverage;
									$hh_nc_position_info{$position}->{$ref_prop_key} = 1;
								}
								
							} else { # Product/position HAVE been seen before
								if($variant_nt ne $reference_nt) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position}->{$variant_nt} += ($count);
									$hh_nc_position_info{$position}->{$reference_nt} -= $coverage-($count);
									$hh_nc_position_info{$position}->{$var_prop_key} += $var_prop;
									$hh_nc_position_info{$position}->{$ref_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 2 of the variant
							if (! exists $hh_nc_position_info{$position2}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position2}->{A} = 0;
								$hh_nc_position_info{$position2}->{C} = 0;
								$hh_nc_position_info{$position2}->{G} = 0;
								$hh_nc_position_info{$position2}->{T} = 0;
								$hh_nc_position_info{$position2}->{A_prop} = 0;
								$hh_nc_position_info{$position2}->{C_prop} = 0;
								$hh_nc_position_info{$position2}->{G_prop} = 0;
								$hh_nc_position_info{$position2}->{T_prop} = 0;
								$hh_nc_position_info{$position2}->{reference} = $reference_nt2;
								$hh_nc_position_info{$position2}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position2}->{coding} = 0;
								$hh_nc_position_info{$position2}->{polymorphic} = 1;
								
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} = ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage-($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} = $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position2}->{$reference_nt2} = $coverage;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt2 ne $reference_nt2) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position2}->{$variant_nt2} += ($count);
									$hh_nc_position_info{$position2}->{$reference_nt2} -= ($count);
									$hh_nc_position_info{$position2}->{$var2_prop_key} += $var_prop;
									$hh_nc_position_info{$position2}->{$ref2_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position2}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position2}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position2\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position2\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 3 of the variant
							if (! exists $hh_nc_position_info{$position3}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position3}->{A} = 0;
								$hh_nc_position_info{$position3}->{C} = 0;
								$hh_nc_position_info{$position3}->{G} = 0;
								$hh_nc_position_info{$position3}->{T} = 0;
								$hh_nc_position_info{$position3}->{A_prop} = 0;
								$hh_nc_position_info{$position3}->{C_prop} = 0;
								$hh_nc_position_info{$position3}->{G_prop} = 0;
								$hh_nc_position_info{$position3}->{T_prop} = 0;
								$hh_nc_position_info{$position3}->{reference} = $reference_nt3;
								$hh_nc_position_info{$position3}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position3}->{coding} = 0;
								$hh_nc_position_info{$position3}->{polymorphic} = 1;
								
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position3}->{$variant_nt3} = ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage-($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} = $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position3}->{$reference_nt3} = $coverage;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt3 ne $reference_nt3) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position3}->{$variant_nt3} += ($count);
									$hh_nc_position_info{$position3}->{$reference_nt3} -= ($count);
									$hh_nc_position_info{$position3}->{$var3_prop_key} += $var_prop;
									$hh_nc_position_info{$position3}->{$ref3_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position3}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position3}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position3\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position3\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 4 of the variant
							if (! exists $hh_nc_position_info{$position4}) { # Product/position HAVEN'T been seen

								$hh_nc_position_info{$position4}->{A} = 0;
								$hh_nc_position_info{$position4}->{C} = 0;
								$hh_nc_position_info{$position4}->{G} = 0;
								$hh_nc_position_info{$position4}->{T} = 0;
								$hh_nc_position_info{$position4}->{A_prop} = 0;
								$hh_nc_position_info{$position4}->{C_prop} = 0;
								$hh_nc_position_info{$position4}->{G_prop} = 0;
								$hh_nc_position_info{$position4}->{T_prop} = 0;
								$hh_nc_position_info{$position4}->{reference} = $reference_nt4;
								$hh_nc_position_info{$position4}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position4}->{coding} = 0;
								$hh_nc_position_info{$position4}->{polymorphic} = 1;
								
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position4}->{$variant_nt4} = ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage-($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} = $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position4}->{$reference_nt4} = $coverage;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt4 ne $reference_nt4) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position4}->{$variant_nt4} += ($count);
									$hh_nc_position_info{$position4}->{$reference_nt4} -= ($count);
									$hh_nc_position_info{$position4}->{$var4_prop_key} += $var_prop;
									$hh_nc_position_info{$position4}->{$ref4_prop_key} -= $var_prop;
								} # NO ACTION REQUIRED if the ref and var are the same
								
								push(@{$hh_nc_position_info{$position4}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position4}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position4\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position4\n".
										"## An averaging has taken place.\n";
								}
							}
							
							# Do for nucleotide 5 of the variant
							if (! exists $hh_nc_position_info{$position5}) { # Product/position HAVEN'T been seen before
								
								$hh_nc_position_info{$position5}->{A} = 0;
								$hh_nc_position_info{$position5}->{C} = 0;
								$hh_nc_position_info{$position5}->{G} = 0;
								$hh_nc_position_info{$position5}->{T} = 0;
								$hh_nc_position_info{$position5}->{A_prop} = 0;
								$hh_nc_position_info{$position5}->{C_prop} = 0;
								$hh_nc_position_info{$position5}->{G_prop} = 0;
								$hh_nc_position_info{$position5}->{T_prop} = 0;
								$hh_nc_position_info{$position5}->{reference} = $reference_nt5;
								$hh_nc_position_info{$position5}->{cov} = $coverage;
								push(@{$hh_nc_position_info{$position5}->{cov_arr}},$coverage);
								#$hh_nc_position_info{$position5}->{coding} = 0;
								$hh_nc_position_info{$position5}->{polymorphic} = 1;
								
								if($variant_nt5 ne $reference_nt5) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position5}->{$variant_nt5} = ($count);
									$hh_nc_position_info{$position5}->{$reference_nt5} = $coverage-($count);
									$hh_nc_position_info{$position5}->{$var5_prop_key} = $var_prop;
									$hh_nc_position_info{$position5}->{$ref5_prop_key} = $ref_prop;
								} else {
									$hh_nc_position_info{$position5}->{$reference_nt5} = $coverage;
									$hh_nc_position_info{$position5}->{$ref5_prop_key} = 1;
								}
							} else { # Product/position HAVE been seen before
								if($variant_nt5 ne $reference_nt5) { # Make sure the reference and variant nts aren't identical
									$hh_nc_position_info{$position5}->{$variant_nt5} += ($count);
									$hh_nc_position_info{$position5}->{$reference_nt5} -= ($count);
									$hh_nc_position_info{$position5}->{$var5_prop_key} += $var_prop;
									$hh_nc_position_info{$position5}->{$ref5_prop_key} -= $var_prop;
								}
								
								push(@{$hh_nc_position_info{$position5}->{cov_arr}},$coverage);
								
								if ($hh_nc_position_info{$position5}->{cov} != $coverage) {
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									print ERROR_FILE "$file_nm\tN/A\t$position5\t".
										"There are conflicting coverages reported. ".
											"An averaging has taken place\n";
									close ERROR_FILE;
									chdir('..');
									
									warn "\n## WARNING: Conflicting coverages reported at ".
										"$file_nm|N/A|$position5\n".
										"## An averaging has taken place.\n";
								}
							}
						} # END MNV of LENGTH 5	
					} # NON-CODING variant, i.e., no product, blank ''
				} # for MNV, ref length == variant length	
			} # MNV
		} # not line 1
	} # end of SNP Report, end of work populating the multidimensional hash
	close $TEMP_SNP_REPORT_HANDLE;
	print "COMPLETED.\n";
	
	print "\nLigating gene segments (computationally, of course!)... ";
	# New segments approach and streamlining for THIS SNP REPORT
	foreach (@product_names_arr) {
		my $product_name = $_;
		
		if(! exists $hh_product_position_info{$product_name}) {
			#print "product name is: $product_name\n";
			
			# Retrieve coordinates
			my @product_coord_arr = @{$product_coordinates_harr{$product_name}->{product_coord_arr}};
			
			my %product_starts;
			my %product_stops;
			
			my $num_segments = (@product_coord_arr / 2);
			for(my $i=1; $i<=scalar(@product_coord_arr); $i++) {
				$product_starts{$i} = $product_coord_arr[2*$i-2];
				$product_stops{$i} = $product_coord_arr[2*$i-1];
			}
			
			# Store start and stop for all segments
			foreach(sort {$a <=> $b} (keys %product_starts)) {
				my $this_start_key = 'start_' . $_;
				my $this_stop_key = 'stop_' . $_;
				$hh_product_position_info{$product_name}->{$this_start_key} = $product_starts{$_};
				$hh_product_position_info{$product_name}->{$this_stop_key} = $product_stops{$_};
			}
			# Store number of segments
			$hh_product_position_info{$product_name}->{num_segments} = $num_segments;
		}
	}
	
	#my @hh_keys = keys %hh_product_position_info;
	#print "\n\n @hh_keys \n\n";
	
	my @curr_products = sort(keys %hh_product_position_info);
	#print "\n@curr_products";
	
	# Ordering sorting products by the first segment's start site in the genome
	my @curr_products_ordered_by_start = 
		sort { $hh_product_position_info{$a}->{start_1} <=> $hh_product_position_info{$b}->{start_1} } keys %hh_product_position_info;
	
	#print "\nProducts sorted by start sites:\n";
	#foreach (@curr_products_ordered_by_start) {
	#	print "$_\t" . $hh_product_position_info{$_}->{start_1} . "\n";
	#}
	
	print "COMPLETED.\n";
	
	my $product_counter = 0;
	my $pop_summary_line = '';
	
	print "\nCalculating and storing non-protein-coding genome and variant data (\"junk\" gets a bad rap!)... ";
	# NONCODING: CALCULATE AVERAGE COVERAGE at each of the sites in this product.
	# We are REPLACING the coverage currently stored with this.
	# FIRST, extract the positions that have already been stored, i.e., were present
	# in the SNP report
	my @stored_positions_sorted = sort {$a <=> $b} (keys %hh_nc_position_info); # only the polymorphic
	#print "\n\nSo far, we have stored: @stored_positions_sorted\n";
	foreach my $curr_spot (@stored_positions_sorted) {
		my $cov_sum = 0;
		my $cov_denom = 0;
		
		foreach my $cov_mm (@{$hh_nc_position_info{$curr_spot}->{cov_arr}}) {
			$cov_sum += $cov_mm;
			$cov_denom += 1;
		}
		
		my $this_avg_cov = ($cov_sum / $cov_denom);
		$hh_nc_position_info{$curr_spot}->{cov} = $this_avg_cov;
		
		# DETERMINE THE MAJORITY NUCLEOTIDE, TOO
		my $A = $hh_nc_position_info{$curr_spot}->{A};
		my $C = $hh_nc_position_info{$curr_spot}->{C};
		my $G = $hh_nc_position_info{$curr_spot}->{G};
		my $T = $hh_nc_position_info{$curr_spot}->{T};
		
		my $majority_nucleotide; # a variant nucleotide may have fixed
		my $curr_majority_count = 0;
		if($A > $curr_majority_count) {
			$curr_majority_count = $A;
			$majority_nucleotide = 'A';
		}
		if($C > $curr_majority_count) {
			$curr_majority_count = $C;
			$majority_nucleotide = 'C';
		} 
		if($G > $curr_majority_count) {
			$curr_majority_count = $G;
			$majority_nucleotide = 'G';
		} 
		if($T > $curr_majority_count) {
			$curr_majority_count = $T;
			$majority_nucleotide = 'T';
		}
		
		$hh_nc_position_info{$curr_spot}->{maj_nt} = $majority_nucleotide;
		
		#print "\n\nAt site $curr_spot, the majority nucleotide is $majority_nucleotide and ".
		#	"the average coverage is $this_avg_cov\n";

		# WARNING if there is a negative number of sites because of conflicting
		# coverages; round negative nucleotide counts to 0
		if($A < 0) {
			$hh_nc_position_info{$curr_spot}->{A} = 0;
			
			if($A < -0.1) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative number of A nucleotides: $A. This most often results from variants ".
					"which are assigned to the wrong site in the SNP Report. Results at this site ".
					"are unreliable; A count set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative number of A nucleotides: $A. This most often results from\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## are unreliable; A count set to 0; proceed with caution.\n";
			}
		}
		if($C < 0) {
			$hh_nc_position_info{$curr_spot}->{C} = 0;
			
			if($C < -0.1) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative number of C nucleotides: $C. This most often results from variants ".
					"which are assigned to the wrong site in the SNP Report. Results at this site ".
					"are unreliable; C count set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative number of C nucleotides: $C. This most often results from\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## are unreliable; C count set to 0; proceed with caution.\n";
			}
		}
		if($G < 0) {
			$hh_nc_position_info{$curr_spot}->{G} = 0;
			
			if($G < -0.1) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative number of G nucleotides: $G. This most often results from variants ".
					"which are assigned to the wrong site in the SNP Report. Results at this site ".
					"are unreliable; G count set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative number of G nucleotides: $G. This most often results from\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## are unreliable; G count set to 0; proceed with caution.\n";
			}
		}
		if($T < 0) {
			$hh_nc_position_info{$curr_spot}->{T} = 0;
			
			if($T < -0.1) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative number of T nucleotides: $T. This most often results from variants ".
					"which are assigned to the wrong site in the SNP Report. Results at this site ".
					"are unreliable; T count set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative number of T nucleotides: $T. This most often results from\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## are unreliable; T count set to 0; proceed with caution.\n";
			}
		}
		
		# Do the same warning and correction for the proportion data
		my $A_prop = $hh_nc_position_info{$curr_spot}->{A_prop};
		my $C_prop = $hh_nc_position_info{$curr_spot}->{C_prop};
		my $G_prop = $hh_nc_position_info{$curr_spot}->{G_prop};
		my $T_prop = $hh_nc_position_info{$curr_spot}->{T_prop};
		
		if($A_prop < 0) {
			$hh_nc_position_info{$curr_spot}->{A_prop} = 0;
			
			if($A_prop < -0.001) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative proportion of A nucleotides of $A_prop. This may result from rounding error, ".
					"in which case the number will be very small in magnitude, or variants ".
					"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
					"are unreliable. In either case, A prop has set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative proportion of A nucleotides: $A_prop.\n".
					"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## may be unreliable; A prop set to 0; proceed with caution.\n";
			}
		}
		if($C_prop < 0) {
			$hh_nc_position_info{$curr_spot}->{C_prop} = 0;
			
			if($C_prop < -0.001) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative proportion of C nucleotides of $C_prop. This may result from rounding error, ".
					"in which case the number will be very small in magnitude, or variants ".
					"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
					"are unreliable. In either case, C prop has set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative proportion of C nucleotides: $C_prop.\n".
					"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## may be unreliable; C prop set to 0; proceed with caution.\n";
			}
		}
		if($G_prop < 0) {
			$hh_nc_position_info{$curr_spot}->{G_prop} = 0;
			
			if($G_prop < -0.001) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative proportion of G nucleotides of $G_prop. This may result from rounding error, ".
					"in which case the number will be very small in magnitude, or variants ".
					"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
					"are unreliable. In either case, G prop has set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative proportion of G nucleotides: $G_prop.\n".
					"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## may be unreliable; G prop set to 0; proceed with caution.\n";
			}
		}
		if($T_prop < 0) {
			$hh_nc_position_info{$curr_spot}->{T_prop} = 0;
			
			if($T_prop < -0.001) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\tVariant data at this site ".
					"imply a negative proportion of T nucleotides of $T_prop. This may result from rounding error, ".
					"in which case the number will be very small in magnitude, or variants ".
					"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
					"are unreliable. In either case, T prop has set to 0; proceed with caution.\n";
				close ERROR_FILE;
				chdir('..');
				
				warn "\n## WARNING: In $file_nm, the variant at site $curr_spot,\n".
					"## the variant data imply a negative proportion of T nucleotides: $T_prop.\n".
					"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
					"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
					"## may be unreliable; T prop set to 0; proceed with caution.\n";
			}
		}
	} # finish calculating mean coverage for noncoding sites
	
	# EXCLUDE sites with 0 actual variants, or else exclude variants BELOW MAF
	my $variants_excluded = 0;
	#my %variants_excluded;
	foreach my $curr_spot (@stored_positions_sorted) {				
		my $cov = $hh_nc_position_info{$curr_spot}->{cov};
		my $min_COUNT = ($minfreq * $cov);
		
		my $A = $hh_nc_position_info{$curr_spot}->{A};
		my $C = $hh_nc_position_info{$curr_spot}->{C};
		my $G = $hh_nc_position_info{$curr_spot}->{G};
		my $T = $hh_nc_position_info{$curr_spot}->{T};
		
		my $A_prop = $hh_nc_position_info{$curr_spot}->{A_prop};
		my $C_prop = $hh_nc_position_info{$curr_spot}->{C_prop};
		my $G_prop = $hh_nc_position_info{$curr_spot}->{G_prop};
		my $T_prop = $hh_nc_position_info{$curr_spot}->{T_prop};
		
		my $majority_nucleotide = $hh_nc_position_info{$curr_spot}->{maj_nt};
		my $majority_prop_key = "$majority_nucleotide\_prop";
		my $curr_majority_count = $hh_nc_position_info{$curr_spot}->{$majority_nucleotide};
		
		#print "\n\nFOR THE MIN MAF TEST: $curr_spot\n".
		#	"\nMajority nucleotide: $majority_nucleotide\n".
		#	"Count: $curr_majority_count\nCov: $cov\nMin count: $min_COUNT\n".
		#	"Min freq: $minfreq\nA: $A\nC: $C\nG: $G\nT: $T\n".
		#	"A prop: $A_prop\nC prop: $C_prop\nG prop: $G_prop\nT prop: $T_prop\n";
		
		my $leftover_prop;
		
		if($curr_majority_count > $min_COUNT) { # Make sure there are site data
			if(($A_prop > 0) && ($A_prop < $minfreq)) {
				$variants_excluded++;
				$hh_nc_position_info{$curr_spot}->{cov} -= $A;
				$hh_nc_position_info{$curr_spot}->{A} = 0;
				$A = 0;
				$hh_nc_position_info{$curr_spot}->{$majority_prop_key} += $A_prop;
				$hh_nc_position_info{$curr_spot}->{A_prop} = 0;
				$leftover_prop += $A_prop;
				$A_prop = 0;
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\t".
					"Variant 'A' excluded from analysis because it falls below the minimum ".
						"minor allele frequency\n";
				close ERROR_FILE;
				chdir('..');
				
				print "\n## Variant 'A' excluded from analysis because it falls below the\n".
					"## minimum minor allele frequency at:\n".
					"## $file_nm|$curr_spot\n";
			} 
			if(($C_prop > 0) && ($C_prop < $minfreq)) {
				$variants_excluded++;
				$hh_nc_position_info{$curr_spot}->{cov} -= $C;
				$hh_nc_position_info{$curr_spot}->{C} = 0;
				$C = 0;
				$hh_nc_position_info{$curr_spot}->{$majority_prop_key} += $C_prop;
				$hh_nc_position_info{$curr_spot}->{C_prop} = 0;
				$leftover_prop += $C_prop;
				$C_prop = 0;
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\t".
					"Variant 'C' excluded from analysis because it falls below the minimum ".
						"minor allele frequency\n";
				close ERROR_FILE;
				chdir('..');
				
				print "\n## Variant 'C' excluded from analysis because it falls below the\n".
					"## minimum minor allele frequency at:\n".
					"## $file_nm|$curr_spot\n";
			} 
			if(($G_prop > 0) && ($G_prop < $minfreq)) {
				$variants_excluded++;
				$hh_nc_position_info{$curr_spot}->{cov} -= $G;
				$hh_nc_position_info{$curr_spot}->{G} = 0;
				$G = 0;
				$hh_nc_position_info{$curr_spot}->{$majority_prop_key} += $G_prop;
				$hh_nc_position_info{$curr_spot}->{G_prop} = 0;
				$leftover_prop += $G_prop;
				$G_prop = 0;
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\t".
					"Variant 'G' excluded from analysis because it falls below the minimum ".
						"minor allele frequency\n";
				close ERROR_FILE;
				chdir('..');
				
				print "\n## Variant 'G' excluded from analysis because it falls below the\n".
					"## minimum minor allele frequency at:\n".
					"## $file_nm|$curr_spot\n";
			} 
			if(($T_prop > 0) && ($T_prop < $minfreq)) {
				$variants_excluded++;
				$hh_nc_position_info{$curr_spot}->{cov} -= $T;
				$hh_nc_position_info{$curr_spot}->{T} = 0;
				$T = 0;
				$hh_nc_position_info{$curr_spot}->{$majority_prop_key} += $T_prop;
				$hh_nc_position_info{$curr_spot}->{T_prop} = 0;
				$leftover_prop += $T_prop;
				$T_prop = 0;
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				print ERROR_FILE "$file_nm\tNA\t$curr_spot\t".
					"Variant 'T' excluded from analysis because it falls below the minimum ".
						"minor allele frequency\n";
				close ERROR_FILE;
				chdir('..');
				
				print "\n## Variant 'T' excluded from analysis because it falls below the\n".
					"## minimum minor allele frequency at:\n".
					"## $file_nm|$curr_spot\n";
			}
		} # sites are further normalized later
		
		if($majority_prop_key eq 'A_prop') {
			$A_prop += $leftover_prop;
		} elsif($majority_prop_key eq 'C_prop') {
			$C_prop += $leftover_prop;
		} elsif($majority_prop_key eq 'G_prop') {
			$G_prop += $leftover_prop;
		} elsif($majority_prop_key eq 'T_prop') {
			$T_prop += $leftover_prop;
		}
		
		my $nt_sum = ($A + $C + $G + $T);
		my $nt_prop_sum = ($A_prop + $C_prop + $G_prop + $T_prop);
		my $updated_cov = $hh_nc_position_info{$curr_spot}->{cov};
		
		#print "New A: $A\nNew C: $C\nNew G: $G\nNew T: $T\n".
		#	"New A prop: $A_prop\nNew C prop: $C_prop\nNew G prop: $G_prop\nNew T prop: $T_prop\n";
		
		# DELETE the hash element if there are no longer any variants here.
		if(($A + $C + $G == 0) || ($A + $C + $T == 0) || ($A + $G + $T == 0) ||
			($C + $G + $T == 0)) {
			delete($hh_nc_position_info{$curr_spot});
		}
		
		#print "New nt sum: $nt_sum\nNew cov: $updated_cov\nNew prop sum: $nt_prop_sum\n\n";
		
		# Round for comparisons
		my $rounded_nt_sum = sprintf("%.3f",$nt_sum);
		my $rounded_updated_cov = sprintf("%.3f",$updated_cov);
		
		# ERRORS
		#if($rounded_nt_sum != $rounded_updated_cov) {
		#	die "\n## WARNING: The MIN. M.A.F. or conflicting coverages have caused an error, ".
		#		"such that the new coverage \n". 
		#		"## does not equal the nucleotide sum. Please contact the author; ".
		#		"SNPGenie terminated.\n\n";
		#}
		
		#if($nt_prop_sum != 1) {
		#	die "\n## WARNING: The MIN. M.A.F. or conflicting coverages have caused an error, ".
		#		"such that the new nucleotide proportion sum does not equal 1.".
		#		"\n## Please contact the author; ".
		#		"SNPGenie terminated.\n\n";
		#}
	}
	
	if($variants_excluded > 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$file_nm\tNA\tN/A\t".
			"A total of $variants_excluded variants have been excluded because they fall below the ".
				"minimum minor allele frequency\n";
		close ERROR_FILE;
		chdir('..');
						
		print "\n## In $file_nm|N/A\n".
			"## A total of $variants_excluded variants have been excluded because they\n".
			"## fall below the minimum minor allele frequency.\n";
	}

	# LABEL AS CODING IF WITHIN A CODING SEGMENT (from GTF)
	# Go through each product, through each segment
	foreach my $curr_product (@curr_products) {
		#print "\nProduct: $curr_product, ";

		# SAVE A num_segments in the hh!
		my $num_segments = $hh_product_position_info{$curr_product}->{num_segments}; 
		#print "\nFor product $curr_product, the num_segments is: $num_segments\n";
		
		# For each segment, step through each site, label coding
		for(my $i = 1; $i <= $num_segments; $i++) {
			my $this_start_key = 'start_' . $i;
			my $this_stop_key = 'stop_' . $i;
			
			my $start = $hh_product_position_info{$curr_product}->{$this_start_key};
			my $stop = $hh_product_position_info{$curr_product}->{$this_stop_key};
			
			for(my $j = $start; $j <= $stop; $j++) { # FOR EACH SITE IN SEGMENT 1 +=1
				$hh_nc_position_info{$j}->{coding} += 1;
			}
		}
	}
	
	# LABEL AS CODING IF CODING from the REVCOM '-' STRAND (from GTF)
	# Go through each product, through each segment
	if($complementmode) {
		foreach my $revcom_product (@curr_compl_products_ordered_by_start) {
			#print "\nProduct: $curr_product, ";

			# SAVE A num_segments in the hh!
			my $num_segments = $hh_compl_position_info{$revcom_product}->{num_segments};
			#print "\nFor revcom product $revcom_product, the num_segments is: $num_segments\n";
			
			# For each segment, step through each site, label coding
			for(my $i = 1; $i <= $num_segments; $i++) {
				my $this_start_key = 'start_' . $i;
				my $this_stop_key = 'stop_' . $i;
				
				my $start = $hh_compl_position_info{$revcom_product}->{$this_start_key};
				my $stop = $hh_compl_position_info{$revcom_product}->{$this_stop_key};
				
				for(my $j = $start; $j <= $stop; $j++) { # FOR EACH SITE IN SEGMENT 1 +=1
					$hh_nc_position_info{$j}->{coding} += 1;
				}
			}
		}
	}
	
	# CALCULATE PI, etc. for NONCODING DNA (i.e., without regarding to product)
	# Variables for calculating averages
	my $sum_mean_pw_diffs = 0;
	my $sum_mean_pw_diffs_coding = 0;
	my $sum_mean_pw_diffs_nc = 0;
	
	my $num_coding_sites = 0;
	my $num_coding_sites_polymorphic = 0;
	my $num_nc_sites = 0;
	my $num_nc_sites_polymorphic = 0;
	
	my $sum_gdiv;
	my $sum_gdiv_coding_poly;
	my $sum_gdiv_nc_poly;
	
	print "COMPLETED.\n";
	
	#comeback
	print "\nProcessing all individual sites (nucleotides take time)... ";
	#print "\nSeq. length is: ".length($seq)."\n";
	for (my $i = 0; $i < length($seq); $i++) { # FOR EACH SITE IN FASTA
		my $position = $i+1;
		my $ref_nt = substr($seq,$i,1);
		
		#print "Site $position\t";
		
		# If polymorphic; we've saved this while going through the SNP report
		if($hh_nc_position_info{$position}->{polymorphic} == 1) { # poly
			# We have already DELETED records that contained no real variants, e.g.,
			# all variant had frequencies 0% or were below the MAF
			#print "\nSite $position is polymorphic\n"; # DOUBLE-CHECK: PASS
			
			my $A = $hh_nc_position_info{$position}->{A};
			my $C = $hh_nc_position_info{$position}->{C};
			my $G = $hh_nc_position_info{$position}->{G};
			my $T = $hh_nc_position_info{$position}->{T};
			my $maj_nt = $hh_nc_position_info{$position}->{maj_nt};
			my $cov = $A + $C + $G + $T;
			my $cov_old = $hh_nc_position_info{$position}->{cov};
			
			my $cov_for_comp = sprintf("%.3f",$cov);
			my $cov_old_for_comp = sprintf("%.3f",$cov_old);
			
			# WARN if the coverage doesn't match the sum of nucleotide counts
			if($cov_for_comp != $cov_old_for_comp) {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\tNA\t$position\tCoverage ($cov_old_for_comp) does not equal ".
					"the nucleotide sum ($cov_for_comp).\n";
				close ERROR_FILE;
				chdir('..');
				
				print "\n## WARNING: In $file_nm, at site $position,\n".
					"## the coverage ($cov_old_for_comp) does not equal the nucleotide sum ($cov_for_comp).\n";
					#."\n\tSITE $position\nA $A\nC $C\nG $G\nT $T\nCOV $cov_old\nSUM $cov\n";
			}
			
			my $A_prop = ($A / $cov);
			my $C_prop = ($C / $cov);
			my $G_prop = ($G / $cov);
			my $T_prop = ($T / $cov);
			
	
			# METAPOPULATION
			if($file_nm =~ /B/) { # BIRD
				#print "\nIt's a bird and we've got A=$A_prop C=$C_prop G=$G_prop T=$T_prop\n";
				push(@{$master_frequencies_hh{$position}->{BIRD}->{A_props_arr}}, $A_prop);
				push(@{$master_frequencies_hh{$position}->{BIRD}->{C_props_arr}}, $C_prop);
				push(@{$master_frequencies_hh{$position}->{BIRD}->{G_props_arr}}, $G_prop);
				push(@{$master_frequencies_hh{$position}->{BIRD}->{T_props_arr}}, $T_prop);
				
##				if($position == 4543) {
##					print "\nFreq of BIRD A at nonpoly $position is $A_prop\n";
##					print "\nFreq of BIRD C at nonpoly $position is $C_prop\n";
##					print "\nFreq of BIRD G at nonpoly $position is $G_prop\n";
##					print "\nFreq of BIRD T at nonpoly $position is $T_prop\n";
##				}
				
			} else { # MOSQ
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{A_props_arr}}, $A_prop);
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{C_props_arr}}, $C_prop);
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{G_props_arr}}, $G_prop);
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{T_props_arr}}, $T_prop);
				
##				if($position == 4543) {
##					print "\nFreq of MOSQ A at nonpoly $position is $A_prop\n";
##					print "\nFreq of MOSQ C at nonpoly $position is $C_prop\n";
##					print "\nFreq of MOSQ G at nonpoly $position is $G_prop\n";
##					print "\nFreq of MOSQ T at nonpoly $position is $T_prop\n";
##				}
			}
			
			# FST purposes
			if($file_nm =~ /^([\w\.]+?)_/) {
				my $sample_name = $1;
				$master_frequencies_hh{$position}->{$sample_name}->{A_freq} = $A_prop;
				$master_frequencies_hh{$position}->{$sample_name}->{C_freq} = $C_prop;
				$master_frequencies_hh{$position}->{$sample_name}->{G_freq} = $G_prop;
				$master_frequencies_hh{$position}->{$sample_name}->{T_freq} = $T_prop;
			}
			
			
			my $gdiv = (1 - ($A_prop * $A_prop) - ($C_prop * $C_prop) - 
				($G_prop * $G_prop) - ($T_prop * $T_prop));
			#print "$position\t$gdiv\n";
			$hh_nc_position_info{$position}->{gdiv} = $gdiv;
			$sum_gdiv += $gdiv;
			
			my $num_pw_diffs = ($A*$C + $A*$G + $A*$T + $C*$G + $C*$T + $G*$T);
			my $total_pw_comps = (($cov * $cov) - $cov)/2;
			my $mean_pw_diffs = ($num_pw_diffs / $total_pw_comps); # this IS pi, div by 1
			$hh_nc_position_info{$position}->{pi} = $mean_pw_diffs;
			$sum_mean_pw_diffs += $mean_pw_diffs;
			
			if($hh_nc_position_info{$position}->{coding} > 0) { # poly-coding site
				$num_coding_sites ++;
				$num_coding_sites_polymorphic ++;
				$sum_mean_pw_diffs_coding += $mean_pw_diffs;
				$sum_gdiv_coding_poly += $gdiv;
				#print "CODING\n";
			} else { # poly-noncoding site
				$num_nc_sites ++;
				$num_nc_sites_polymorphic ++;
				$sum_mean_pw_diffs_nc += $mean_pw_diffs;
				$sum_gdiv_nc_poly += $gdiv;
				#print "NC\n";
			}
			
			#print "\nMy pi at site $position is $mean_pw_diffs\n";

		} else { # nonpoly
			# Save the reference nucleotide
			#print "\nSite $position is not polymorphic\n";
			
			my $A_prop;
			my $C_prop;
			my $G_prop;
			my $T_prop;
			
			if($ref_nt eq 'A') {
				$A_prop = 1;
				$C_prop = 0;
				$G_prop = 0;
				$T_prop = 0;
			} elsif($ref_nt eq 'C') {
				$A_prop = 0;
				$C_prop = 1;
				$G_prop = 0;
				$T_prop = 0;
			} elsif($ref_nt eq 'G') {
				$A_prop = 0;
				$C_prop = 0;
				$G_prop = 1;
				$T_prop = 0;
			} elsif($ref_nt eq 'T') {
				$A_prop = 0;
				$C_prop = 0;
				$G_prop = 0;
				$T_prop = 1;
			} else {
				#warn "\n\nThere's an N in the reference sequence!\n\n";
			}
			
			
			# METAPOPULATION
			if($file_nm =~ /B/) { # BIRD
				push(@{$master_frequencies_hh{$position}->{BIRD}->{A_props_arr}}, $A_prop);
				push(@{$master_frequencies_hh{$position}->{BIRD}->{C_props_arr}}, $C_prop);
				push(@{$master_frequencies_hh{$position}->{BIRD}->{G_props_arr}}, $G_prop);
				push(@{$master_frequencies_hh{$position}->{BIRD}->{T_props_arr}}, $T_prop);
				
##				if($position == 4543) {
##					print "\nFreq of BIRD A at nonpoly $position is $A_prop\n";
##					print "\nFreq of BIRD C at nonpoly $position is $C_prop\n";
##					print "\nFreq of BIRD G at nonpoly $position is $G_prop\n";
##					print "\nFreq of BIRD T at nonpoly $position is $T_prop\n";
##				}
				
				
				
			} else { # MOSQ
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{A_props_arr}}, $A_prop);
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{C_props_arr}}, $C_prop);
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{G_props_arr}}, $G_prop);
				push(@{$master_frequencies_hh{$position}->{MOSQ}->{T_props_arr}}, $T_prop);
				
##				if($position == 4543) {
##					print "\nFreq of MOSQ A at nonpoly $position is $A_prop\n";
##					print "\nFreq of MOSQ C at nonpoly $position is $C_prop\n";
##					print "\nFreq of MOSQ G at nonpoly $position is $G_prop\n";
##					print "\nFreq of MOSQ T at nonpoly $position is $T_prop\n";
##				}
			}
			
			# FST purposes
			if($file_nm =~ /^([\w\.]+?)_/) {
				my $sample_name = $1;
				$master_frequencies_hh{$position}->{$sample_name}->{A_freq} = $A_prop;
				$master_frequencies_hh{$position}->{$sample_name}->{C_freq} = $C_prop;
				$master_frequencies_hh{$position}->{$sample_name}->{G_freq} = $G_prop;
				$master_frequencies_hh{$position}->{$sample_name}->{T_freq} = $T_prop;
			}
			
			
			if($hh_nc_position_info{$position}->{coding} > 0) { # nonpoly-coding site
				$num_coding_sites ++;
				#print "CODING\n";
			} else { # nonpoly-noncoding site
				$num_nc_sites ++;
				#print "NC\n";
			}	
		}
	}

	# Calculate total pi values
	my $num_sites_total = ($num_coding_sites + $num_nc_sites);
	
	my $pi_total;
	if($num_sites_total > 0) {
		$pi_total = ($sum_mean_pw_diffs / $num_sites_total);
	} else {
		$pi_total = '*';
	}
	
	my $pi_coding;
	if($num_coding_sites > 0) {
		$pi_coding = ($sum_mean_pw_diffs_coding / $num_coding_sites);
	} else {
		$pi_coding = '*';
	}
	
	my $pi_nc;
	if($num_nc_sites > 0) {
		$pi_nc = ($sum_mean_pw_diffs_nc / $num_nc_sites);
	} else {
		$pi_nc = '*';
	}
	
	my $mean_gdiv;
	if($num_sites_total > 0) {
		$mean_gdiv = ($sum_gdiv / $num_sites_total);
	} else {
		$mean_gdiv = '*';
	}
	
	my $mean_gdiv_coding_poly;
	if($num_coding_sites_polymorphic > 0) {
		$mean_gdiv_coding_poly = ($sum_gdiv_coding_poly / $num_coding_sites_polymorphic);
	} else {
		$mean_gdiv_coding_poly = '*';
	}
	
	my $mean_gdiv_nc_poly;
	if($num_nc_sites_polymorphic > 0) {
		$mean_gdiv_nc_poly = ($sum_gdiv_nc_poly / $num_nc_sites_polymorphic);
	} else {
		$mean_gdiv_nc_poly = '*';
	}
	
	$h_nc_results{num_poly_coding_sites} = $num_coding_sites_polymorphic;
	$h_nc_results{num_poly_nc_sites} = $num_nc_sites_polymorphic;
	
	$h_nc_results{mean_gdiv} = $mean_gdiv;
	$h_nc_results{mean_gdiv_coding_poly} = $mean_gdiv_coding_poly;
	$h_nc_results{mean_gdiv_nc_poly} = $mean_gdiv_nc_poly;
	
	$h_nc_results{num_coding_sites} = $num_coding_sites;
	$h_nc_results{num_nc_sites} = $num_nc_sites;
	$h_nc_results{num_sites_total} = $num_sites_total;
	$h_nc_results{sum_mean_pw_diffs} = $sum_mean_pw_diffs;
	$h_nc_results{sum_mean_pw_diffs_coding} = $sum_mean_pw_diffs_coding;
	$h_nc_results{sum_mean_pw_diffs_nc} = $sum_mean_pw_diffs_nc;
	$h_nc_results{pi_total} = $pi_total;
	$h_nc_results{pi_coding} = $pi_coding;
	$h_nc_results{pi_nc} = $pi_nc;
	
	#print "\nMy total pi is: $pi_total\nMy coding pi is: $pi_coding\n".
	#	"My noncoding sum diffs is: $sum_mean_pw_diffs_nc\n".
	#	"My noncoding sites is: $num_nc_sites\nMy noncoding pi is: $pi_nc\n".
	#	"My coding sites is: $num_coding_sites\n";
	
	$pop_summary_line .= "$file_nm\t$num_sites_total\t$num_coding_sites\t".
		"$num_nc_sites\t".
		"$pi_total\t$pi_coding\t$pi_nc\t";
	
	print "COMPLETED.\n";
	
	# NON-CODING WORK HAPPENS HERE
	# We must assign EVERY SITE in the FASTA, labeling as polymorphic=1/0 and as
	# coding=1/0. We will determine that it is polymorphic if there has been a variant 
	# stored at the site. We will determine it is coding simply if it falls within the
	# range of some product.
	
	my $pop_sum_ndiffs = 0;
	my $pop_sum_sdiffs = 0;
	my $pop_sum_ndiffs_vRef = 0;
	my $pop_sum_sdiffs_vRef = 0;
	my $pop_sum_nsites = 0;
	my $pop_sum_ssites = 0;
	my $pop_sum_nsites_Ref = 0;
	my $pop_sum_ssites_Ref = 0;
	my $pop_sum_gdiv = 0;
	my $pop_sum_poly_sites = 0;
	my $pop_sum_Nsites_gdiv = 0;
	my $pop_sum_Ssites_gdiv = 0;
	my $pop_sum_Asites_gdiv = 0;
	my $pop_sum_N_gdiv = 0;
	my $pop_sum_S_gdiv = 0;
	my $pop_sum_A_gdiv = 0;
	
	print "\nProcessing population genetic estimates codon-by-codon (beware stochasticity!)... ";
	# Here we loop through each product present in this SNP Report
	# This is where the major CODON work takes place
	foreach my $curr_product (@curr_products) {
		
		# CALCULATE AVERAGE COVERAGE at each of the sites in this product.
		# We are REPLACING the coverage currently stored with this.
		my @stored_positions_sorted = sort {$a <=> $b} (keys %{$hh_product_position_info{$curr_product}});
		foreach my $curr_spot (@stored_positions_sorted) {
			#print "\ncurr_product is: $curr_product; curr_spot is: $curr_spot\n";
			if (!($curr_spot =~ 'start') && !($curr_spot =~ 'stop') && 
				!($curr_spot =~ 'fasta') && !($curr_spot =~ 'num_segments')) {
				#print "\n@{$hh_product_position_info{$curr_product}->{$curr_spot}->{cov_arr}} is: ".@{$hh_product_position_info{$curr_product}->{$curr_spot}->{cov_arr}}."\n\n";
				#if(@{$hh_product_position_info{$curr_product}->{$curr_spot}->{cov_arr}}) {
				my $cov_sum = 0;
				my $cov_denom = 0;
				
				#print "\n curr spot is: $curr_spot\n";
				
				foreach my $cov_mm (@{$hh_product_position_info{$curr_product}->{$curr_spot}->{cov_arr}}) {
					$cov_sum += $cov_mm;
					$cov_denom += 1;
				}
				my $this_avg_cov = ($cov_sum / $cov_denom);
				$hh_product_position_info{$curr_product}->{$curr_spot}->{cov} = $this_avg_cov;
				
				# DETERMINE THE MAJORITY NUCLEOTIDE, TOO
				my $A = $hh_product_position_info{$curr_product}->{$curr_spot}->{A};
				my $C = $hh_product_position_info{$curr_product}->{$curr_spot}->{C};
				my $G = $hh_product_position_info{$curr_product}->{$curr_spot}->{G};
				my $T = $hh_product_position_info{$curr_product}->{$curr_spot}->{T};
				
				my $majority_nucleotide; # a variant nucleotide may have fixed
				my $curr_majority_count = 0;
				if($A > $curr_majority_count) {
					$curr_majority_count = $A;
					$majority_nucleotide = 'A';
				}
				if($C > $curr_majority_count) {
					$curr_majority_count = $C;
					$majority_nucleotide = 'C';
				} 
				if($G > $curr_majority_count) {
					$curr_majority_count = $G;
					$majority_nucleotide = 'G';
				} 
				if($T > $curr_majority_count) {
					$curr_majority_count = $T;
					$majority_nucleotide = 'T';
				}
				
				$hh_product_position_info{$curr_product}->{$curr_spot}->{maj_nt} = $majority_nucleotide;
				
				# WARNING if there is a negative number of sites because of conflicting
				# coverages; round negative nucleotide counts to 0
				if($A < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{A} = 0;
					
					if($A < -0.1) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative number of A nucleotides: $A. This most often results from variants ".
							"which are assigned to the wrong site in the SNP Report. Results at this site ".
							"are unreliable; A count set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative number of A nucleotides: $A. This most often results from\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## are unreliable; A count set to 0; proceed with caution.\n";
					}
				}
				if($C < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{C} = 0;
					
					if($C < -0.1) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative number of C nucleotides: $C. This most often results from variants ".
							"which are assigned to the wrong site in the SNP Report. Results at this site ".
							"are unreliable; C count set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative number of C nucleotides: $C. This most often results from\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## are unreliable; C count set to 0; proceed with caution.\n";
					}
				}
				if($G < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{G} = 0;
					
					if($G < -0.1) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative number of G nucleotides: $G. This most often results from variants ".
							"which are assigned to the wrong site in the SNP Report. Results at this site ".
							"are unreliable; G count set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative number of G nucleotides: $G. This most often results from\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## are unreliable; G count set to 0; proceed with caution.\n";
					}
				}
				if($T < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{T} = 0;
					
					if($T < -0.1) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative number of T nucleotides: $T. This most often results from variants ".
							"which are assigned to the wrong site in the SNP Report. Results at this site ".
							"are unreliable; T count set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative number of T nucleotides: $T. This most often results from\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## are unreliable; T count set to 0; proceed with caution.\n";
					}
				}
				
				# Do the same warning and correction for the proportion data
				my $A_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{A_prop};
				my $C_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{C_prop};
				my $G_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{G_prop};
				my $T_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{T_prop};
				
				if($A_prop < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{A_prop} = 0;
					
					if($A_prop < -0.001) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative proportion of A nucleotides of $A_prop. This may result from rounding error, ".
							"in which case the number will be very small in magnitude, or variants ".
							"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
							"are unreliable. In either case, A prop has set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative proportion of A nucleotides: $A_prop.\n".
							"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## may be unreliable; A prop set to 0; proceed with caution.\n";
					}
				}
				if($C_prop < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{C_prop} = 0;
					
					if($C_prop < -0.001) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative proportion of C nucleotides of $C_prop. This may result from rounding error, ".
							"in which case the number will be very small in magnitude, or variants ".
							"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
							"are unreliable. In either case, C prop has set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative proportion of C nucleotides: $C_prop.\n".
							"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## may be unreliable; C prop set to 0; proceed with caution.\n";
					}
				}
				if($G_prop < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{G_prop} = 0;
					
					if($G_prop < -0.001) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative proportion of G nucleotides of $G_prop. This may result from rounding error, ".
							"in which case the number will be very small in magnitude, or variants ".
							"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
							"are unreliable. In either case, G prop has set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative proportion of G nucleotides: $G_prop.\n".
							"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## may be unreliable; G prop set to 0; proceed with caution.\n";
					}
				}
				if($T_prop < 0) {
					$hh_product_position_info{$curr_product}->{$curr_spot}->{T_prop} = 0;
					
					if($T_prop < -0.001) {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\tVariant data at this site ".
							"imply a negative proportion of T nucleotides of $T_prop. This may result from rounding error, ".
							"in which case the number will be very small in magnitude, or variants ".
							"which are assigned to the wrong site in the SNP Report. If the latter, results at this site ".
							"are unreliable. In either case, T prop has set to 0; proceed with caution.\n";
						close ERROR_FILE;
						chdir('..');
						
						warn "\n## WARNING: In $file_nm, $curr_product, the codon at site $curr_spot,\n".
							"## the variant data imply a negative proportion of T nucleotides: $T_prop.\n".
							"## This may result from rounding error, in which case the number will be very small in magnitude, or\n".
							"## variants which are assigned to the wrong site in the SNP Report. Results at this site\n".
							"## may be unreliable; T prop set to 0; proceed with caution.\n";
					}
				}
			}
		}
		
		# FILTER by MIN. MINOR ALLELE FREQ., if there is one.
		#
		# We USED TO rely on defining nucleotide counts as we went along, then calculate proportions;
		# this was done as: 
		#		my $A_prop = ($num_1_A / $site_cov)
		#
		# NOW we are relying on proportions as we go along instead, then calculating more accurate counts;
		# this is done as: 
		#		$prop_1_A = $hh_product_position_info{$curr_product}->{$site_pos_1}->{A_prop};
		#		my $A_prop = $prop_1_A;
		#		$num_1_A = ($A_prop * $site_1_cov);
		#
		# I believe we should redefine BOTH to ensure future flexibility. Thus, we will
		# ELIMINATE any nucleotides with frequencies falling below the cutoff (i.e., set 
		# both their frequency and count to 0), SUBTRACT them from the coverage (i.e.,
		# subtract the nucleotide count), and keep the majority count; also ADD
		# the nucleotide frequency from the reference frequency
		#
		# This will include modifying:
		#		$hh_product_position_info{$product_name}->{$position5}->{A} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{C} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{G} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{T} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{A_prop} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{C_prop} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{G_prop} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{T_prop} = 0;
		#		$hh_product_position_info{$product_name}->{$position5}->{cov}

		my $variants_excluded = 0;
		#my %variants_excluded;
		foreach my $curr_spot (@stored_positions_sorted) {
			if (!($curr_spot =~ 'start') && !($curr_spot =~ 'stop') && 
				!($curr_spot =~ 'fasta') && !($curr_spot =~ 'num_segments')) {
				
				my $cov = $hh_product_position_info{$curr_product}->{$curr_spot}->{cov};
				my $min_COUNT = ($minfreq * $cov);
				
				my $A = $hh_product_position_info{$curr_product}->{$curr_spot}->{A};
				my $C = $hh_product_position_info{$curr_product}->{$curr_spot}->{C};
				my $G = $hh_product_position_info{$curr_product}->{$curr_spot}->{G};
				my $T = $hh_product_position_info{$curr_product}->{$curr_spot}->{T};
				
				my $A_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{A_prop};
				my $C_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{C_prop};
				my $G_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{G_prop};
				my $T_prop = $hh_product_position_info{$curr_product}->{$curr_spot}->{T_prop};
				
				my $majority_nucleotide = $hh_product_position_info{$curr_product}->{$curr_spot}->{maj_nt};
				my $majority_prop_key = "$majority_nucleotide\_prop";
				my $curr_majority_count = $hh_product_position_info{$curr_product}->{$curr_spot}->{$majority_nucleotide};
				
				#print "\n\nFOR THE MIN MAF TEST: $curr_product $curr_spot\n".
				#	"\nMajority nucleotide: $majority_nucleotide\n".
				#	"Count: $curr_majority_count\nCov: $cov\nMin count: $min_COUNT\n".
				#	"Min freq: $minfreq\nA: $A\nC: $C\nG: $G\nT: $T\n".
				#	"A prop: $A_prop\nC prop: $C_prop\nG prop: $G_prop\nT prop: $T_prop\n";
				
				my $leftover_prop;
				
				if($curr_majority_count > $min_COUNT) { # Make sure there are site data
					if(($A_prop > 0) && ($A_prop < $minfreq)) {
						$variants_excluded++;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{cov} -= $A;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{A} = 0;
						$A = 0;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{$majority_prop_key} += $A_prop;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{A_prop} = 0;
						$leftover_prop += $A_prop;
						$A_prop = 0;
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\t".
							"Variant 'A' excluded from analysis because it falls below the minimum ".
								"minor allele frequency\n";
						close ERROR_FILE;
						chdir('..');
						
						print "\n## Variant 'A' excluded from analysis because it falls below the\n".
							"## minimum minor allele frequency at:\n".
							"## $file_nm|$curr_product|$curr_spot\n";
					} 
					if(($C_prop > 0) && ($C_prop < $minfreq)) {
						$variants_excluded++;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{cov} -= $C;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{C} = 0;
						$C = 0;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{$majority_prop_key} += $C_prop;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{C_prop} = 0;
						$leftover_prop += $C_prop;
						$C_prop = 0;
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\t".
							"Variant 'C' excluded from analysis because it falls below the minimum ".
								"minor allele frequency\n";
						close ERROR_FILE;
						chdir('..');
						
						print "\n## Variant 'C' excluded from analysis because it falls below the\n".
							"## minimum minor allele frequency at:\n".
							"## $file_nm|$curr_product|$curr_spot\n";
					} 
					if(($G_prop > 0) && ($G_prop < $minfreq)) {
						$variants_excluded++;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{cov} -= $G;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{G} = 0;
						$G = 0;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{$majority_prop_key} += $G_prop;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{G_prop} = 0;
						$leftover_prop += $G_prop;
						$G_prop = 0;
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\t".
							"Variant 'G' excluded from analysis because it falls below the minimum ".
								"minor allele frequency\n";
						close ERROR_FILE;
						chdir('..');
						
						print "\n## Variant 'G' excluded from analysis because it falls below the\n".
							"## minimum minor allele frequency at:\n".
							"## $file_nm|$curr_product|$curr_spot\n";
					} 
					if(($T_prop > 0) && ($T_prop < $minfreq)) {
						$variants_excluded++;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{cov} -= $T;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{T} = 0;
						$T = 0;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{$majority_prop_key} += $T_prop;
						$hh_product_position_info{$curr_product}->{$curr_spot}->{T_prop} = 0;
						$leftover_prop += $T_prop;
						$T_prop = 0;
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						print ERROR_FILE "$file_nm\t$curr_product\t$curr_spot\t".
							"Variant 'T' excluded from analysis because it falls below the minimum ".
								"minor allele frequency\n";
						close ERROR_FILE;
						chdir('..');
						
						print "\n## Variant 'T' excluded from analysis because it falls below the\n".
							"## minimum minor allele frequency at:\n".
							"## $file_nm|$curr_product|$curr_spot\n";
					}
				} # sites are further normalized later
				
				if($majority_prop_key eq 'A_prop') {
					$A_prop += $leftover_prop;
				} elsif($majority_prop_key eq 'C_prop') {
					$C_prop += $leftover_prop;
				} elsif($majority_prop_key eq 'G_prop') {
					$G_prop += $leftover_prop;
				} elsif($majority_prop_key eq 'T_prop') {
					$T_prop += $leftover_prop;
				}
				
				my $nt_sum = ($A + $C + $G + $T);
				my $nt_prop_sum = ($A_prop + $C_prop + $G_prop + $T_prop);
				my $updated_cov = $hh_product_position_info{$curr_product}->{$curr_spot}->{cov};
				
				#print "New A: $A\nNew C: $C\nNew G: $G\nNew T: $T\n".
				#	"New A prop: $A_prop\nNew C prop: $C_prop\nNew G prop: $G_prop\nNew T prop: $T_prop\n";
				
				# DELETE the hash element if there are no longer any variants here.
				if(($A + $C + $G == 0) || ($A + $C + $T == 0) || ($A + $G + $T == 0) ||
					($C + $G + $T == 0)) {
					delete($hh_product_position_info{$curr_product}->{$curr_spot});
				}
				
				#print "New nt sum: $nt_sum\nNew cov: $updated_cov\nNew prop sum: $nt_prop_sum\n\n";
				
				# Round for comparisons
				my $rounded_nt_sum = sprintf("%.3f",$nt_sum);
				my $rounded_updated_cov = sprintf("%.3f",$updated_cov);
				
				# ERRORS
				#if($rounded_nt_sum != $rounded_updated_cov) {
				#	die "\n## WARNING: The MIN. M.A.F. or conflicting coverages have caused an error, ".
				#		"such that the new coverage \n". 
				#		"## does not equal the nucleotide sum. Please contact the author; ".
				#		"SNPGenie terminated.\n\n";
				#}
				
				#if($nt_prop_sum != 1) {
				#	die "\n## WARNING: The MIN. M.A.F. or conflicting coverages have caused an error, ".
				#		"such that the new nucleotide proportion sum does not equal 1.".
				#		"\n## Please contact the author; ".
				#		"SNPGenie terminated.\n\n";
				#}
				
			}
		}
		
		if($variants_excluded > 0) {
			chdir('SNPGenie_Results');
			open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
			print ERROR_FILE "$file_nm\t$curr_product\tN/A\t".
				"A total of $variants_excluded variants have been excluded because they fall below the ".
					"minimum minor allele frequency\n";
			close ERROR_FILE;
			chdir('..');
							
			print "\n## In $file_nm|$curr_product|N/A\n".
				"## A total of $variants_excluded variants have been excluded because they\n".
				"## fall below the minimum minor allele frequency.\n";
		}
		
		#print "\n\nProduct $curr_product should correspond to the fasta: $fasta_to_open\n\n";
		
		#print "\n\nFor the fasta $fasta_to_open I get this sequence:\n$seq\n\n";
		
		#if (length($seq) != (1 + $hh_product_position_info{$curr_product}->{stop_1} - 
		#	$hh_product_position_info{$curr_product}->{start_1})) {
		#	print "\nThe length does not equal stop-start+1 in $curr_product\n";
		#}
		
		# The $seq contains the ENTIRE fasta file, not just the product's sites
		
		# WARN if the reference in the FASTA doesn't match the SNP Report
		for (my $i = 0; $i < length($seq); $i++) {
			my $position = $i+1;
			
			if($hh_product_position_info{$curr_product}->{$position}) { # If the variant exists
				#print "YES, there IS a variant stored at $position\n";
				my $snp_report_reference_nt = $hh_product_position_info{$curr_product}->{$position}->{reference};
				my $fasta_nt = $seq_by_index_arr[$position-1];
				
				if($snp_report_reference_nt ne $fasta_nt) {
					chdir('SNPGenie_Results');
					open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
					# FILE | PRODUCT | SITE | WARNING
					print ERROR_FILE "$file_nm\t$curr_product\t$position\tThe reference nucleotide ".
						"in the SNP Report does not match the FASTA file. Results at this site ".
						"are unreliable: proceed with caution, and consider re-calling SNPs (or ".
						"fixing the SNP Reports) to ".
						"address the issue.\n";
					close ERROR_FILE;
					chdir('..');
					
					print "\n## WARNING: In $file_nm, $curr_product, site $position,\n".
						"## the reference nucleotide in the SNP Report does not match the FASTA file.\n".
						"## Results at this site are unreliable: proceed with caution, and consider\n".
						"## re-calling SNPs (or fixing the SNP Reports) to address the issue.\n";
				}
			}
		}
		
		# New segments approach
		my $num_segments = $hh_product_position_info{$curr_product}->{num_segments};
		#print "\nFor product $curr_product, the num_segments is: $num_segments\n";
		my $coding_length_sum;
		my $coding_seq;
		my @coding_seq_pos_arr;
		my %segment_lengths; # to store the length of each of the $num_segments segments
		#my %segment_point_remainders;
		
		# New segments approach
		my %codon_by_pos_hash;
		my $curr_partial_codon_pos = 0;
		my $curr_partial_codon = '';
		
		# For each segment, step through each site, add to length and sequence
		for(my $i = 1; $i <= $num_segments; $i++) {
			my $this_start_key = 'start_' . $i;
			my $this_stop_key = 'stop_' . $i;
			
			my $start_site = $hh_product_position_info{$curr_product}->{$this_start_key};
			my $start_index = ($start_site-1);
			my $stop_site = $hh_product_position_info{$curr_product}->{$this_stop_key};
			my $stop_index = ($stop_site-1);
			my $coding_length = ($stop_site - $start_site + 1);
			$segment_lengths{$i} = $coding_length;
			$coding_length_sum += $coding_length;
			#$segment_point_remainders{$i} = ($coding_length_sum % 3); # this will save CURRENTS for each break point
			$coding_seq .= substr($seq,$start_index,$coding_length);
			
			if($curr_partial_codon_pos == 0) { # Could be FIRST or any %3==0 segment
				$curr_partial_codon_pos = $start_site;
			}
			
			for (my $j = $start_site; $j <= $stop_site; $j++) { # $j goes from start to end POSITION (not INDEX) of segment $i in FASTA-derived sequence
				my $j_index = $j-1;
				my $this_nt = substr($seq,$j_index,1);
				$curr_partial_codon .= $this_nt;
				push(@coding_seq_pos_arr,$j);
				
				if(length($curr_partial_codon) == 3) { # once $curr_partial_codon reaches 3-nt
					# SAVE THE CODON
					$codon_by_pos_hash{$curr_partial_codon_pos} = $curr_partial_codon;
					
					if ( $i == 1 && # FIRST segment
						$j == ($start_site+2) && # end of FIRST codon
						$curr_partial_codon ne 'ATG' && $curr_partial_codon ne 'AUG' && # it ISN'T a START codon
						(! exists $seen_no_start_hash{$curr_product})) { # we haven't yet noticed there's no START
						
						print "\n## WARNING: Please be aware that the product ".
								"$curr_product does not begin with a START (ATG) codon at site $start_site, but rather $curr_partial_codon.".
								"\n## If this was unexpected, please check your annotations.\n";
								
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | CODON | WARNING
						print ERROR_FILE "N/A\t$curr_product\t$start_site\tDoes not begin with a ".
						"START (ATG) codon, but rather $curr_partial_codon. If this was unexpected, please check your annotations\n";
						close ERROR_FILE;
						chdir('..');
								
						$seen_no_start_hash{$curr_product} = 1;
					}
					
					# WARNING if first codon isn't ATG, but proceed #comeback this necessary?
					if (($i == 1) && ($j == $start_site) && # first segment, first site
						($curr_partial_codon ne 'ATG' && $curr_partial_codon ne 'AUG') && # it ISN'T a START codon, and
						(! exists $seen_no_start_hash{$curr_product})) { # we haven't yet noticed there's no START
						
						print "\n## WARNING: Please be aware that the product ".
								"$curr_product does not begin with a START (ATG) codon at site $start_site, but rather $curr_partial_codon.".
								"\n## If this was unexpected, please check your annotations.\n";
								
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | CODON | WARNING
						print ERROR_FILE "N/A\t$curr_product\t$start_site\tDoes not begin with a ".
						"START (ATG) codon, but rather $curr_partial_codon. If this was unexpected, please check your annotations\n";
						close ERROR_FILE;
						chdir('..');
								
						$seen_no_start_hash{$curr_product} = 1;
					}
					
					# WARNING if there is a mid-sequence STOP codon, but proceed
					if ($i < $num_segments && # NOT the last segment
						($curr_partial_codon eq "TAA" || # it IS a stop codon
						$curr_partial_codon eq "TAG" ||
						$curr_partial_codon eq "TGA" ||
						$curr_partial_codon eq "UAA" ||
						$curr_partial_codon eq "UAG" ||
						$curr_partial_codon eq "UGA")
						) { 
					
							$seen_product_early_stop_hash{$curr_product} = 1;
							print "\n## WARNING: Please be aware that there is a ".
								"mid-sequence STOP codon in the sequence ".
								"$curr_product|$j.\n## Please check your annotations for: ".
								"(1) incorrect frame; or (2) incorrect starting or ending coordinates.".
								"\n## A premature STOP codon may also indicate a pseudogene, for which ".
								"piN vs. piS analysis may not be appropriate.\n";
								
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							# FILE | PRODUCT | SITE | CODON | WARNING
							print ERROR_FILE "N/A\t$curr_product\t$j\tMid-sequence ".
							"STOP codon. Please check your annotations for: (1) incorrect ".
							"frame; or (2) incorrect starting or ending coordinates. A ".
							"premature STOP codon may also indicate a pseudogene, for which ".
							"piN vs. piS analysis may not be appropriate.\n";
							close ERROR_FILE;
							chdir('..');
						
					} elsif ( (($j < $stop_site) && # we're in last segment and BEFORE last site
						(! exists $seen_product_early_stop_hash{$curr_product})) && # we haven't seen the product yet &&
						($curr_partial_codon eq "TAA" || # it IS a stop codon
						$curr_partial_codon eq "TAG" ||
						$curr_partial_codon eq "TGA" ||
						$curr_partial_codon eq "UAA" ||
						$curr_partial_codon eq "UAG" ||
						$curr_partial_codon eq "UGA")
						) {
							$seen_product_early_stop_hash{$curr_product} = 1;
							print "\n## WARNING: Please be aware that there is a ".
								"mid-sequence STOP codon in the sequence ".
								"$curr_product|$j.\n## Please check your annotations for: ".
								"(1) incorrect frame; or (2) incorrect starting or ending coordinates.".
								"\n## A premature STOP codon may also indicate a pseudogene, for which ".
								"piN vs. piS analysis may not be appropriate.\n";
								
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							# FILE | PRODUCT | SITE | CODON | WARNING
							print ERROR_FILE "N/A\t$curr_product\t$j\tMid-sequence ".
							"STOP codon. Please check your annotations for: (1) incorrect ".
							"frame; or (2) incorrect starting or ending coordinates. A ".
							"premature STOP codon may also indicate a pseudogene, for which ".
							"piN vs. piS analysis may not be appropriate.\n";
							close ERROR_FILE;
							chdir('..');
					}
					
					# RESET CODON
					$curr_partial_codon = '';
					
					# RESET POSITION
					if($j == $stop_site) {
						$curr_partial_codon_pos = 0;
					} else {
						$curr_partial_codon_pos = ($j+1);
					}
				}
			}
		}
		
		# WARNING and DIE if the total coding length wasn't a multiple of 3
		if(($coding_length_sum % 3) != 0) {
			chdir('SNPGenie_Results');
			open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
			# FILE | PRODUCT | SITE | CODON | WARNING
			print ERROR_FILE "N/A\t$curr_product\tN/A\tTotal length is not a multiple".
			" of 3. Please check your annotations; a complete codon set is required.".
			" SNPGenie terminated.\n";
			close ERROR_FILE;
			chdir('..');
			
			die "\n## WARNING: For product $curr_product, the total length is not a multiple\n".
				"## of 3. Please check your annotations; a complete codon set is required.\n".
				"## SNPGenie terminated.\n";
		}
		
		my @codon_positions_sorted = sort {$a <=> $b} (keys %codon_by_pos_hash);
		
		my $total_N_diffs = 0;
		my $total_S_diffs = 0;
		
		my $total_N_sites = 0;
		my $total_S_sites = 0;
		
		# CHANGE DIRECTORY TO RESULTS FOR WRITING
		chdir('SNPGenie_Results');
		
		my $codon_counter = 0;
		my $contiguous_zero_diffs_counter;
		
		my @coding_positions_sorted = sort {$a <=> $b} @coding_seq_pos_arr;
		
		# New segments approach implemented here
		# Go through all codons and sites within, checking for variants, storing the 
		# number of nonsynonymous and synonymous changes.
		for(my $arr_index = 0; $arr_index < scalar(@coding_positions_sorted); $arr_index+=3) { # $arr_index gives index in array for sorted CODING start genome coordinates			
			my $curr_site = $coding_positions_sorted[$arr_index];
			my $site_pos_1 = $curr_site;
			my $site_pos_2 = $coding_positions_sorted[$arr_index+1];
			my $site_pos_3 = $coding_positions_sorted[$arr_index+2];
			
			my $N_diffs_per_comparison = 0;
			my $S_diffs_per_comparison = 0;
			
			my $N_sites = 0;
			my $S_sites = 0;
			
			# GET THE CODON
			my $curr_codon = $codon_by_pos_hash{$curr_site};
			
			# GET THE CODON'S AMINO ACID
			my $curr_amino_acid = &get_amino_acid($curr_codon);
			
			# FETCH THE SITES' AVERAGE COVERAGES, which have already been calculated
			my $site_1_cov;
			my $site_2_cov;
			my $site_3_cov;
			if (exists $hh_product_position_info{$curr_product}->{$site_pos_1}) { 
				# If there was a variant(s) stored at position 1
				$site_1_cov = $hh_product_position_info{$curr_product}->{$site_pos_1}->{cov};
			} else { # no variant at site 1
				$site_1_cov = 1;
			}
			if (exists $hh_product_position_info{$curr_product}->{$site_pos_2}) { 
				# If there was a variant(s) stored at position 2
				$site_2_cov = $hh_product_position_info{$curr_product}->{$site_pos_2}->{cov};
			} else { # no variant at site 2
				$site_2_cov = 1;
			}
			if (exists $hh_product_position_info{$curr_product}->{$site_pos_3}) { 
				# If there was a variant(s) stored at position 3
				$site_3_cov = $hh_product_position_info{$curr_product}->{$site_pos_3}->{cov};
			} else { # no variant at site 3
				$site_3_cov = 1;
			}
			
			my $num_1_A = 0;
			my $num_1_C = 0;
			my $num_1_G = 0;
			my $num_1_T = 0;
			
			my $prop_1_A = 0;
			my $prop_1_C = 0;
			my $prop_1_G = 0;
			my $prop_1_T = 0;
			
			my $num_2_A = 0;
			my $num_2_C = 0;
			my $num_2_G = 0;
			my $num_2_T = 0;
			
			my $prop_2_A = 0;
			my $prop_2_C = 0;
			my $prop_2_G = 0;
			my $prop_2_T = 0;
			
			my $num_3_A = 0;
			my $num_3_C = 0;
			my $num_3_G = 0;
			my $num_3_T = 0;
			
			my $prop_3_A = 0;
			my $prop_3_C = 0;
			my $prop_3_G = 0;
			my $prop_3_T = 0;
			
			my %site_1_num_hash; # contains {A,C,G,T}=>NUM # need info from above
			my %site_2_num_hash; # contains {A,C,G,T}=>NUM # need info from above
			my %site_3_num_hash; # contains {A,C,G,T}=>NUM # need info from above
			
			# IF NO VARIANTS in codon, contribute to the POSSIBILITY of a WARNING
			if ((! exists $hh_product_position_info{$curr_product}->{$site_pos_1}) &&
				(! exists $hh_product_position_info{$curr_product}->{$site_pos_2}) &&
				(! exists $hh_product_position_info{$curr_product}->{$site_pos_3})) {
				$contiguous_zero_diffs_counter++;
			} else {
				$contiguous_zero_diffs_counter = 0;
			}
			
			if ($contiguous_zero_diffs_counter == 1000) {
				print "\n## WARNING: In file $file_nm, in ".
							"$curr_product, near site $curr_site,\n## ".
							"there is a stretch of ≥1000 codons with zero ".
							"SNPs called. If this was unexpected,\n## ".
							"you may need to specify Unix (\\n) newline characters! See ".
							"Troubleshooting.\n";
							
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | WARNING
				print ERROR_FILE "$file_nm\t$curr_product\t$curr_site\t".
					"No SNPs in >=1000 contiguous codons. If this was unexpected, ".
					"you may need to specify Unix (\\n) newline characters! See ".
					"Troubleshooting\n";
				close ERROR_FILE;
			}
			
			####
			## For purposes of calculating AVERAGE dN and dS versus the REFERENCE
			my $this_codon_num_N_diffs = 0;
			my $this_codon_num_S_diffs = 0;
			my $ref_amino_acid = &get_amino_acid($curr_codon);
			my $codon_cov_sum;
			my $codon_cov_denom_sum;
			##
			####
			
			####
			## GENE DIVERSITY stuff
			my $this_codon_site1_gene_diversity;
			my $this_codon_site2_gene_diversity;
			my $this_codon_site3_gene_diversity;
			
			my $this_codon_site1_polymorphic;
			my $this_codon_site2_polymorphic;
			my $this_codon_site3_polymorphic;
			
			my $site1_siteb_category;
			my $site2_siteb_category;
			my $site3_siteb_category;
			
			my $site1_codonb_category;
			my $site2_codonb_category;
			my $site3_codonb_category;
			## END GENE DIVERSITY stuff
			####
			
			########################
			### CODON POSITION 1 ###
			########################
			if (exists $hh_product_position_info{$curr_product}->{$site_pos_1}) { # If there was a variant(s) stored at position 1			
				# call &get_amino_acid(CODON)
				# call &get_number_of_site(CODON,POSITION) to return @arr=(#N,#S)
				
				$num_1_A = $hh_product_position_info{$curr_product}->{$site_pos_1}->{A};
				$num_1_C = $hh_product_position_info{$curr_product}->{$site_pos_1}->{C};
				$num_1_G = $hh_product_position_info{$curr_product}->{$site_pos_1}->{G};
				$num_1_T = $hh_product_position_info{$curr_product}->{$site_pos_1}->{T};
				
				if ($num_1_A != 0) {
					$site_1_num_hash{'A'} = $num_1_A;
				}
				
				if ($num_1_C != 0) {
					$site_1_num_hash{'C'} = $num_1_C;
				}
				
				if ($num_1_G != 0) {
					$site_1_num_hash{'G'} = $num_1_G;
				}
				
				if ($num_1_T != 0) {
					$site_1_num_hash{'T'} = $num_1_T;
				}
				
				#print "\nnum_1_A is $num_1_A\n";
				
				my $site_cov = $hh_product_position_info{$curr_product}->{$site_pos_1}->{cov};
				$site_1_cov = $site_cov;
				
				my $num_pairwise_comparisons = ((($site_cov * $site_cov) - $site_cov)/2);
				
				# NEW PROPORTION APPROACH
				$prop_1_A = $hh_product_position_info{$curr_product}->{$site_pos_1}->{A_prop};
				$prop_1_C = $hh_product_position_info{$curr_product}->{$site_pos_1}->{C_prop};
				$prop_1_G = $hh_product_position_info{$curr_product}->{$site_pos_1}->{G_prop};
				$prop_1_T = $hh_product_position_info{$curr_product}->{$site_pos_1}->{T_prop};
				
#				my $A_prop = ($num_1_A / $site_cov);
#				my $C_prop = ($num_1_C / $site_cov);
#				my $G_prop = ($num_1_G / $site_cov);
#				my $T_prop = ($num_1_T / $site_cov);
				
				# NEW PROPORTION APPROACH
				my $A_prop = $prop_1_A;
				my $C_prop = $prop_1_C;
				my $G_prop = $prop_1_G;
				my $T_prop = $prop_1_T;
				$num_1_A = ($A_prop * $site_1_cov);
				$num_1_C = ($C_prop * $site_1_cov);
				$num_1_G = ($G_prop * $site_1_cov);
				$num_1_T = ($T_prop * $site_1_cov);
				
				#print "\nA_prop is $A_prop\n";
				
				# Warn if the total proportion doesn't equal 1 at a site
				my $total_prop = $A_prop + $C_prop + $G_prop + $T_prop;
				if (($total_prop < 0.999999) || ($total_prop > 1.000001)) {
					my $total_percentage = $total_prop*100;
					my $tot_perc_rounded = sprintf("%.2f",$total_percentage);
					
					print "\n## WARNING: In $file_nm|$curr_product|".
							"$curr_site,\n## the nucleotide total (which should ".
							"be 100.00%) is instead: $tot_perc_rounded\%.\n## This should occur ".
							"only when conflicting coverages have been reported.\n";
					
					open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
					# FILE | PRODUCT | SITE | WARNING
					print ERROR_FILE "$file_nm\t$curr_product\t$curr_site\t".
						"Nucleotide total does not equal 100.00% of coverage, but instead ".
						"$tot_perc_rounded\%. This should occur only when conflicting ".
						"coverages have been reported\n";
					close ERROR_FILE;
				}
				
				# Construct the codon with different nucleotides in position 1
				my $codon_1_A = "A".substr($curr_codon,1,2);
				my $codon_1_C = "C".substr($curr_codon,1,2);
				my $codon_1_G = "G".substr($curr_codon,1,2);
				my $codon_1_T = "T".substr($curr_codon,1,2);
				
				my $amino_acid_1_A = &get_amino_acid($codon_1_A);
				my $amino_acid_1_C = &get_amino_acid($codon_1_C);
				my $amino_acid_1_G = &get_amino_acid($codon_1_G);
				my $amino_acid_1_T = &get_amino_acid($codon_1_T);
				
				#print "\n\nWith A: $codon_1_A\nWith C: $codon_1_C\nWith G: $codon_1_G\n".
				#	"With T: $codon_1_T\n\n";
				
				# Get the number of sites with A in pos 1
				my @num_sites_1_A_arr = &get_number_of_sites($codon_1_A,1); #to return @arr=(#N,#S)
				my @num_sites_1_C_arr = &get_number_of_sites($codon_1_C,1); #to return @arr=(#N,#S)
				my @num_sites_1_G_arr = &get_number_of_sites($codon_1_G,1); #to return @arr=(#N,#S)
				my @num_sites_1_T_arr = &get_number_of_sites($codon_1_T,1); #to return @arr=(#N,#S)
				
				#print "\nnum_sites_1_A_arr is: @num_sites_1_A_arr\n";
				
				$N_sites += (($num_sites_1_A_arr[0] * $A_prop) + ($num_sites_1_C_arr[0] * $C_prop) +
					($num_sites_1_G_arr[0] * $G_prop) + ($num_sites_1_T_arr[0] * $T_prop));
					
				#print "\nN_sites is $N_sites\n";
					
				$S_sites += (($num_sites_1_A_arr[1] * $A_prop) + ($num_sites_1_C_arr[1] * $C_prop) +
					($num_sites_1_G_arr[1] * $G_prop) + ($num_sites_1_T_arr[1] * $T_prop));
				
				my $N_diffs = 0;
				my $S_diffs = 0;
				
				my $AC_type;
				my $AG_type;
				my $AT_type;
				my $CG_type;
				my $CT_type;
				my $GT_type;
				
				my $AC = ($num_1_A * $num_1_C);
				my $AG = ($num_1_A * $num_1_G);
				my $AT = ($num_1_A * $num_1_T);
				my $CG = ($num_1_C * $num_1_G);
				my $CT = ($num_1_C * $num_1_T);
				my $GT = ($num_1_G * $num_1_T);
				
				if ($amino_acid_1_A eq $amino_acid_1_C) {
					$AC_type = 'S';
					$S_diffs += $AC;
				} else {
					$AC_type = 'N';
					$N_diffs += $AC;
				}
				
				if ($amino_acid_1_A eq $amino_acid_1_G) {
					$AG_type = 'S';
					$S_diffs += $AG;
				} else {
					$AG_type = 'N';
					$N_diffs += $AG;
				}
				
				if ($amino_acid_1_A eq $amino_acid_1_T) {
					$AT_type = 'S';
					$S_diffs += $AT;
				} else {
					$AT_type = 'N';
					$N_diffs += $AT;
				}
				
				if ($amino_acid_1_C eq $amino_acid_1_G) {
					$CG_type = 'S';
					$S_diffs += $CG;
				} else {
					$CG_type = 'N';
					$N_diffs += $CG;
				}
				
				if ($amino_acid_1_C eq $amino_acid_1_T) {
					$CT_type = 'S';
					$S_diffs += $CT;
				} else {
					$CT_type = 'N';
					$N_diffs += $CT;
				}
				
				if ($amino_acid_1_G eq $amino_acid_1_T) {
					$GT_type = 'S';
					$S_diffs += $GT;
				} else {
					$GT_type = 'N';
					$N_diffs += $GT;
				}
				
				$N_diffs_per_comparison += ($N_diffs / $num_pairwise_comparisons);
				$S_diffs_per_comparison += ($S_diffs / $num_pairwise_comparisons);
				
				####
				## AVERAGE dN and dS versus REFERENCE stuff
				# Need to make sure the reference isn't a STOP codon
				my $ref_nt = $hh_product_position_info{$curr_product}->{$site_pos_1}->{reference};
				if (($num_1_A != 0) && ($ref_nt ne 'A')) {
					if($amino_acid_1_A eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_1_A;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_1_A ne 'TAA') && ($codon_1_A ne 'TAG') && ($codon_1_A ne 'TGA') &&
						($codon_1_A ne 'UAA') && ($codon_1_A ne 'UAG') && ($codon_1_A ne 'UGA')) {
						$this_codon_num_N_diffs += $num_1_A;
					}
				}
				
				if (($num_1_C != 0) && ($ref_nt ne 'C')) {
					if($amino_acid_1_C eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_1_C;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_1_C ne 'TAA') && ($codon_1_C ne 'TAG') && ($codon_1_C ne 'TGA') &&
						($codon_1_C ne 'UAA') && ($codon_1_C ne 'UAG') && ($codon_1_C ne 'UGA')) {
						$this_codon_num_N_diffs += $num_1_C;
					}
				}
				
				if (($num_1_G != 0) && ($ref_nt ne 'G')) {
					if($amino_acid_1_G eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_1_G;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_1_G ne 'TAA') && ($codon_1_G ne 'TAG') && ($codon_1_G ne 'TGA') &&
						($codon_1_G ne 'UAA') && ($codon_1_G ne 'UAG') && ($codon_1_G ne 'UGA')) {
						$this_codon_num_N_diffs += $num_1_G;
					}
				}
				
				if (($num_1_T != 0) && ($ref_nt ne 'T')) {
					if($amino_acid_1_T eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_1_T;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_1_T ne 'TAA') && ($codon_1_T ne 'TAG') && ($codon_1_T ne 'TGA') &&
						($codon_1_T ne 'UAA') && ($codon_1_T ne 'UAG') && ($codon_1_T ne 'UGA')) {
						$this_codon_num_N_diffs += $num_1_T;
					}
				}
				$codon_cov_sum += $site_1_cov;
				$codon_cov_denom_sum += 1;
				## END AVERAGE dN and dS versus REFERENCE stuff
				####
				
				####
				## GENE DIVERSITY stuff
				#my $ref_nt = $hh_product_position_info{$curr_product}->{$site_pos_1}->{reference}; # already done above
				my $this_site_category;
				my $this_site_N;
				my $this_site_S;
				my $this_site_sq_var_freq_sum;
				if ($num_1_A != 0) {
					if ($ref_nt ne 'A') {
						my $A_freq = ($num_1_A / $site_1_cov);
						my $A_freq_sq = ($A_freq ** 2);
						$this_site_sq_var_freq_sum += $A_freq_sq;
						if($amino_acid_1_A eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $A_freq = ($num_1_A / $site_1_cov);
						my $A_freq_sq = ($A_freq ** 2);
						$this_site_sq_var_freq_sum += $A_freq_sq;
					}
				}
				
				if ($num_1_C != 0) {
					if ($ref_nt ne 'C') {
						my $C_freq = ($num_1_C / $site_1_cov);
						my $C_freq_sq = ($C_freq ** 2);
						$this_site_sq_var_freq_sum += $C_freq_sq;
						if($amino_acid_1_C eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $C_freq = ($num_1_C / $site_1_cov);
						my $C_freq_sq = ($C_freq ** 2);
						$this_site_sq_var_freq_sum += $C_freq_sq;
					}
				}
				
				if ($num_1_G != 0) {
					if ($ref_nt ne 'G') {
						my $G_freq = ($num_1_G / $site_1_cov);
						my $G_freq_sq = ($G_freq ** 2);
						$this_site_sq_var_freq_sum += $G_freq_sq;
						if($amino_acid_1_G eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $G_freq = ($num_1_G / $site_1_cov);
						my $G_freq_sq = ($G_freq ** 2);
						$this_site_sq_var_freq_sum += $G_freq_sq;
					}
				}
				
				if ($num_1_T != 0) {
					if ($ref_nt ne 'T') {
						my $T_freq = ($num_1_T / $site_1_cov);
						my $T_freq_sq = ($T_freq ** 2);
						$this_site_sq_var_freq_sum += $T_freq_sq;
						if($amino_acid_1_T eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $T_freq = ($num_1_T / $site_1_cov);
						my $T_freq_sq = ($T_freq ** 2);
						$this_site_sq_var_freq_sum += $T_freq_sq;
					}
				}
				
				my $this_site_category_sum = ($this_site_N + $this_site_S);
				if($this_site_category_sum == 0) {
					$site1_siteb_category = 0;
				} elsif ($this_site_category_sum == 2) {
					$site1_siteb_category = 3;
				} else {
					if($this_site_S == 1) {
						$site1_siteb_category = 1;
					} else {
						$site1_siteb_category = 2;
					}
				}
				
				my $this_site_gene_diversity = (1 - $this_site_sq_var_freq_sum);
				$this_codon_site1_gene_diversity = $this_site_gene_diversity;
				$this_codon_site1_polymorphic = 1;
				## END GENE DIVERSITY stuff
				####
					
			} else {
				my @num_sites_arr = &get_number_of_sites($curr_codon,1);
				$N_sites += $num_sites_arr[0];
				$S_sites += $num_sites_arr[1];
				
				my $N_diffs = 0;
				my $S_diffs = 0;
				
				$N_diffs_per_comparison += 0;
				$S_diffs_per_comparison += 0;
				
				# Find out what nucleotide it is and give it a num of 1
				if ($seq_by_index_arr[$site_pos_1-1] eq 'A') {
					$num_1_A = 1;
					$site_1_num_hash{'A'} = $num_1_A;
				} elsif ($seq_by_index_arr[$site_pos_1-1] eq 'C') {
					$num_1_C = 1;
					$site_1_num_hash{'C'} = $num_1_C;
				} elsif ($seq_by_index_arr[$site_pos_1-1] eq 'G') {
					$num_1_G = 1;
					$site_1_num_hash{'G'} = $num_1_G;
				} elsif ($seq_by_index_arr[$site_pos_1-1] eq 'T') {
					$num_1_T = 1;
					$site_1_num_hash{'T'} = $num_1_T;
				}
				
				# Set this site's cov to 1, allowing site/cov to be 100%
				$site_1_cov = 1;
				
				####
				## GENE DIVERSITY stuff
				$this_codon_site1_gene_diversity = 0;
				$this_codon_site1_polymorphic = 0;
				$site1_siteb_category = 0;
				$site1_codonb_category = 0;
				## END GENE DIVERSITY stuff
				####
			}
			
			#print "In: $curr_product, site $curr_site, the nt is ".$seq_by_pos_hash{$site_pos_1}.
			#	", and A=". $site_1_num_hash{'A'} . " / C=".$site_1_num_hash{'C'} .
			#	" / G=".$site_1_num_hash{'G'} ." / T=".$site_1_num_hash{'T'}."\n";
			
			########################
			### CODON POSITION 2 ###
			########################
			if (exists $hh_product_position_info{$curr_product}->{$site_pos_2}) { # If there was a variant(s) stored at position 2
				$num_2_A = $hh_product_position_info{$curr_product}->{$site_pos_2}->{A};
				$num_2_C = $hh_product_position_info{$curr_product}->{$site_pos_2}->{C};
				$num_2_G = $hh_product_position_info{$curr_product}->{$site_pos_2}->{G};
				$num_2_T = $hh_product_position_info{$curr_product}->{$site_pos_2}->{T};
				
				if ($num_2_A != 0) {
					$site_2_num_hash{'A'} = $num_2_A;
				}
				
				if ($num_2_C != 0) {
					$site_2_num_hash{'C'} = $num_2_C;
				}
				
				if ($num_2_G != 0) {
					$site_2_num_hash{'G'} = $num_2_G;
				}
				
				if ($num_2_T != 0) {
					$site_2_num_hash{'T'} = $num_2_T;
				}
				
				#print "\nnum_2_A is $num_2_A\n";
				
				my $site_cov = $hh_product_position_info{$curr_product}->{$site_pos_2}->{cov};
				$site_2_cov = $site_cov;
				
				my $num_pairwise_comparisons = ((($site_cov * $site_cov) - $site_cov)/2);
				
				# NEW PROPORTION APPROACH
				$prop_2_A = $hh_product_position_info{$curr_product}->{$site_pos_2}->{A_prop};
				$prop_2_C = $hh_product_position_info{$curr_product}->{$site_pos_2}->{C_prop};
				$prop_2_G = $hh_product_position_info{$curr_product}->{$site_pos_2}->{G_prop};
				$prop_2_T = $hh_product_position_info{$curr_product}->{$site_pos_2}->{T_prop};
				
#				my $A_prop = ($num_2_A / $site_cov);
#				my $C_prop = ($num_2_C / $site_cov);
#				my $G_prop = ($num_2_G / $site_cov);
#				my $T_prop = ($num_2_T / $site_cov);
				
				# NEW PROPORTION APPROACH
				my $A_prop = $prop_2_A;
				my $C_prop = $prop_2_C;
				my $G_prop = $prop_2_G;
				my $T_prop = $prop_2_T;
				$num_2_A = ($A_prop * $site_2_cov);
				$num_2_C = ($C_prop * $site_2_cov);
				$num_2_G = ($G_prop * $site_2_cov);
				$num_2_T = ($T_prop * $site_2_cov);
				
				#print "\nA_prop is $A_prop\n";
				
				# Warn if the total proportion doesn't equal 1 at a site
				my $total_prop = $A_prop + $C_prop + $G_prop + $T_prop;
				if (($total_prop < 0.999999) || ($total_prop > 1.000001)) {
					my $total_percentage = $total_prop*100;
					my $tot_perc_rounded = sprintf("%.2f",$total_percentage);
					
					print "\n## WARNING: In $file_nm|$curr_product|".
							"$curr_site,\n## the nucleotide total (which should ".
							"be 100.00%) is instead: $tot_perc_rounded\%.\n## This should occur ".
							"only when conflicting coverages have been reported.\n";
					
					#chdir('SNPGenie_Results');
					open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
					# FILE | PRODUCT | SITE | CODON | WARNING
					print ERROR_FILE "$file_nm\t$curr_product\t$curr_site\t".
						"Nucleotide total does not equal 100.00% of coverage, but instead ".
						"$tot_perc_rounded\%. This should occur only when conflicting ".
						"coverages have been reported\n";
					close ERROR_FILE;
					#chdir('..');
				}
				
				# Construct the codon with different nucleotides in position 2
				my $codon_2_A = "".substr($curr_codon,0,1)."A".substr($curr_codon,2,1);
				my $codon_2_C = "".substr($curr_codon,0,1)."C".substr($curr_codon,2,1);
				my $codon_2_G = "".substr($curr_codon,0,1)."G".substr($curr_codon,2,1);
				my $codon_2_T = "".substr($curr_codon,0,1)."T".substr($curr_codon,2,1);
				
				my $amino_acid_2_A = &get_amino_acid($codon_2_A);
				my $amino_acid_2_C = &get_amino_acid($codon_2_C);
				my $amino_acid_2_G = &get_amino_acid($codon_2_G);
				my $amino_acid_2_T = &get_amino_acid($codon_2_T);
				
				#print "\n\nWith A: $codon_2_A\nWith C: $codon_2_C\nWith G: $codon_2_G\n".
				#	"With T: $codon_2_T\n\n";
				
				# Get the number of sites with A in pos 2
				my @num_sites_2_A_arr = &get_number_of_sites($codon_2_A,2); #to return @arr=(#N,#S)
				my @num_sites_2_C_arr = &get_number_of_sites($codon_2_C,2); #to return @arr=(#N,#S)
				my @num_sites_2_G_arr = &get_number_of_sites($codon_2_G,2); #to return @arr=(#N,#S)
				my @num_sites_2_T_arr = &get_number_of_sites($codon_2_T,2); #to return @arr=(#N,#S)
				
				#print "\nnum_sites_2_A_arr is: @num_sites_2_A_arr\n";
				
				$N_sites += (($num_sites_2_A_arr[0] * $A_prop) + ($num_sites_2_C_arr[0] * $C_prop) +
					($num_sites_2_G_arr[0] * $G_prop) + ($num_sites_2_T_arr[0] * $T_prop));
					
				#print "\nN_sites is $N_sites\n";
					
				$S_sites += (($num_sites_2_A_arr[1] * $A_prop) + ($num_sites_2_C_arr[1] * $C_prop) +
					($num_sites_2_G_arr[1] * $G_prop) + ($num_sites_2_T_arr[1] * $T_prop));
				
				my $N_diffs = 0;
				my $S_diffs = 0;
				
				my $AC_type;
				my $AG_type;
				my $AT_type;
				my $CG_type;
				my $CT_type;
				my $GT_type;
				
				my $AC = ($num_2_A * $num_2_C);
				my $AG = ($num_2_A * $num_2_G);
				my $AT = ($num_2_A * $num_2_T);
				my $CG = ($num_2_C * $num_2_G);
				my $CT = ($num_2_C * $num_2_T);
				my $GT = ($num_2_G * $num_2_T);
				
				if ($amino_acid_2_A eq $amino_acid_2_C) {
					$AC_type = 'S';
					$S_diffs += $AC;
				} else {
					$AC_type = 'N';
					$N_diffs += $AC;
				}
				
				if ($amino_acid_2_A eq $amino_acid_2_G) {
					$AG_type = 'S';
					$S_diffs += $AG;
				} else {
					$AG_type = 'N';
					$N_diffs += $AG;
				}
				
				if ($amino_acid_2_A eq $amino_acid_2_T) {
					$AT_type = 'S';
					$S_diffs += $AT;
				} else {
					$AT_type = 'N';
					$N_diffs += $AT;
				}
				
				if ($amino_acid_2_C eq $amino_acid_2_G) {
					$CG_type = 'S';
					$S_diffs += $CG;
				} else {
					$CG_type = 'N';
					$N_diffs += $CG;
				}
				
				if ($amino_acid_2_C eq $amino_acid_2_T) {
					$CT_type = 'S';
					$S_diffs += $CT;
				} else {
					$CT_type = 'N';
					$N_diffs += $CT;
				}
				
				if ($amino_acid_2_G eq $amino_acid_2_T) {
					$GT_type = 'S';
					$S_diffs += $GT;
				} else {
					$GT_type = 'N';
					$N_diffs += $GT;
				}
				
				$N_diffs_per_comparison += ($N_diffs / $num_pairwise_comparisons);
				$S_diffs_per_comparison += ($S_diffs / $num_pairwise_comparisons);
				
				## AVERAGE dN and dS versus REFERENCE stuff
				my $ref_nt = $hh_product_position_info{$curr_product}->{$site_pos_2}->{reference};
				if (($num_2_A != 0) && ($ref_nt ne 'A')) {
					if($amino_acid_2_A eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_2_A;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_2_A ne 'TAA') && ($codon_2_A ne 'TAG') && ($codon_2_A ne 'TGA') &&
						($codon_2_A ne 'UAA') && ($codon_2_A ne 'UAG') && ($codon_2_A ne 'UGA')) {
						$this_codon_num_N_diffs += $num_2_A;
					}
				}
				
				if (($num_2_C != 0) && ($ref_nt ne 'C')) {
					if($amino_acid_2_C eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_2_C;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_2_C ne 'TAA') && ($codon_2_C ne 'TAG') && ($codon_2_C ne 'TGA') &&
						($codon_2_C ne 'UAA') && ($codon_2_C ne 'UAG') && ($codon_2_C ne 'UGA')) {
						$this_codon_num_N_diffs += $num_2_C;
					}
				}
				
				if (($num_2_G != 0) && ($ref_nt ne 'G')) {
					if($amino_acid_2_G eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_2_G;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_2_G ne 'TAA') && ($codon_2_G ne 'TAG') && ($codon_2_G ne 'TGA') &&
						($codon_2_G ne 'UAA') && ($codon_2_G ne 'UAG') && ($codon_2_G ne 'UGA')) {
						$this_codon_num_N_diffs += $num_2_G;
					}
				}
				
				if (($num_2_T != 0) && ($ref_nt ne 'T')) {
					if($amino_acid_2_T eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_2_T;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_2_T ne 'TAA') && ($codon_2_T ne 'TAG') && ($codon_2_T ne 'TGA') &&
						($codon_2_T ne 'UAA') && ($codon_2_T ne 'UAG') && ($codon_2_T ne 'UGA')) {
						$this_codon_num_N_diffs += $num_2_T;
					}
				}
				$codon_cov_sum += $site_2_cov;
				$codon_cov_denom_sum += 1;
				## END AVERAGE dN and dS versus REFERENCE stuff
				
				####
				## GENE DIVERSITY stuff
				#my $ref_nt = $hh_product_position_info{$curr_product}->{$site_pos_2}->{reference}; # already done above
				my $this_site_category;
				my $this_site_N;
				my $this_site_S;
				my $this_site_sq_var_freq_sum;
				if ($num_2_A != 0) {
					if ($ref_nt ne 'A') {
						my $A_freq = ($num_2_A / $site_2_cov);
						my $A_freq_sq = ($A_freq ** 2);
						$this_site_sq_var_freq_sum += $A_freq_sq;
						if($amino_acid_2_A eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $A_freq = ($num_2_A / $site_2_cov);
						my $A_freq_sq = ($A_freq ** 2);
						$this_site_sq_var_freq_sum += $A_freq_sq;
					}
				}
				
				if ($num_2_C != 0) {
					if ($ref_nt ne 'C') {
						my $C_freq = ($num_2_C / $site_2_cov);
						my $C_freq_sq = ($C_freq ** 2);
						$this_site_sq_var_freq_sum += $C_freq_sq;
						if($amino_acid_2_C eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $C_freq = ($num_2_C / $site_2_cov);
						my $C_freq_sq = ($C_freq ** 2);
						$this_site_sq_var_freq_sum += $C_freq_sq;
					}
				}
				
				if ($num_2_G != 0) {
					if ($ref_nt ne 'G') {
						my $G_freq = ($num_2_G / $site_2_cov);
						my $G_freq_sq = ($G_freq ** 2);
						$this_site_sq_var_freq_sum += $G_freq_sq;
						if($amino_acid_2_G eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $G_freq = ($num_2_G / $site_2_cov);
						my $G_freq_sq = ($G_freq ** 2);
						$this_site_sq_var_freq_sum += $G_freq_sq;
					}
				}
				
				if ($num_2_T != 0) {
					if ($ref_nt ne 'T') {
						my $T_freq = ($num_2_T / $site_2_cov);
						my $T_freq_sq = ($T_freq ** 2);
						$this_site_sq_var_freq_sum += $T_freq_sq;
						if($amino_acid_2_T eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $T_freq = ($num_2_T / $site_2_cov);
						my $T_freq_sq = ($T_freq ** 2);
						$this_site_sq_var_freq_sum += $T_freq_sq;
					}
				}
				
				my $this_site_category_sum = ($this_site_N + $this_site_S);
				if($this_site_category_sum == 0) {
					$site2_siteb_category = 0;
				} elsif ($this_site_category_sum == 2) {
					$site2_siteb_category = 3;
				} else {
					if($this_site_S == 1) {
						$site2_siteb_category = 1;
					} else {
						$site2_siteb_category = 2;
					}
				}
				
				my $this_site_gene_diversity = (1 - $this_site_sq_var_freq_sum);
				$this_codon_site2_gene_diversity = $this_site_gene_diversity;
				$this_codon_site2_polymorphic = 1;
				## END GENE DIVERSITY stuff
				####
							
			} else {
				my @num_sites_arr = &get_number_of_sites($curr_codon,2);
				$N_sites += $num_sites_arr[0];
				$S_sites += $num_sites_arr[1];
				
				$N_diffs_per_comparison += 0;
				$S_diffs_per_comparison += 0;	
				
				# Find out what nucleotide it is and give it a num of 1
				if ($seq_by_index_arr[$site_pos_2-1] eq 'A') {
					$num_2_A = 1;
					$site_2_num_hash{'A'} = $num_2_A;
				} elsif ($seq_by_index_arr[$site_pos_2-1] eq 'C') {
					$num_2_C = 1;
					$site_2_num_hash{'C'} = $num_2_C;
				} elsif ($seq_by_index_arr[$site_pos_2-1] eq 'G') {
					$num_2_G = 1;
					$site_2_num_hash{'G'} = $num_2_G;
				} elsif ($seq_by_index_arr[$site_pos_2-1] eq 'T') {
					$num_2_T = 1;
					$site_2_num_hash{'T'} = $num_2_T;
				}
				
				# Set this site's cov to 1, allowing site/cov to be 100%
				$site_2_cov = 1;
				
				####
				## GENE DIVERSITY stuff
				$this_codon_site2_gene_diversity = 0;
				$this_codon_site2_polymorphic = 0;
				$site2_siteb_category = 0;
				$site2_codonb_category = 0;
				## END GENE DIVERSITY stuff
				####	
			}
			
			#print "in $curr_product, site $site_pos_2, the nt is ".$seq_by_pos_hash{$site_pos_2}.", and  site_2_num_hash A is ". $site_2_num_hash{'A'} . "\n";
			
			#print "In: $curr_product, site $site_pos_2, the nt is ".$seq_by_pos_hash{$site_pos_2}.
			#	", and A=". $site_2_num_hash{'A'} . " / C=".$site_2_num_hash{'C'} .
			#	" / G=".$site_2_num_hash{'G'} ." / T=".$site_2_num_hash{'T'}."\n";
			
			########################
			### CODON POSITION 3 ###
			########################
			if (exists $hh_product_position_info{$curr_product}->{$site_pos_3}) { # If there was a variant(s) stored at position 1			
				# call &get_amino_acid(CODON)
				# call &get_number_of_site(CODON,POSITION) to return @arr=(#N,#S)
				
				$num_3_A = $hh_product_position_info{$curr_product}->{$site_pos_3}->{A};
				$num_3_C = $hh_product_position_info{$curr_product}->{$site_pos_3}->{C};
				$num_3_G = $hh_product_position_info{$curr_product}->{$site_pos_3}->{G};
				$num_3_T = $hh_product_position_info{$curr_product}->{$site_pos_3}->{T};
				
				#print "\nnum_1_A is $num_1_A\n";
				
				if ($num_3_A != 0) {
					$site_3_num_hash{'A'} = $num_3_A;
				}
				
				if ($num_3_C != 0) {
					$site_3_num_hash{'C'} = $num_3_C;
				}
				
				if ($num_3_G != 0) {
					$site_3_num_hash{'G'} = $num_3_G;
				}
				
				if ($num_3_T != 0) {
					$site_3_num_hash{'T'} = $num_3_T;
				}
				
				my $site_cov = $hh_product_position_info{$curr_product}->{$site_pos_3}->{cov};
				$site_3_cov = $site_cov;
				
				my $num_pairwise_comparisons = ((($site_cov * $site_cov) - $site_cov)/2);
				
				# NEW PROPORTION APPROACH
				$prop_3_A = $hh_product_position_info{$curr_product}->{$site_pos_3}->{A_prop};
				$prop_3_C = $hh_product_position_info{$curr_product}->{$site_pos_3}->{C_prop};
				$prop_3_G = $hh_product_position_info{$curr_product}->{$site_pos_3}->{G_prop};
				$prop_3_T = $hh_product_position_info{$curr_product}->{$site_pos_3}->{T_prop};		
				
#				my $A_prop = ($num_3_A / $site_cov);
#				my $C_prop = ($num_3_C / $site_cov);
#				my $G_prop = ($num_3_G / $site_cov);
#				my $T_prop = ($num_3_T / $site_cov);

				# NEW PROPORTION APPROACH
				my $A_prop = $prop_3_A;
				my $C_prop = $prop_3_C;
				my $G_prop = $prop_3_G;
				my $T_prop = $prop_3_T;
				$num_3_A = ($A_prop * $site_3_cov);
				$num_3_C = ($C_prop * $site_3_cov);
				$num_3_G = ($G_prop * $site_3_cov);
				$num_3_T = ($T_prop * $site_3_cov);
				
				#print "\nA_prop is $A_prop\n";
				
				# Warn if the total proportion doesn't equal 1 at a site
				my $total_prop = $A_prop + $C_prop + $G_prop + $T_prop;
				if (($total_prop < 0.999999) || ($total_prop > 1.000001)) {
					my $total_percentage = $total_prop*100;
					my $tot_perc_rounded = sprintf("%.2f",$total_percentage);
					
					print "\n## WARNING: In $file_nm|$curr_product|".
							"$curr_site,\n## the nucleotide total (which should ".
							"be 100.00%) is instead: $tot_perc_rounded\%.\n## This should occur ".
							"only when conflicting coverages have been reported.\n";
					
					#chdir('SNPGenie_Results');
					open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
					# FILE | PRODUCT | SITE | CODON | WARNING
					print ERROR_FILE "$file_nm\t$curr_product\t$curr_site\t".
						"Nucleotide total does not equal 100.00% of coverage, but instead ".
						"$tot_perc_rounded\%. This should occur only when conflicting ".
						"coverages have been reported\n";
					close ERROR_FILE;
					#chdir('..');
				}
				
				# Construct the codon with different nucleotides in position 1
				my $codon_3_A = "".substr($curr_codon,0,2)."A";
				my $codon_3_C = "".substr($curr_codon,0,2)."C";
				my $codon_3_G = "".substr($curr_codon,0,2)."G";
				my $codon_3_T = "".substr($curr_codon,0,2)."T";
				
				my $amino_acid_3_A = &get_amino_acid($codon_3_A);
				my $amino_acid_3_C = &get_amino_acid($codon_3_C);
				my $amino_acid_3_G = &get_amino_acid($codon_3_G);
				my $amino_acid_3_T = &get_amino_acid($codon_3_T);
				
				#print "\n\nWith A: $codon_3_A\nWith C: $codon_3_C\nWith G: $codon_3_G\n".
				#	"With T: $codon_3_T\n\n";
				
				# Get the number of sites with A in pos 1
				my @num_sites_3_A_arr = &get_number_of_sites($codon_3_A,3); #to return @arr=(#N,#S)
				my @num_sites_3_C_arr = &get_number_of_sites($codon_3_C,3); #to return @arr=(#N,#S)
				my @num_sites_3_G_arr = &get_number_of_sites($codon_3_G,3); #to return @arr=(#N,#S)
				my @num_sites_3_T_arr = &get_number_of_sites($codon_3_T,3); #to return @arr=(#N,#S)
				
				#print "\nnum_sites_3_A_arr is: @num_sites_3_A_arr\n";
				
				$N_sites += (($num_sites_3_A_arr[0] * $A_prop) + ($num_sites_3_C_arr[0] * $C_prop) +
					($num_sites_3_G_arr[0] * $G_prop) + ($num_sites_3_T_arr[0] * $T_prop));
					
				#print "\nN_sites is $N_sites\n";
					
				$S_sites += (($num_sites_3_A_arr[1] * $A_prop) + ($num_sites_3_C_arr[1] * $C_prop) +
					($num_sites_3_G_arr[1] * $G_prop) + ($num_sites_3_T_arr[1] * $T_prop));
				
				my $N_diffs = 0;
				my $S_diffs = 0;
				
				#my $AC_type;
				#my $AG_type;
				#my $AT_type;
				#my $CG_type;
				#my $CT_type;
				#my $GT_type;
				
				my $AC = ($num_3_A * $num_3_C);
				my $AG = ($num_3_A * $num_3_G);
				my $AT = ($num_3_A * $num_3_T);
				my $CG = ($num_3_C * $num_3_G);
				my $CT = ($num_3_C * $num_3_T);
				my $GT = ($num_3_G * $num_3_T);
				
				# SITE-BASED
				if ($amino_acid_3_A eq $amino_acid_3_C) {
					#$AC_type = 'S';
					$S_diffs += $AC;
				} else {
					#$AC_type = 'N';
					$N_diffs += $AC;
				}
				
				if ($amino_acid_3_A eq $amino_acid_3_G) {
					#$AG_type = 'S';
					$S_diffs += $AG;
				} else {
					#$AG_type = 'N';
					$N_diffs += $AG;
				}
				
				if ($amino_acid_3_A eq $amino_acid_3_T) {
					#$AT_type = 'S';
					$S_diffs += $AT;
				} else {
					#$AT_type = 'N';
					$N_diffs += $AT;
				}
				
				if ($amino_acid_3_C eq $amino_acid_3_G) {
					#$CG_type = 'S';
					$S_diffs += $CG;
				} else {
					#$CG_type = 'N';
					$N_diffs += $CG;
				}
				
				if ($amino_acid_3_C eq $amino_acid_3_T) {
					#$CT_type = 'S';
					$S_diffs += $CT;
				} else {
					#$CT_type = 'N';
					$N_diffs += $CT;
				}
				
				if ($amino_acid_3_G eq $amino_acid_3_T) {
					#$GT_type = 'S';
					$S_diffs += $GT;
				} else {
					#$GT_type = 'N';
					$N_diffs += $GT;
				}
				
				$N_diffs_per_comparison += ($N_diffs / $num_pairwise_comparisons);
				$S_diffs_per_comparison += ($S_diffs / $num_pairwise_comparisons);
				
				## AVERAGE dN and dS versus REFERENCE stuff
				my $ref_nt = $hh_product_position_info{$curr_product}->{$site_pos_3}->{reference};
				if (($num_3_A != 0) && ($ref_nt ne 'A')) {
					if($amino_acid_3_A eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_3_A;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_3_A ne 'TAA') && ($codon_3_A ne 'TAG') && ($codon_3_A ne 'TGA') &&
						($codon_3_A ne 'UAA') && ($codon_3_A ne 'UAG') && ($codon_3_A ne 'UGA')) {
						$this_codon_num_N_diffs += $num_3_A;
					}
				}
				
				if (($num_3_C != 0) && ($ref_nt ne 'C')) {
					if($amino_acid_3_C eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_3_C;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_3_C ne 'TAA') && ($codon_3_C ne 'TAG') && ($codon_3_C ne 'TGA') &&
						($codon_3_C ne 'UAA') && ($codon_3_C ne 'UAG') && ($codon_3_C ne 'UGA')) {
						$this_codon_num_N_diffs += $num_3_C;
					}
				}
				
				if (($num_3_G != 0) && ($ref_nt ne 'G')) {
					if($amino_acid_3_G eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_3_G;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_3_G ne 'TAA') && ($codon_3_G ne 'TAG') && ($codon_3_G ne 'TGA') &&
						($codon_3_G ne 'UAA') && ($codon_3_G ne 'UAG') && ($codon_3_G ne 'UGA')) {
						$this_codon_num_N_diffs += $num_3_G;
					}
				}
				
				if (($num_3_T != 0) && ($ref_nt ne 'T')) {
					if($amino_acid_3_T eq $ref_amino_acid) {
						$this_codon_num_S_diffs += $num_3_T;
					} elsif(($curr_codon ne 'TAA') && ($curr_codon ne 'TAG') && ($curr_codon ne 'TGA') &&
						($curr_codon ne 'UAA') && ($curr_codon ne 'UAG') && ($curr_codon ne 'UGA') &&
						($codon_3_T ne 'TAA') && ($codon_3_T ne 'TAG') && ($codon_3_T ne 'TGA') &&
						($codon_3_T ne 'UAA') && ($codon_3_T ne 'UAG') && ($codon_3_T ne 'UGA')) {
						$this_codon_num_N_diffs += $num_3_T;
					}
				}
				$codon_cov_sum += $site_3_cov;
				$codon_cov_denom_sum += 1;
				## END AVERAGE dN and dS versus REFERENCE stuff
				
				####
				## GENE DIVERSITY stuff
				my $this_site_category;
				my $this_site_N;
				my $this_site_S;
				my $this_site_sq_var_freq_sum;
				if ($num_3_A != 0) {
					if ($ref_nt ne 'A') {
						my $A_freq = ($num_3_A / $site_3_cov);
						my $A_freq_sq = ($A_freq ** 2);
						$this_site_sq_var_freq_sum += $A_freq_sq;
						if($amino_acid_3_A eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $A_freq = ($num_3_A / $site_3_cov);
						my $A_freq_sq = ($A_freq ** 2);
						$this_site_sq_var_freq_sum += $A_freq_sq;
					}
				}
				
				if ($num_3_C != 0) {
					if ($ref_nt ne 'C') {
						my $C_freq = ($num_3_C / $site_3_cov);
						my $C_freq_sq = ($C_freq ** 2);
						$this_site_sq_var_freq_sum += $C_freq_sq;
						if($amino_acid_3_C eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $C_freq = ($num_3_C / $site_3_cov);
						my $C_freq_sq = ($C_freq ** 2);
						$this_site_sq_var_freq_sum += $C_freq_sq;
					}
				}
				
				if ($num_3_G != 0) {
					if ($ref_nt ne 'G') {
						my $G_freq = ($num_3_G / $site_3_cov);
						my $G_freq_sq = ($G_freq ** 2);
						$this_site_sq_var_freq_sum += $G_freq_sq;
						if($amino_acid_3_G eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $G_freq = ($num_3_G / $site_3_cov);
						my $G_freq_sq = ($G_freq ** 2);
						$this_site_sq_var_freq_sum += $G_freq_sq;
					}
				}
				
				if ($num_3_T != 0) {
					if ($ref_nt ne 'T') {
						my $T_freq = ($num_3_T / $site_3_cov);
						my $T_freq_sq = ($T_freq ** 2);
						$this_site_sq_var_freq_sum += $T_freq_sq;
						if($amino_acid_3_T eq $ref_amino_acid) {
							$this_site_S = 1;
						} else {
							$this_site_N = 1;
						}
					} else {
						my $T_freq = ($num_3_T / $site_3_cov);
						my $T_freq_sq = ($T_freq ** 2);
						$this_site_sq_var_freq_sum += $T_freq_sq;
					}
				}
				
				my $this_site_category_sum = ($this_site_N + $this_site_S);
				if($this_site_category_sum == 0) {
					$site3_siteb_category = 0;
				} elsif ($this_site_category_sum == 2) {
					$site3_siteb_category = 3;
				} else {
					if($this_site_S == 1) {
						$site3_siteb_category = 1;
					} else {
						$site3_siteb_category = 2;
					}
				}
				
				my $this_site_gene_diversity = (1 - $this_site_sq_var_freq_sum);
				$this_codon_site3_gene_diversity = $this_site_gene_diversity;
				$this_codon_site3_polymorphic = 1;
				## END GENE DIVERSITY stuff
				####
				
			} else {
				my @num_sites_arr = &get_number_of_sites($curr_codon,3);
				$N_sites += $num_sites_arr[0];
				$S_sites += $num_sites_arr[1];
				
				$N_diffs_per_comparison += 0;
				$S_diffs_per_comparison += 0;
				
				# Find out what nucleotide it is and give it a num of 1
				if ($seq_by_index_arr[$site_pos_3-1] eq 'A') {
					$num_3_A = 1;
					$site_3_num_hash{'A'} = $num_3_A;
				} elsif ($seq_by_index_arr[$site_pos_3-1] eq 'C') {
					$num_3_C = 1;
					$site_3_num_hash{'C'} = $num_3_C;
				} elsif ($seq_by_index_arr[$site_pos_3-1] eq 'G') {
					$num_3_G = 1;
					$site_3_num_hash{'G'} = $num_3_G;
				} elsif ($seq_by_index_arr[$site_pos_3-1] eq 'T') {
					$num_3_T = 1;
					$site_3_num_hash{'T'} = $num_3_T;
				}
				
				# Set this site's cov to 1, allowing site/cov to be 100%
				$site_3_cov = 1;
				
				####
				## GENE DIVERSITY stuff
				$this_codon_site3_gene_diversity = 0;
				$this_codon_site3_polymorphic = 0;
				$site3_siteb_category = 0;
				$site3_codonb_category = 0;
				## END GENE DIVERSITY stuff
				####
			}
			#print "in $curr_product, site $site_pos_3, the nt is ".$seq_by_pos_hash{$site_pos_3}.", and site_3_num_hash A is ". $site_3_num_hash{'A'} . "\n";
			
			#print "In: $curr_product, site $site_pos_3, the nt is ".$seq_by_pos_hash{$site_pos_3}.
			#	", and A=". $site_3_num_hash{'A'} . " / C=".$site_3_num_hash{'C'} .
			#	" / G=".$site_3_num_hash{'G'} ." / T=".$site_3_num_hash{'T'}."\n";
			
			
			#############################################################################
			#############################################################################
			#####################  SUMMATION FOR CODON  #################################
			#############################################################################
			#############################################################################
			
			
			#if($progress_period_count < 80) {
			#	print ".";
			#	$progress_period_count++;
			#} else {
			#	print "\n";
			#	$progress_period_count = 0;
			#}
			
			# These will be SITE-BASED
			my $pi_N;
			my $pi_S;
			
			my $tester_site = -1; #2077
			my $tester_product = 'PB2';
			my $tester_report = 'Virus_36_subsampled-0.001%.txt';
			
			# ABOVE, if it's not present (not a variant), give it a coverage of 1 and a NUM of 1
		
			my %this_codon_freq_hash; # contains {NNN,NNN,NNN,...}=>PROPORTION
			#my %this_codon_num_hash; # contain same, but with numbers
			my $stop_codon_variants_possible = 0;
			
			my $codon_based_freq_sum = 0;
			
			# Build all possible codons at this site, along with their relative frequencies
			foreach my $first_nt (keys %site_1_num_hash) {
				my $working_site_1_num = $site_1_num_hash{$first_nt};
				my $working_codon_1 = $first_nt;
			
				foreach my $second_nt (keys %site_2_num_hash) {
					my $working_site_2_num = $site_2_num_hash{$second_nt};
					my $working_codon_2 = $working_codon_1 . $second_nt;

					foreach my $third_nt (keys %site_3_num_hash) {
						my $working_site_3_num = $site_3_num_hash{$third_nt};
						my $working_codon_3 = $working_codon_2 . $third_nt;
						
						#print "Working codon at site $curr_site: $working_codon_3\n";
							
						my $nt_1_freq = ($working_site_1_num / $site_1_cov);
						my $nt_2_freq = ($working_site_2_num / $site_2_cov);
						my $nt_3_freq = ($working_site_3_num / $site_3_cov);
	
						$this_codon_freq_hash{$working_codon_3} = ($nt_1_freq * $nt_2_freq * $nt_3_freq);
						$codon_based_freq_sum += ($nt_1_freq * $nt_2_freq * $nt_3_freq);
						
						if ($curr_site == $tester_site && $curr_product eq $tester_product && $curr_snp_report_name eq $tester_report) {
							print "\nIn codon $working_codon_3, \nNUM nt1: $working_site_1_num || COV1: $site_1_cov || freq: $nt_1_freq\n".
								"NUM nt2: $working_site_2_num || COV2: $site_2_cov || freq: $nt_2_freq\n".
								"NUM nt3: $working_site_3_num || COV3: $site_3_cov || freq: $nt_3_freq\n".
								"CODON FREQ: ".$this_codon_freq_hash{$working_codon_3};
						}
						
						if (($working_codon_3 eq 'TAA') || ($working_codon_3 eq 'TAG') || ($working_codon_3 eq 'TGA') || 
							($working_codon_3 eq 'UAA') || ($working_codon_3 eq 'UAG') || ($working_codon_3 eq 'UGA')) {
							$stop_codon_variants_possible += 1;
						}
					}
				}
			}
			
			# Now %this_codon_freq_hash has keys of all the possible codons at this
			# position, with values as their frequencies. THIS INCLUDES STOP CODONS.
			# We have decided NOT to renormalize below. When we do renormalize, the number
			# of sites for codons that are predominantly STOP returns to 3. This is
			# most consistent with Nei-Gojobori.
			
			# Coverage is set to 1 if there are no variants at the site. If the actual 
			# coverage is ever 1, there's no way to tell pairwise DIFFERENCES anyway
			my $this_codon_avg_cov;		
			my $site_1_cov_to_add;
			my $site_2_cov_to_add;
			my $site_3_cov_to_add;
			my $num_covered_sites = 0;
			
			if ($site_1_cov > 1) { # i.e., NOT 1
				$num_covered_sites++;
				$site_1_cov_to_add = $site_1_cov;
			} else {
				$site_1_cov_to_add = 0;
			}
			
			if ($site_2_cov > 1) { # i.e., NOT 1
				$num_covered_sites++;
				$site_2_cov_to_add = $site_2_cov;
			} else {
				$site_2_cov_to_add = 0;
			}
			
			if ($site_3_cov > 1) { # i.e., NOT 1
				$num_covered_sites++;
				$site_3_cov_to_add = $site_3_cov;
			} else {
				$site_3_cov_to_add = 0;
			}
			
			if ($num_covered_sites > 0) {
				my $avg_cov = (($site_1_cov_to_add + $site_2_cov_to_add + $site_3_cov_to_add) / $num_covered_sites);
				#my $this_codon_num = ($nt_1_freq * $nt_2_freq * $nt_3_freq);#*$cov;
				#$this_codon_num_hash{$working_codon_3} = $avg_cov;
				$this_codon_avg_cov = $avg_cov;
			} else {
				#$this_codon_num_hash{$working_codon_3} = 0;
				$this_codon_avg_cov = 0;
			}
			
			my $codon_based_N_sites = 0;
			my $codon_based_S_sites = 0;
			my $seen_possible_stop = 0;
			my $seen_possible_nonstop = 0;
			
			# NUMBER OF SITES
			foreach my $possible_codon (keys %this_codon_freq_hash) { # This will
				# return 0 for all STOP codon sites, so we don't have to test it here
				my @num_sites_1_arr = &get_number_of_sites($possible_codon,1);
				my @num_sites_2_arr = &get_number_of_sites($possible_codon,2);
				my @num_sites_3_arr = &get_number_of_sites($possible_codon,3);
				
				my $possible_codon_num_N_sites = ($num_sites_1_arr[0] + $num_sites_2_arr[0] + $num_sites_3_arr[0]);
				my $possible_codon_num_S_sites = ($num_sites_1_arr[1] + $num_sites_2_arr[1] + $num_sites_3_arr[1]);
				
				# even if it IS a STOP codon (in which case it has 0 sites)
				$codon_based_N_sites += ($this_codon_freq_hash{$possible_codon} * $possible_codon_num_N_sites);
				$codon_based_S_sites += ($this_codon_freq_hash{$possible_codon} * $possible_codon_num_S_sites);
				
#				if(($seen_possible_stop == 0) && (($possible_codon eq 'TAA') || ($possible_codon eq 'TAG') || ($possible_codon eq 'TGA') || 
#				($possible_codon eq 'UAA') || ($possible_codon eq 'UAG') || ($possible_codon eq 'UGA'))) {
#					$seen_possible_stop = 1;
#				}
#				
#				if(($seen_possible_nonstop == 0) && (($possible_codon ne 'TAA') && ($possible_codon ne 'TAG') && ($possible_codon ne 'TGA') && 
#				($possible_codon ne 'UAA') && ($possible_codon ne 'UAG') && ($possible_codon ne 'UGA'))) {
#					$seen_possible_nonstop = 1;
#				}
			}
			
#			# RENORMALIZE CODON-based NUMBER of SITES, with 0 sites for exclusive STOP codons
#			if(($seen_possible_nonstop == 0) && ($seen_possible_stop == 1)) { # ONLY STOP codons seen
#				$codon_based_N_sites = 0;
#				$codon_based_S_sites = 0;
#			} else {
#				my $codon_based_normalization_sum = ($codon_based_N_sites + $codon_based_S_sites);
#				$codon_based_N_sites = (3 * $codon_based_N_sites / $codon_based_normalization_sum);
#				$codon_based_S_sites = (3 * $codon_based_S_sites / $codon_based_normalization_sum);
#			}
			
			# PERFECT!
			#if ($curr_snp_report_name eq $tester_report && $curr_product eq $tester_product && $curr_site == $tester_site) {
			#	print "\n".$this_codon_freq_hash{'ATG'}."\n";
			#	print "\nN sites here are $tester_N_sites\nS sites here are $tester_S_sites\n";
			#}
			
			# If the REFERENCE is a STOP codon, it should have 0 sites; otherwise normal
			my $ref_based_N_sites;
			my $ref_based_S_sites;
			
			if(($curr_codon eq 'TAA') || ($curr_codon eq 'TAG') || ($curr_codon eq 'TGA') || 
				($curr_codon eq 'UAA') || ($curr_codon eq 'UAG') || ($curr_codon eq 'UGA')) {
				
				$ref_based_N_sites = 0;
				$ref_based_S_sites = 0;
				
			} else {
				my @ref_based_sites_1 = &get_number_of_sites($curr_codon,1);
				my @ref_based_sites_2 = &get_number_of_sites($curr_codon,2);
				my @ref_based_sites_3 = &get_number_of_sites($curr_codon,3);
				
				$ref_based_N_sites = $ref_based_sites_1[0] + $ref_based_sites_2[0] + 
					$ref_based_sites_3[0];
				$ref_based_S_sites = $ref_based_sites_1[1] + $ref_based_sites_2[1] + 
					$ref_based_sites_3[1];
			}
			
			# These are SITE-BASED
			if ($N_sites > 0) {
				$pi_N = ($N_diffs_per_comparison / $N_sites);
			} else {
				$pi_N = "*";
			}
			
			if ($S_sites > 0) {
				$pi_S = ($S_diffs_per_comparison / $S_sites);
			} else {
				$pi_S = "*";
			}
			
			my $pi_N_over_S; #($pi_N / $pi_S);
			
			if ($pi_N eq "*" || $pi_S eq "*" || $N_sites == 0 || $S_sites == 0 || $pi_S == 0) {
				$pi_N_over_S = "*";
			} else {
				$pi_N_over_S = ($pi_N / $pi_S);
			}
			
			# CODON-BASED DIFFERENCES
			# Use our amazing new subroutine, return_avg_diffs, to return a BETTER
			# estime of the number of differences. We will also use
			# %this_codon_freq_hash, whose KEYS are the codons are this site,
			# and whose VALUES are their frequencies (proportions).
			# We must multiply the frequencies of all pairwise comparisons of codons,
			# and divide those products by the sum of the products for the frequency
			# (weight) of the pair, and its contribution to the number of differences
			# at this site.
			
			my $avg_N_diffs = 0;
			my $avg_S_diffs = 0;
			
			my @this_site_codons = sort (keys %this_codon_freq_hash);
			#my @this_site_codons = sort (keys %this_codon_num_hash);
			my $num_codons_here = scalar (@this_site_codons); 
			my $this_codon_pw_comps = ((($this_codon_avg_cov*$this_codon_avg_cov) - 
				$this_codon_avg_cov) / 2); # this will INCLUDE comparisons involving
				# STOP codons, which will later simply not contribute to numbers of
				# differences.
			my $this_codon_pw_comps_minus_STOP = $this_codon_pw_comps;
			
			my %product_freq_hh;
			my %product_num_hh;
			my $product_freq_sum;
			
			# Get the sum of the frequency products, and store each individual product,
			# INCLUDING comparisons between identical codons (there will be many, and the
			# sum needs to include these for proper weights!)
			# That is, getting the proportion of comparisons which are, e.g., 
			# ACG versus ACA, etc.
			if ($num_codons_here > 1) { # e.g., (ATT,ACT) -- two codons here
				for (my $i=0; $i<$num_codons_here; $i++) {
				
					my $codon1 = $this_site_codons[$i];
					
					for (my $j=$i; $j<$num_codons_here; $j++) {
						my $codon2 = $this_site_codons[$j];
						my $freq1 = $this_codon_freq_hash{$codon1};
						my $freq2 = $this_codon_freq_hash{$codon2};
						my $freq_product = ($freq1 * $freq2);
						$product_freq_sum += $freq_product;
						$product_freq_hh{$codon1}->{$codon2} = $freq_product;
						
						#print "freq1 is: $freq1\nfreq2 is: $freq2\n";
						
						if ($curr_site == $tester_site && $curr_product eq $tester_product && $curr_snp_report_name eq $tester_report) {
							print "\n".
							"CODON 1: $codon1\nCODON 2: $codon2\n".
							"FREQ1: $freq1\nFREQ2: $freq2\n";
						}
						
						if ($i != $j) {
							my $num1 = $freq1 * $this_codon_avg_cov;
							my $num2 = $freq2 * $this_codon_avg_cov;
							my $num_product = $num1 * $num2;
							$product_num_hh{$codon1}->{$codon2} = $num_product;
							
							if ($curr_site == $tester_site && $curr_product eq $tester_product && $curr_snp_report_name eq $tester_report) {
								print "NUM1: $num1\nNUM2: $num2\n\n".
								"";
							}
						}
					}
				}
			} # else do nothing; there are 0 differences
			
			# CODON-BASED DIFFERENCES
			my $new_avg_N_diffs = 0;
			my $new_avg_S_diffs = 0;

			# Now do the same loop, but now add to avg_N_diffs and avg_S_diffs, weighting
			# by the product frequency over the sum of such frequencies
			# That is, CODON-BASED DIFFERENCES
			if ($num_codons_here > 1) {
				for (my $i=0; $i<$num_codons_here; $i++) {
					my $codon1 = $this_site_codons[$i];
					
					for (my $j=$i; $j<$num_codons_here; $j++) {
						my $codon2 = $this_site_codons[$j];
						
						#print "The product_freq_sum is $product_freq_sum\n";
						
						my $weight = ($product_freq_hh{$codon1}->{$codon2} / $product_freq_sum);
						
						if ($i ne $j) {
							my @diffs_arr = &return_avg_diffs($codon1,$codon2); # Here we
								# return ZERO (0) DIFFERENCES if a STOP CODON is involved
							my $curr_N_diffs = $diffs_arr[0];
							my $curr_S_diffs = $diffs_arr[1];
							
							my $N_to_add = ($curr_N_diffs * $weight);
							my $S_to_add = ($curr_S_diffs * $weight);
							
							my $new_weight = ($product_num_hh{$codon1}->{$codon2} / $this_codon_pw_comps);
							my $new_N_to_add = ($curr_N_diffs * $new_weight);
							my $new_S_to_add = ($curr_S_diffs * $new_weight);
							
							$avg_N_diffs += $N_to_add;
							$avg_S_diffs += $S_to_add;
							
							$new_avg_N_diffs += $new_N_to_add;
							$new_avg_S_diffs += $new_S_to_add;
							
							#print "The number of N diffs are $curr_N_diffs\nThe number of ".
							#	"S diffs are $curr_S_diffs\n$avg_N_diffs are being added to N\n".
							#	"$avg_S_diffs are being added to S\n\n\n";
							if ($curr_site == $tester_site && $curr_product eq $tester_product && $curr_snp_report_name eq $tester_report) {
								print "".
								"CODON 1: $codon1\nCODON 2: $codon2\n".#"NEW WT: $new_weight\n".
								"".
								"NEW WT: $new_weight\nNEW N ADD: $new_N_to_add\nNEW S ADD: $new_S_to_add\n\n";
							}
							
						} #else there are no differences; $i and $j are identical; add nothing
					}
					
					if ($curr_site == $tester_site && $curr_product eq $tester_product && $curr_snp_report_name eq $tester_report) {
						print "AVG COV: $this_codon_avg_cov\nNUM PW COMPS: $this_codon_pw_comps\n\n";
					}
					
				}
			} # else do nothing; there are 0 differences
			
			## AVERAGE dN and dS versus REFERENCE stuff
			my $avg_dN_vs_ref;
			my $avg_dS_vs_ref;
			my $codon_test_avg_cov = 0;
			my $this_codon_mean_N_diffs_vs_ref = 0;
			my $this_codon_mean_S_diffs_vs_ref = 0;
			if($codon_cov_denom_sum > 0) { # Then we've seen at least one variant
				my $this_codon_avg_cov = ($codon_cov_sum / $codon_cov_denom_sum);
				my $this_codon_avg_N_sites = ($this_codon_avg_cov * $codon_based_N_sites);
				my $this_codon_avg_S_sites = ($this_codon_avg_cov * $codon_based_S_sites);
				$codon_test_avg_cov = $this_codon_avg_cov;
				$this_codon_mean_N_diffs_vs_ref = ($this_codon_num_N_diffs / 
					$this_codon_avg_cov);
				$this_codon_mean_S_diffs_vs_ref = ($this_codon_num_S_diffs / 
					$this_codon_avg_cov);
				$hh_product_position_info{$curr_product}->{mean_N_diffs_ref} +=
					$this_codon_mean_N_diffs_vs_ref;
				$hh_product_position_info{$curr_product}->{mean_S_diffs_ref} +=
					$this_codon_mean_S_diffs_vs_ref;
				
				if($this_codon_avg_N_sites > 0) {
					$avg_dN_vs_ref = ($this_codon_num_N_diffs / $this_codon_avg_N_sites);
				} else {
					$avg_dN_vs_ref = '*';
				}
				
				if($this_codon_avg_S_sites > 0) {
					$avg_dS_vs_ref = ($this_codon_num_S_diffs / $this_codon_avg_S_sites);
				} else {
					$avg_dS_vs_ref = '*';
				}
			} else { # There are NO variants
				if($codon_based_N_sites > 0) {
					$avg_dN_vs_ref = 0;
				} else {
					$avg_dN_vs_ref = '*';
				}
				
				if($codon_based_S_sites > 0) {
					$avg_dS_vs_ref = 0;
				} else {
					$avg_dS_vs_ref = '*';
				}
			}
			##
			
			###########################
			#### OVERLAPPING ORF NUCLEOTIDE SITES STUFF
			# Figure out if these sites overlap other ORFs OVERLAP OVERLAPPING
			my $site1_num_overlap = 0;
			my $site2_num_overlap = 0;
			my $site3_num_overlap = 0;
			
			my $site1_overlap = 0;
			my $site2_overlap = 0;
			my $site3_overlap = 0;
			
			# Recorded NUM OVERLAPPING PRODUCTS: this will give us how many products are at the site
			# THIS INCLUDES the COMPLEMENT '-' strand
			# New segments approach 
			# SITE 1
			if($hh_nc_position_info{$site_pos_1}->{coding} > 1) {
				$site1_overlap = 1;
				$site1_num_overlap = $hh_nc_position_info{$site_pos_1}->{coding} - 1;
			}
			
			# SITE 2
			if($hh_nc_position_info{$site_pos_2}->{coding} > 1) {
				$site2_overlap = 1;
				$site2_num_overlap = $hh_nc_position_info{$site_pos_2}->{coding} - 1;
			}
			
			# SITE 3
			if($hh_nc_position_info{$site_pos_3}->{coding} > 1) {
				$site3_overlap = 1;
				$site3_num_overlap = $hh_nc_position_info{$site_pos_3}->{coding} - 1;
			}
			
			my $codon_overlap_nts = ($site1_overlap+$site2_overlap+$site3_overlap);
			
			#if($sitebasedmode == 1) {
			#	$ntd_headers_to_print .= "πN/πS\t";
			#}
			
			if($sepfiles == 1) { # PRINT TO INDIVIDUAL NUCLEOTIDE DIVERSITY FILE
				open(OUTFILE_NT_DIV_IND,">>$new_file_prefix\_results\.txt");
				# PRINT HEADERS TO INDIVIDUAL NUCLEOTIDE DIVERSITY FILE, IF IT'S THE FIRST CDS
				if ($product_counter == 0) { 
					my $ntd_headers_to_print_ind = "file\tproduct\tsite\tcodon\tnum_overlap_ORF_nts\t".
						"mean_nonsyn_diffs\tmean_syn_diffs\t";
					if($sitebasedmode == 1) {
						$ntd_headers_to_print_ind .= "mean_nonsyn_diffs_site_based\tmean_syn_diffs_site_based\t";
					}
					$ntd_headers_to_print_ind .= "nonsyn_sites\tsyn_sites\t";
					if($sitebasedmode == 1) {
						$ntd_headers_to_print_ind .= "nonsyn_sites_site_based\tsyn_sites_site_based\t";
					}
					$ntd_headers_to_print_ind .= "nonsyn_sites_ref\tsyn_sites_ref\t";
					#$ntd_headers_to_print_ind .= "piN\tpiS\t";
					if($ratiomode == 1) {
						$ntd_headers_to_print_ind .= "piN/piS\t";
					}
					#$ntd_headers_to_print_ind .= "mean_dN_vs_ref\tmean_dS_vs_ref\n"; # \tAverage_cov
					$ntd_headers_to_print_ind .= "mean_nonsyn_diffs_vs_ref\tmean_syn_diffs_vs_ref\t".
						"mean_gdiv\tmean_nonsyn_gdiv\tmean_syn_gdiv\n";
				
					$product_counter += 1; # Increment so header is only printed once
					
					print OUTFILE_NT_DIV_IND "$ntd_headers_to_print_ind";
				}
				
				# PRINT DATA LINE TO INDIVIDUAL NUCLEOTIDE DIVERSITY FILE
				my $ntd_line_to_print_ind = "$file_nm\t$curr_product\t$curr_site\t".
					"$curr_codon\t$codon_overlap_nts\t$new_avg_N_diffs\t$new_avg_S_diffs\t";
				if($sitebasedmode == 1) {
					$ntd_line_to_print_ind .= "$N_diffs_per_comparison\t$S_diffs_per_comparison\t";
				}	
				$ntd_line_to_print_ind .= "$codon_based_N_sites\t$codon_based_S_sites\t";
				if($sitebasedmode == 1) {
					$ntd_line_to_print_ind .= "$N_sites\t$S_sites\t";
				}
				$ntd_line_to_print_ind .= "$ref_based_N_sites\t$ref_based_S_sites\t";
				#$ntd_line_to_print_ind .= "$pi_N\t$pi_S\t";
				if($ratiomode == 1) {
					$ntd_line_to_print_ind .= "$pi_N_over_S\t"; # SITE-BASED
				}
				#$ntd_line_to_print_ind .= "$avg_dN_vs_ref\t$avg_dS_vs_ref\n"; # \t$codon_test_avg_cov
				#$ntd_line_to_print_ind .= "$this_codon_num_N_diffs\t".
				#	"$this_codon_num_S_diffs\n";
				$ntd_line_to_print_ind .= "$this_codon_mean_N_diffs_vs_ref\t".
					"$this_codon_mean_S_diffs_vs_ref\t";
				
				print OUTFILE_NT_DIV_IND "$ntd_line_to_print_ind";
				close OUTFILE_NT_DIV_IND;
			}
			
			# PRINT DATA LINE TO COMPILED NUCLEOTIDE DIVERSITY FILE
			open(OUTFILE_NT_DIV,">>codon\_results\.txt");
			
			my $ntd_line_to_print = "$file_nm\t$curr_product\t$curr_site\t".
				"$curr_codon\t$codon_overlap_nts\t$new_avg_N_diffs\t$new_avg_S_diffs\t";
			if($sitebasedmode == 1) {
				$ntd_line_to_print .= "$N_diffs_per_comparison\t$S_diffs_per_comparison\t";
			}
			$ntd_line_to_print .= "$codon_based_N_sites\t$codon_based_S_sites\t";
			if($sitebasedmode == 1) {
				$ntd_line_to_print .= "$N_sites\t$S_sites\t";
			}
			$ntd_line_to_print .= "$ref_based_N_sites\t$ref_based_S_sites\t";
			#$ntd_line_to_print .= "$pi_N\t$pi_S\t";
			if($ratiomode == 1) {
				$ntd_line_to_print .= "$pi_N_over_S\t"; # SITE-BASED
			}
			#$ntd_line_to_print .= "$avg_dN_vs_ref\t$avg_dS_vs_ref\n"; # \t$codon_test_avg_cov
			#$ntd_line_to_print .= "$this_codon_num_N_diffs\t$this_codon_num_S_diffs\n";
			$ntd_line_to_print .= "$this_codon_mean_N_diffs_vs_ref\t".
				"$this_codon_mean_S_diffs_vs_ref\t";
			
			print OUTFILE_NT_DIV "$ntd_line_to_print";
			close OUTFILE_NT_DIV;
			
			# Finally, add PRODUCT TOTALS data to the %hh_product_position_info
			$hh_product_position_info{$curr_product}->{N_diffs_codonb} += $N_diffs_per_comparison;
			$hh_product_position_info{$curr_product}->{S_diffs_codonb} += $S_diffs_per_comparison;
			$hh_product_position_info{$curr_product}->{N_sites_codonb} += $codon_based_N_sites;
			$hh_product_position_info{$curr_product}->{S_sites_codonb} += $codon_based_S_sites;
			$hh_product_position_info{$curr_product}->{N_diffs_vs_ref} += $this_codon_num_N_diffs;
			$hh_product_position_info{$curr_product}->{S_diffs_vs_ref} += $this_codon_num_S_diffs;
			$hh_product_position_info{$curr_product}->{mean_N_diffs_vs_ref} += $this_codon_mean_N_diffs_vs_ref;
			$hh_product_position_info{$curr_product}->{mean_S_diffs_vs_ref} += $this_codon_mean_S_diffs_vs_ref;
			
			
			####
			## GENE DIVERSITY stuff
			my $curr_nt_1 = substr($curr_codon,0,1);
			my $curr_nt_2 = substr($curr_codon,1,1);
			my $curr_nt_3 = substr($curr_codon,2,1);
			
			# DETERMINE CODON-BASED SITE CATEGORIES
			# Categories: 0=no variants; 1=S; 2=N; 3=both/ambiguous
			# All codons are stored, sorted, in @this_site_codons
			# The nucleotides present at each site are stored in %site_1_num_hash,
			# %site_2_num_hash, and %site_3_num_hash, with nucleotides as keys
			# and the number of those nucleotides as values. Will be 0 if not present.
			
			my $this_site1_S = 0;
			my $this_site1_N = 0;
			my $this_site2_S = 0;
			my $this_site2_N = 0;
			my $this_site3_S = 0;
			my $this_site3_N = 0;
			
			foreach my $possible_codon (@this_site_codons) {
				my $curr_possible_amino_acid = &get_amino_acid($possible_codon);
				
				my $site1_given_nt = substr($possible_codon,0,1);
				foreach my $possible_site1_nt (keys %site_1_num_hash) {
					if(($site_1_num_hash{$possible_site1_nt} > 0) && ($site1_given_nt ne $possible_site1_nt)) {
						my $curr_variant_codon = "" . $possible_site1_nt . substr($possible_codon,1,2);
						my $curr_variant_amino_acid = &get_amino_acid($curr_variant_codon);
						
						if($curr_possible_amino_acid eq $curr_variant_amino_acid) { # this is a SYNONYMOUS difference
							$this_site1_S = 1;
							# COULD DO A KIND OF WEIGHTING HERE
						} else { # this is a NONSYNONYMOUS difference
							$this_site1_N = 1;
						}
					}
				}
				
				my $site2_given_nt = substr($possible_codon,1,1);
				foreach my $possible_site2_nt (keys %site_2_num_hash) {
					if(($site_2_num_hash{$possible_site2_nt} > 0) && ($site2_given_nt ne $possible_site2_nt)) {
						my $curr_variant_codon = "" . substr($possible_codon,0,1) . $possible_site2_nt . substr($possible_codon,2,1);
						my $curr_variant_amino_acid = &get_amino_acid($curr_variant_codon);
						
						if($curr_possible_amino_acid eq $curr_variant_amino_acid) { # this is a SYNONYMOUS difference
							$this_site2_S = 1;
							# COULD DO A KIND OF WEIGHTING HERE
						} else { # this is a NONSYNONYMOUS difference
							$this_site2_N = 1;
						}
					}
				}
				
				my $site3_given_nt = substr($possible_codon,2,1);
				foreach my $possible_site3_nt (keys %site_3_num_hash) {
					if(($site_3_num_hash{$possible_site3_nt} > 0) && ($site3_given_nt ne $possible_site3_nt)) {
						my $curr_variant_codon = "" . substr($possible_codon,0,2) . $possible_site3_nt;
						my $curr_variant_amino_acid = &get_amino_acid($curr_variant_codon);
						
						if($curr_possible_amino_acid eq $curr_variant_amino_acid) { # this is a SYNONYMOUS difference
							$this_site3_S = 1;
							# COULD DO A KIND OF WEIGHTING HERE
						} else { # this is a NONSYNONYMOUS difference
							$this_site3_N = 1;
						}
					}
				}
			}
			
			# CLASSIFY SITES
			my $site1_category_sum = ($this_site1_N + $this_site1_S);
			if($site1_category_sum == 0) {
				$site1_codonb_category = 0;
			} elsif ($site1_category_sum == 2) {
				$site1_codonb_category = 3;
			} else {
				if($this_site1_S == 1) {
					$site1_codonb_category = 1;
				} else {
					$site1_codonb_category = 2;
				}
			}
			
			my $site2_category_sum = ($this_site2_N + $this_site2_S);
			if($site2_category_sum == 0) {
				$site2_codonb_category = 0;
			} elsif ($site2_category_sum == 2) {
				$site2_codonb_category = 3;
			} else {
				if($this_site2_S == 1) {
					$site2_codonb_category = 1;
				} else {
					$site2_codonb_category = 2;
				}
			}
			
			my $site3_category_sum = ($this_site3_N + $this_site3_S);
			if($site3_category_sum == 0) {
				$site3_codonb_category = 0;
			} elsif ($site3_category_sum == 2) {
				$site3_codonb_category = 3;
			} else {
				if($this_site3_S == 1) {
					$site3_codonb_category = 1;
				} else {
					$site3_codonb_category = 2;
				}
			}
			
			# If the reference is a STOP codon, set the site category to NA
			if($curr_codon eq 'TAA' || $curr_codon eq 'TAG' || $curr_codon eq 'TGA' || 
				$curr_codon eq 'UAA' || $curr_codon eq 'UAG' || $curr_codon eq 'UGA') {
					$site1_siteb_category = 'NA';
					$site2_siteb_category = 'NA';
					$site3_siteb_category = 'NA';
					$site1_codonb_category = 'NA';
					$site2_codonb_category = 'NA';
					$site3_codonb_category = 'NA';
			}
			
			# PRINT GENE DIVERSITY DATA TO FILE
			open(OUTFILE_GENE_DIV,">>site\_results\.txt");			
			
			# For codon's FIRST SITE -- print only if there is a polymorphism
			if($this_codon_site1_polymorphic == 1) {
				my $this_cov = $hh_product_position_info{$curr_product}->{$site_pos_1}->{cov};
				my $maj_nt = $hh_product_position_info{$curr_product}->{$site_pos_1}->{maj_nt};
				
				my $siteb_category;
				if ($site1_siteb_category == 1) {
					$siteb_category = 'Synonymous';
				} elsif ($site1_siteb_category == 2) {
					$siteb_category = 'Nonsynonymous';
				} elsif ($site1_siteb_category == 3) {
					$siteb_category = 'Ambiguous';
				} elsif ($site1_siteb_category eq 'NA') {
					$siteb_category = 'NA';
				} else {
					$siteb_category = 'ERROR';
				}
				
				my $codonb_category;
				if ($site1_codonb_category == 1) {
					$codonb_category = 'Synonymous';
				} elsif ($site1_codonb_category == 2) {
					$codonb_category = 'Nonsynonymous';
				} elsif ($site1_codonb_category == 3) {
					$codonb_category = 'Ambiguous';
				} elsif ($site1_codonb_category eq 'NA') {
					$codonb_category = 'NA';
				} else {
					$codonb_category = 'ERROR';
				}
		
				my $pi;
				if($hh_nc_position_info{$site_pos_1}->{pi} > 0) {
					$pi = $hh_nc_position_info{$site_pos_1}->{pi};
				} else {
					$pi = 0;
				}
				
				# Round nucleotide counts to 0 for rounding error to the 9th decimal
#				if($num_1_A < 0.000000001 && $num_1_A > -0.000000001) {
#					$num_1_A = 0;
#				}
#				if($num_1_C < 0.000000001 && $num_1_C > -0.000000001) {
#					$num_1_C = 0;
#				}
#				if($num_1_G < 0.000000001 && $num_1_G > -0.000000001) {
#					$num_1_G = 0;
#				}
#				if($num_1_T < 0.000000001 && $num_1_T > -0.000000001) {
#					$num_1_T = 0;
#				}
				
				print OUTFILE_GENE_DIV "$file_nm\t$curr_product\t$site_pos_1\t".
					"$curr_nt_1\t$maj_nt\t".
					"1\t$site1_num_overlap\t$curr_site\t$curr_codon\t$pi\t".
					#"$this_codon_site1_polymorphic\t".
					"$this_codon_site1_gene_diversity\t$siteb_category\t$codonb_category\t".
					"$this_cov\t".
					"$num_1_A\t$num_1_C\t$num_1_G\t$num_1_T\n";
			}
				
			# For codon's SECOND SITE -- print only if there is a polymorphism
			if($this_codon_site2_polymorphic == 1) {
				my $this_cov = $hh_product_position_info{$curr_product}->{$site_pos_2}->{cov};
				my $maj_nt = $hh_product_position_info{$curr_product}->{$site_pos_2}->{maj_nt};
				
				my $siteb_category;
				if ($site2_siteb_category == 1) {
					$siteb_category = 'Synonymous';
				} elsif ($site2_siteb_category == 2) {
					$siteb_category = 'Nonsynonymous';
				} elsif ($site2_siteb_category == 3) {
					$siteb_category = 'Ambiguous';
				} elsif ($site2_siteb_category eq 'NA') {
					$siteb_category = 'NA';
				} else {
					$siteb_category = 'ERROR';
				}
				
				my $codonb_category;
				if ($site2_codonb_category == 1) {
					$codonb_category = 'Synonymous';
				} elsif ($site2_codonb_category == 2) {
					$codonb_category = 'Nonsynonymous';
				} elsif ($site2_codonb_category == 3) {
					$codonb_category = 'Ambiguous';
				} elsif ($site2_codonb_category eq 'NA') {
					$codonb_category = 'NA';
				} else {
					$codonb_category = 'ERROR';
				}
			
				my $pi;
				if($hh_nc_position_info{$site_pos_2}->{pi} > 0) {
					$pi = $hh_nc_position_info{$site_pos_2}->{pi};
				} else {
					$pi = 0;
				}
				
				# Round nucleotide counts to 0 for rounding error to the 9th decimal
#				if($num_2_A < 0.000000001 && $num_2_A > -0.000000001) {
#					$num_2_A = 0;
#				}
#				if($num_2_C < 0.000000001 && $num_2_C > -0.000000001) {
#					$num_2_C = 0;
#				}
#				if($num_2_G < 0.000000001 && $num_2_G > -0.000000001) {
#					$num_2_G = 0;
#				}
#				if($num_2_T < 0.000000001 && $num_2_T > -0.000000001) {
#					$num_2_T = 0;
#				}
				
				print OUTFILE_GENE_DIV "$file_nm\t$curr_product\t$site_pos_2\t".
					"$curr_nt_2\t$maj_nt\t".
					"2\t$site2_num_overlap\t$curr_site\t$curr_codon\t$pi\t".
					#"$this_codon_site2_polymorphic\t".
					"$this_codon_site2_gene_diversity\t$siteb_category\t$codonb_category\t".
					"$this_cov\t".
					"$num_2_A\t$num_2_C\t$num_2_G\t$num_2_T\n";
			}
				
			# For codon's THIRD SITE -- print only if there is a polymorphism
			if($this_codon_site3_polymorphic == 1) {
				my $this_cov = $hh_product_position_info{$curr_product}->{$site_pos_3}->{cov};
				my $maj_nt = $hh_product_position_info{$curr_product}->{$site_pos_3}->{maj_nt};
				
				my $siteb_category;
				if ($site3_siteb_category == 1) {
					$siteb_category = 'Synonymous';
				} elsif ($site3_siteb_category == 2) {
					$siteb_category = 'Nonsynonymous';
				} elsif ($site3_siteb_category == 3) {
					$siteb_category = 'Ambiguous';
				} elsif ($site3_siteb_category eq 'NA') {
					$siteb_category = 'NA';
				} else {
					$siteb_category = 'ERROR';
				}
				
				my $codonb_category;
				if ($site3_codonb_category == 1) {
					$codonb_category = 'Synonymous';
				} elsif ($site3_codonb_category == 2) {
					$codonb_category = 'Nonsynonymous';
				} elsif ($site3_codonb_category == 3) {
					$codonb_category = 'Ambiguous';
				} elsif ($site3_codonb_category eq 'NA') {
					$codonb_category = 'NA';
				} else {
					$codonb_category = 'ERROR';
				}
					
				my $pi;
				if($hh_nc_position_info{$site_pos_3}->{pi} > 0) {
					$pi = $hh_nc_position_info{$site_pos_3}->{pi};
				} else {
					$pi = 0;
				}
				
				# Round nucleotide counts to 0 for rounding error to the 9th decimal
#				if($num_3_A < 0.000000001 && $num_3_A > -0.000000001) {
#					$num_3_A = 0;
#				}
#				if($num_3_C < 0.000000001 && $num_3_C > -0.000000001) {
#					$num_3_C = 0;
#				}
#				if($num_3_G < 0.000000001 && $num_3_G > -0.000000001) {
#					$num_3_G = 0;
#				}
#				if($num_3_T < 0.000000001 && $num_3_T > -0.000000001) {
#					$num_3_T = 0;
#				}
				
				print OUTFILE_GENE_DIV "$file_nm\t$curr_product\t$site_pos_3\t".
					"$curr_nt_3\t$maj_nt\t".
					"3\t$site3_num_overlap\t$curr_site\t$curr_codon\t$pi\t".
					#"$this_codon_site3_polymorphic\t".
					"$this_codon_site3_gene_diversity\t$siteb_category\t$codonb_category\t".
					"$this_cov\t".
					"$num_3_A\t$num_3_C\t$num_3_G\t$num_3_T\n";
			}
			
			close OUTFILE_GENE_DIV;
			
			# To get AVERAGE GENE DIVERSITY VALUES for site
			my $this_codon_avg_gene_diversity = 0;
			my $this_codon_avg_gene_diversity_S = 0;
			my $this_codon_avg_gene_diversity_N = 0;
			my $this_codon_numerator_gene_diversity;
			my $this_codon_numerator_gene_diversity_S;
			my $this_codon_numerator_gene_diversity_N;
			my $this_codon_denominator_gene_diversity = 0;
			my $this_codon_denominator_gene_diversity_S = 0;
			my $this_codon_denominator_gene_diversity_N = 0;
			
			if($this_codon_site1_polymorphic == 1) {
				$this_codon_numerator_gene_diversity += $this_codon_site1_gene_diversity;
				$this_codon_denominator_gene_diversity ++;
			}
			if($this_codon_site2_polymorphic == 1) {
				$this_codon_numerator_gene_diversity += $this_codon_site2_gene_diversity;
				$this_codon_denominator_gene_diversity ++;
			}
			if($this_codon_site3_polymorphic == 1) {
				$this_codon_numerator_gene_diversity += $this_codon_site3_gene_diversity;
				$this_codon_denominator_gene_diversity ++;
			}
			
			if($this_codon_denominator_gene_diversity > 0) {
				$this_codon_avg_gene_diversity = ($this_codon_numerator_gene_diversity /
					$this_codon_denominator_gene_diversity);
			}
			
			# Categories: 0=no variants; 1=S; 2=N; 3=both/ambiguous
			
			# Add data to the %hh_product_position_info
			$hh_product_position_info{$curr_product}->{sum_gene_diversity} += 
				($this_codon_site1_gene_diversity + $this_codon_site2_gene_diversity +
				$this_codon_site3_gene_diversity);
			$hh_product_position_info{$curr_product}->{num_poly_sites} += 
				($this_codon_site1_polymorphic + $this_codon_site2_polymorphic +
				$this_codon_site3_polymorphic);
			
			# Add SITE 1 information to sums by category
			if($site1_codonb_category == 1) { # 1 is Synonymous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_S} +=
					$this_codon_site1_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_S} += 1;
				$this_codon_numerator_gene_diversity_S += $this_codon_site1_gene_diversity;
				$this_codon_denominator_gene_diversity_S ++;
			} elsif($site1_codonb_category == 2) { # 2 is Nonsynonymous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_N} +=
					$this_codon_site1_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_N} += 1;
				$this_codon_numerator_gene_diversity_N += $this_codon_site1_gene_diversity;
				$this_codon_denominator_gene_diversity_N ++;
			} elsif($site1_codonb_category == 3) { # 3 is Ambiguous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_A} +=
					$this_codon_site1_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_A} += 1;
			}
			
			# Add SITE 2 information to sums by category
			if($site2_codonb_category == 1) { # 1 is Synonymous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_S} +=
					$this_codon_site2_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_S} += 1;
				$this_codon_numerator_gene_diversity_S += $this_codon_site2_gene_diversity;
				$this_codon_denominator_gene_diversity_S ++;
			} elsif($site2_codonb_category == 2) { # 2 is Nonsynonymous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_N} +=
					$this_codon_site2_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_N} += 1;
				$this_codon_numerator_gene_diversity_N += $this_codon_site2_gene_diversity;
				$this_codon_denominator_gene_diversity_N ++;
			} elsif($site2_codonb_category == 3) { # 3 is Ambiguous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_A} +=
					$this_codon_site2_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_A} += 1;
			}
			
			# Add SITE 3 information to sums by category
			if($site3_codonb_category == 1) { # 1 is Synonymous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_S} +=
					$this_codon_site3_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_S} += 1;
				$this_codon_numerator_gene_diversity_S += $this_codon_site3_gene_diversity;
				$this_codon_denominator_gene_diversity_S ++;
			} elsif($site3_codonb_category == 2) { # 2 is Nonsynonymous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_N} +=
					$this_codon_site3_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_N} += 1;
				$this_codon_numerator_gene_diversity_N += $this_codon_site3_gene_diversity;
				$this_codon_denominator_gene_diversity_N ++;
			} elsif($site3_codonb_category == 3) { # 3 is Ambiguous
				$hh_product_position_info{$curr_product}->{sum_gene_diversity_A} +=
					$this_codon_site3_gene_diversity;
				$hh_product_position_info{$curr_product}->{num_category_A} += 1;
			}
			
			# Calculate average gene diversity by site category if there are variants
			if($this_codon_denominator_gene_diversity_S > 0) {
				$this_codon_avg_gene_diversity_S = ($this_codon_numerator_gene_diversity_S /
					$this_codon_denominator_gene_diversity_S);
			} else {
				$this_codon_avg_gene_diversity_S = '*';
			}
			
			if($this_codon_denominator_gene_diversity_N > 0) {
				$this_codon_avg_gene_diversity_N = ($this_codon_numerator_gene_diversity_N /
					$this_codon_denominator_gene_diversity_N);
			} else {
				$this_codon_avg_gene_diversity_N = '*';
			}
			
			# PRINT AVERAGE GENE DIVERSITY VALUES for this codon to the nucleotide 
			# diversity file(s)
			if($sepfiles == 1) { # PRINT TO INDIVIDUAL NUCLEOTIDE DIVERSITY FILE
				open(OUTFILE_NT_DIV_IND,">>$new_file_prefix\_results\.txt");
				print OUTFILE_NT_DIV_IND "$this_codon_avg_gene_diversity\t".
					"$this_codon_avg_gene_diversity_N\t$this_codon_avg_gene_diversity_S\n";
				close OUTFILE_NT_DIV_IND;
			}
			
			open(OUTFILE_NT_DIV,">>codon\_results\.txt");
			print OUTFILE_NT_DIV "$this_codon_avg_gene_diversity\t".
					"$this_codon_avg_gene_diversity_N\t$this_codon_avg_gene_diversity_S\n";
			close OUTFILE_NT_DIV;
			
			## END GENE DIVERSITY stuff
			####
		} # End going through all codons
		
		####### PRINT PRODUCT TOTALS HERE BEFORE FINISHING WITH CURRENT PRODUCT #########
		## 1 ## NUCLEOTIDE DIVERSITY totals for PRODUCT
		my $sum_N_diffs = $hh_product_position_info{$curr_product}->{N_diffs_codonb};
		my $sum_S_diffs = $hh_product_position_info{$curr_product}->{S_diffs_codonb};
		my $sum_N_sites = $hh_product_position_info{$curr_product}->{N_sites_codonb};
		my $sum_S_sites = $hh_product_position_info{$curr_product}->{S_sites_codonb};
		#my $sum_N_diffs_vs_ref = $hh_product_position_info{$curr_product}->{N_diffs_vs_ref};
		#	# This, EARLIER, must be divided by COVERAGE at each site before being added
		my $sum_mean_N_diffs_vs_ref = $hh_product_position_info{$curr_product}->{mean_N_diffs_vs_ref};
		#my $sum_S_diffs_vs_ref = $hh_product_position_info{$curr_product}->{S_diffs_vs_ref};
		#	# This, EARLIER, must be divided by COVERAGE at each site before being added
		my $sum_mean_S_diffs_vs_ref = $hh_product_position_info{$curr_product}->{mean_S_diffs_vs_ref};
		
		my $product_piN;
		if($sum_N_sites > 0) {
			$product_piN = ($sum_N_diffs/$sum_N_sites);
		} else {
			$product_piN = '*';
		}
		
		my $product_piS;
		if($sum_S_sites > 0) {
			$product_piS = ($sum_S_diffs/$sum_S_sites);
		} else {
			$product_piS = '*';
		}
		
		my $product_dNvREF;
		if($sum_N_sites > 0) {
			#$product_dNvREF = ($sum_N_diffs_vs_ref/$sum_N_sites);
			$product_dNvREF = ($sum_mean_N_diffs_vs_ref/$sum_N_sites);
		} else {
			$product_dNvREF = '*';
		}
		
		my $product_dSvREF;
		if($sum_S_sites > 0) {
			#$product_dSvREF = ($sum_S_diffs_vs_ref/$sum_S_sites);
			$product_dSvREF = ($sum_mean_S_diffs_vs_ref/$sum_S_sites);
		} else {
			$product_dSvREF = '*';
		}
		
		# PRINT TO FILE
		#print "\n\n$file_nm\nPRODUCT: $curr_product\npiN: $product_piN\n".
		#	"piS: $product_piS\ndN vs. Ref: $product_dNvREF\n".
		#	"dS vs. Ref: $product_dSvREF\n\n";
		
		## 2 ## GENE DIVERSITY totals for PRODUCT
		my $sum_gene_diversity = $hh_product_position_info{$curr_product}->{sum_gene_diversity};
		my $num_poly_sites = $hh_product_position_info{$curr_product}->{num_poly_sites};
		my $sum_gene_diversity_S = $hh_product_position_info{$curr_product}->{sum_gene_diversity_S};
		my $sum_gene_diversity_N = $hh_product_position_info{$curr_product}->{sum_gene_diversity_N};
		my $sum_gene_diversity_A = $hh_product_position_info{$curr_product}->{sum_gene_diversity_A};
		my $num_category_S = $hh_product_position_info{$curr_product}->{num_category_S};
		my $num_category_N = $hh_product_position_info{$curr_product}->{num_category_N};
		my $num_category_A = $hh_product_position_info{$curr_product}->{num_category_A};
		
		my $mean_gene_diversity_polymorphic;
		if($num_poly_sites > 0) {
			$mean_gene_diversity_polymorphic = ($sum_gene_diversity / $num_poly_sites);
		} else {
			$mean_gene_diversity_polymorphic = '*';
		}
		
		my $mean_gene_diversity_S;
		if($num_category_S > 0) {
			$mean_gene_diversity_S = ($sum_gene_diversity_S / $num_category_S);
		} else {
			$mean_gene_diversity_S = '*';
		}
		
		my $mean_gene_diversity_N;
		if($num_category_N > 0) {
			$mean_gene_diversity_N = ($sum_gene_diversity_N / $num_category_N);
		} else {
			$mean_gene_diversity_N = '*';
		}
		
		my $mean_gene_diversity_A;
		if($num_category_A > 0) {
			$mean_gene_diversity_A = ($sum_gene_diversity_A / $num_category_A);
		} else {
			$mean_gene_diversity_A = '*';
		}
		
		# PRINT TO FILE
		#print "\n\n$file_nm\nPRODUCT: $curr_product\nMean GD: $mean_gene_diversity_polymorphic\n".
		#	"MEAN GD_S: $mean_gene_diversity_S\nMEAN GD_N: $mean_gene_diversity_N\n".
		#	"MEAN GD_A: $mean_gene_diversity_A\n\n";
		
		### PRODUCT SUMMARY FILE ## WORK IN PROGRESS
		open(PRODUCT_SUMMARY,">>product\_results\.txt");
		print PRODUCT_SUMMARY "$file_nm\t$curr_product\t$sum_N_diffs\t$sum_S_diffs\t".
			#"$sum_N_diffs_vs_ref\t$sum_S_diffs_vs_ref\t".
			"$sum_mean_N_diffs_vs_ref\t$sum_mean_S_diffs_vs_ref\t".
			"$sum_N_sites\t$sum_S_sites\t".
			"$product_piN\t$product_piS\t$product_dNvREF\t$product_dSvREF\t".
			"$mean_gene_diversity_polymorphic\t$mean_gene_diversity_N\t$mean_gene_diversity_S\n";
		close PRODUCT_SUMMARY;
		
		# Go back to working directory, in which SNPGenie_Results resides
		chdir('..');
		
		# ADD TO POPULATION SUMMARY VARIABLES
		$pop_sum_ndiffs += $sum_N_diffs;
		$pop_sum_sdiffs += $sum_S_diffs;
		$pop_sum_ndiffs_vRef += $sum_mean_N_diffs_vs_ref;
		$pop_sum_sdiffs_vRef += $sum_mean_S_diffs_vs_ref;
		$pop_sum_nsites += $sum_N_sites;
		$pop_sum_ssites += $sum_S_sites;
		##$pop_sum_nsites_Ref = 0;
		##$pop_sum_ssites_Ref = 0;
		$pop_sum_gdiv += $sum_gene_diversity;
		#$pop_sum_poly_sites += $num_poly_sites; # Can't do this, because ORFs may overlap
		$pop_sum_Nsites_gdiv += $num_category_N; # These are okay, because may be diff. cats
		$pop_sum_Ssites_gdiv += $num_category_S;
		$pop_sum_Asites_gdiv += $num_category_A;
		$pop_sum_N_gdiv += $sum_gene_diversity_N;
		$pop_sum_S_gdiv += $sum_gene_diversity_S;
		$pop_sum_A_gdiv += $sum_gene_diversity_A;
	} # END LOOP THROUGH ALL PRODUCTS
	print "COMPLETED.\n";
	
	print "\nPerforming final calculations, noncoding overlap analysis, and output... ";
	# DEAL WITH NONCODING SITES for site file and population summary file	
	chdir('SNPGenie_Results');
	open(OUTFILE_GENE_DIV,">>site\_results\.txt");
	foreach my $curr_site (sort {$a <=> $b} (keys %hh_nc_position_info)) { # only poly were stored #comeback; we've got an array already
		if($hh_nc_position_info{$curr_site}->{polymorphic} == 1) { # poly
			if(! $hh_nc_position_info{$curr_site}->{coding}) { # poly-noncoding
				my $ref_nt = $hh_nc_position_info{$curr_site}->{reference};
				my $maj_nt = $hh_nc_position_info{$curr_site}->{maj_nt};
				my $overlapping_ORFs = 0;
				my $pi = $hh_nc_position_info{$curr_site}->{pi};
				my $gdiv = $hh_nc_position_info{$curr_site}->{gdiv};
				my $cov = $hh_nc_position_info{$curr_site}->{cov};
				
				my $A = $hh_nc_position_info{$curr_site}->{A};
				my $C = $hh_nc_position_info{$curr_site}->{C};
				my $G = $hh_nc_position_info{$curr_site}->{G};
				my $T = $hh_nc_position_info{$curr_site}->{T};				
				
				# If REVERSE COMPLEMENT MODE, we include these products as overlapping
				# Go through each product, through each segment
				if($complementmode) {
					foreach my $revcom_product (@curr_compl_products_ordered_by_start) {
						#print "\nProduct: $curr_product, ";
						
						my $num_segments = $hh_compl_position_info{$revcom_product}->{num_segments};				
						
						#COMEBACK this is computationally very intense.
						# For each segment, step through each site, label coding
						for(my $i = 1; $i <= $num_segments; $i++) {
							my $this_start_key = 'start_' . $i;
							my $this_stop_key = 'stop_' . $i;
							
							my $this_start = $hh_compl_position_info{$revcom_product}->{$this_start_key};
							my $this_stop = $hh_compl_position_info{$revcom_product}->{$this_stop_key};
							
							if(($curr_site >= $this_start) && ($curr_site <= $this_stop)) {
								$overlapping_ORFs += 1;
							}
							
							# Assumes that they have been ordered
							if($this_start > $curr_site) {
								last;
							}
						}
					}
				}
				
				# Round nucleotide counts to 0 for rounding error to the 9th decimal
				if($A < 0.000000001 && $A > -0.000000001) {
					$A = 0;
				}
				if($C < 0.000000001 && $C > -0.000000001) {
					$C = 0;
				}
				if($G < 0.000000001 && $G > -0.000000001) {
					$G = 0;
				}
				if($T < 0.000000001 && $T > -0.000000001) {
					$T = 0;
				}
				
				print OUTFILE_GENE_DIV "$file_nm\tnoncoding\t$curr_site\t".
					"$ref_nt\t$maj_nt\t".
					"NA\t$overlapping_ORFs\tNA\tNA\t".
					#"$this_codon_site3_polymorphic\t".
					"$pi\t$gdiv\tnoncoding\tnoncoding\t".
					"$cov\t".
					"$A\t$C\t$G\t$T\n";
			}
		}
	}
	close OUTFILE_GENE_DIV;
	chdir('..');
	
	my $pop_piN;
	if($pop_sum_nsites > 0) {
		$pop_piN = ($pop_sum_ndiffs/$pop_sum_nsites);
	} else {
		$pop_piN = '*';
	}
	
	my $pop_piS;
	if($pop_sum_ssites > 0) {
		$pop_piS = ($pop_sum_sdiffs/$pop_sum_ssites);
	} else {
		$pop_piS = '*';
	}
	
	my $pop_dNvREF;
	if($pop_sum_nsites > 0) {
		#$pop_dNvREF = ($sum_N_diffs_vs_ref/$sum_N_sites);
		$pop_dNvREF = ($pop_sum_ndiffs_vRef/$pop_sum_nsites);
	} else {
		$pop_dNvREF = '*';
	}
	
	my $pop_dSvREF;
	if($pop_sum_ssites > 0) {
		#$pop_dSvREF = ($sum_S_diffs_vs_ref/$sum_S_sites);
		$pop_dSvREF = ($pop_sum_sdiffs_vRef/$pop_sum_ssites);
	} else {
		$pop_dSvREF = '*';
	}
	
	my $num_poly_coding_sites = $h_nc_results{num_poly_coding_sites};
	my $num_poly_nc_sites = $h_nc_results{num_poly_nc_sites};
	$pop_sum_poly_sites = ($num_poly_coding_sites + $num_poly_nc_sites);
	
	my $pop_gdiv_poly;
	if($pop_sum_poly_sites > 0) {
		$pop_gdiv_poly = ($pop_sum_gdiv / $pop_sum_poly_sites);
	} else {
		$pop_gdiv_poly = '*';
	}
	
	my $pop_N_gdiv;
	if($pop_sum_Nsites_gdiv > 0) {
		$pop_N_gdiv = ($pop_sum_N_gdiv / $pop_sum_Nsites_gdiv);
	} else {
		$pop_N_gdiv = '*';
	}
	
	my $pop_S_gdiv;
	if($pop_sum_Ssites_gdiv > 0) {
		$pop_S_gdiv = ($pop_sum_S_gdiv / $pop_sum_Ssites_gdiv);
	} else {
		$pop_S_gdiv = '*';
	}
	
	my $mean_gdiv = $h_nc_results{mean_gdiv};
	my $mean_gdiv_coding_poly = $h_nc_results{mean_gdiv_coding_poly};
	my $mean_gdiv_nc_poly = $h_nc_results{mean_gdiv_nc_poly};
	
#	print "\npop_sum_gdiv: $pop_sum_gdiv\npop_sum_poly_sites: $pop_sum_poly_sites\n".
#		"pop_sum_N_gdiv: $pop_sum_N_gdiv\tpop_sum_Nsites_gdiv: $pop_sum_Nsites_gdiv\n".
#		"pop_sum_S_gdiv: $pop_sum_S_gdiv\tpop_sum_Ssites_gdiv: $pop_sum_Ssites_gdiv\n";
	
	$pop_summary_line .= #"$pop_sum_ndiffs\t$pop_sum_sdiffs\t$pop_sum_ndiffs_vRef\t".
		#"$pop_sum_sdiffs_vRef\t" .
		"$pop_sum_nsites\t$pop_sum_ssites\t$pop_piN\t$pop_piS\t".
		"$pop_dNvREF\t$pop_dSvREF\t$pop_gdiv_poly\t$pop_N_gdiv\t$pop_S_gdiv\t".
		"$mean_gdiv\t".
		"$pop_sum_poly_sites\t".
		"$mean_gdiv_coding_poly\t".
		"$num_poly_coding_sites\t".
		"$mean_gdiv_nc_poly\t".
		"$num_poly_nc_sites\n";
	#}
	
	#$progress_period_count = 0;
	
	####### PRINT SNP REPORT TOTALS HERE BEFORE FINISHING WITH CURRENT FILE #########
	
	chdir('SNPGenie_Results');
	open(POP_SUMMARY,">>population\_summary\.txt");
	print POP_SUMMARY "$pop_summary_line";
	close POP_SUMMARY;
	chdir('..');
	
	print "$file_nm COMPLETED.\n";
	
	# Remove temp files
	#&remove_tempfiles;
	unlink $temp_snp_report_name;
	
	
	
	####### Generate random fastas: create all possible FASTA files	
	# This is for later use in between-population mutation rate estimation
	if($generate_random_fastas) {
		print "\nForming randomly generated representative sequences... ";
		my @all_possible_seqs;
		
		# Store data in variables to save room on screen and verify
		my $max_cov = 0;
		
		# Get maximum coverage, which will be the number of sequences we make
		foreach (sort {$a <=> $b} keys %hh_nc_position_info) {
			my $coverage = $hh_nc_position_info{$_}->{cov};
			if($coverage > $max_cov) {
				$max_cov = $coverage;
			}
		}
		
		foreach (sort {$a <=> $b} keys %hh_nc_position_info) {
			if($hh_nc_position_info{$_}->{polymorphic}) {
			
				my $this_A_prop = $hh_nc_position_info{$_}->{A_prop};
				my $this_C_prop = $hh_nc_position_info{$_}->{C_prop};
				my $this_G_prop = $hh_nc_position_info{$_}->{G_prop};
				my $this_T_prop = $hh_nc_position_info{$_}->{T_prop};
				my $this_ref = $hh_nc_position_info{$_}->{reference};
				my $this_cov = $hh_nc_position_info{$_}->{cov};
				
				my $this_A = ($this_A_prop * $max_cov);
				my $this_C = ($this_C_prop * $max_cov);
				my $this_G = ($this_G_prop * $max_cov);
				my $this_T = ($this_T_prop * $max_cov);
				
				#print "\nsite=$_ A=$this_A C=$this_C G=$this_G T=$this_T";
				
				# Round to nearest int
				$this_A = sprintf("%.0f", $this_A);
				$this_C = sprintf("%.0f", $this_C);
				$this_G = sprintf("%.0f", $this_G);
				$this_T = sprintf("%.0f", $this_T);
				
				#print "\nsite=$_ A=$this_A C=$this_C G=$this_G T=$this_T";
				
				# Find maj_nt
				my $maj_nt; # a variant nucleotide may have fixed
				my $curr_maj_count = 0;
				if($this_A > $curr_maj_count) {
					$curr_maj_count = $this_A;
					$maj_nt = 'A';
				}
				if($this_C > $curr_maj_count) {
					$curr_maj_count = $this_C;
					$maj_nt = 'C';
				} 
				if($this_G > $curr_maj_count) {
					$curr_maj_count = $this_G;
					$maj_nt = 'G';
				} 
				if($this_T > $curr_maj_count) {
					$curr_maj_count = $this_T;
					$maj_nt = 'T';
				}
				
				my $leftover_number = $max_cov - ($this_A + $this_C + $this_G + $this_T);
				
				if($leftover_number != 0) { # it's negative, our count is too large, take away from maj_nt
					# or it's positive, our count is too small, add to maj_nt
					if($maj_nt eq 'A') {
						$this_A += $leftover_number;
					} elsif($maj_nt eq 'C') {
						$this_C += $leftover_number;
					} elsif($maj_nt eq 'G') {
						$this_G += $leftover_number;
					} elsif($maj_nt eq 'T') {
						$this_T += $leftover_number;
					}	
				}
				
				# Now we build an array of all the nucleotides to add at the current site
				my @new_nts_arr;
				for(my $i = 1; $i <= $this_A; $i++) {
					push(@new_nts_arr,'A');
				}
				for(my $i = 1; $i <= $this_C; $i++) {
					push(@new_nts_arr,'C');
				}
				for(my $i = 1; $i <= $this_G; $i++) {
					push(@new_nts_arr,'G');
				}
				for(my $i = 1; $i <= $this_T; $i++) {
					push(@new_nts_arr,'T');
				}
				
				for(my $i = 1; $i <= $max_cov; $i++) {
					my $j = $i - 1;
					
					srand();
					my $rand_nt_index = int(rand(scalar(@new_nts_arr)));
					my $rand_nt = $new_nts_arr[$rand_nt_index];
					
					$all_possible_seqs[$j] .= $rand_nt;
					
					# Delete element from array
					splice(@new_nts_arr, $rand_nt_index, 1);
					
				}
				#print "\nsite=$_ A=$this_A C=$this_C G=$this_G T=$this_T ref=$this_ref cov=$this_cov";
				
			} else { # not polymorphic, so just add the reference nucleotide to all
				my $this_nt = $seq_by_index_arr[$_ - 1];
				
				for(my $i = 1; $i <= $max_cov; $i++) {
					my $j = $i - 1;
					$all_possible_seqs[$j] .= $this_nt;
				}
			}
		}
		
		#print "\nMaximum coverage is $max_cov\n\n";
		
		my $counter = 1;
		open(OUT_FASTAS, ">>rand_seqs.fasta");
		foreach(@all_possible_seqs) {
			my $seq_length = length($_);
			
			print OUT_FASTAS ">seq_num_$counter\n";
			
			for(my $i=0; $i<$seq_length; $i+=60) {
				if($i>($seq_length - 60)) {
					my $line = substr($_, $i);
					print OUT_FASTAS "$line\n";
					last;
				} else {
					my $line = substr($_, $i, 60);
					print OUT_FASTAS "$line\n";
				}
			}
			
			$counter++;
		}
		close OUT_FASTAS;
		print "COMPLETED.\n";
	} #######	End generate random fastas
	
} # Finished with SNP Report

# Perform a SLIDING WINDOW if asked
if($slidingwindow) {
	print "\n\nPerforming sliding window for all files, length $slidingwindow codons...\n\n";
	&sliding_window($slidingwindow);
}


# GROUP and METAPOPULATION ANALYSES
### my %site_sample_props;
### 
### # MANUALLY CHANGE THIS
### my @sorted_groups = qw/BIRD MOSQ/;
### print "\nSorted groups are: @sorted_groups\n";
### 
### my @sorted_samples; # = sort (keys %{$master_frequencies_hh{1}}); # grab first site's
### #print "\nSorted samples are: @sorted_samples\n";
### 
### foreach my $sample_name (sort (keys %{$master_frequencies_hh{1}})) {
### 	my $group_name_match = 0;
### 	
### 	foreach my $group_name (@sorted_groups) {
### 		if($group_name eq $sample_name) {
### 			$group_name_match = 1;
### 		}
### 	}
### 	
### 	unless($group_name_match == 1) {
### 		push(@sorted_samples, $sample_name);
### 	}
### }
### 
### print "\nCleansed sorted samples are: @sorted_samples\n";
### 
### 
### 
### # METAPOPULATION ANALYSIS # push(@{$master_frequencies_hh{$position}-> {BIRD} -> {A_props_arr}} , $A_prop);
### open(OUT_METAPOP, ">>snpgenie_metapop_results.txt");
### my $metapop_header = "site\t";
### my @nucleotides = qw(A C G T);
### 
### foreach my $sample (@sorted_groups) {
### 	$metapop_header .= "$sample\_n\t";
### }
### 
### foreach my $nucleotide (@nucleotides) {
### 	foreach my $sample (@sorted_groups) {
### 		$metapop_header .= "$sample\_$nucleotide\t";
### 	}
### }
### 
### $metapop_header .= "curr_max_diff_nt\tcurr_max_diff";
### 
### print OUT_METAPOP "$metapop_header\n";
### 		
### for (my $i = 0; $i < length($seq); $i++) {
### 	my $position = $i+1;
### 	
### 	my $output_line = "$position\t";
### 	
### 	foreach my $sample (@sorted_groups) {
### 	
### 		my $sample_size = scalar(@{$master_frequencies_hh{$position}->{$sample}->{A_props_arr}}); # The A's will do for sample size
### 		$output_line .= "$sample_size\t";
### 		
### 	}
### 	
### 	# MUST MAINTAIN ORDER OF NUCLEOTIDES HERE
### 	#my $interesting_site = '';
### 	my $curr_max_diff = 0;
### 	my $curr_max_diff_nt;
### 	
### 	my $last_A_prop = '';
	
##	# QUICK FIX
##	my @A_freqs_BIRD = @{$master_frequencies_hh{$position}->{BIRD}->{A_props_arr}};
##	my @A_freqs_MOSQ = @{$master_frequencies_hh{$position}->{MOSQ}->{A_props_arr}};
##	
##	my $A_BIRD_numerator = sum(@A_freqs_BIRD);
##	my $A_BIRD_denominator = scalar(@A_freqs_BIRD);
##	my $A_BIRD_mean_prop = $A_BIRD_numerator / $A_BIRD_denominator;
##	
##	
##	my $A_MOSQ_numerator = sum(@A_freqs_MOSQ);
##	my $A_MOSQ_denominator = scalar(@A_freqs_MOSQ);
##	my $A_MOSQ_mean_prop = $A_MOSQ_numerator / $A_MOSQ_denominator;
##	
##	
##	
##	my @C_freqs_BIRD = @{$master_frequencies_hh{$position}->{BIRD}->{C_props_arr}};
##	my @C_freqs_MOSQ = @{$master_frequencies_hh{$position}->{MOSQ}->{C_props_arr}};
##	
##	my $C_BIRD_numerator = sum(@C_freqs_BIRD);
##	my $C_BIRD_denominator = scalar(@C_freqs_BIRD);
##	my $C_BIRD_mean_prop = $C_BIRD_numerator / $C_BIRD_denominator;
##	
##	
##	my $C_MOSQ_numerator = sum(@C_freqs_MOSQ);
##	my $C_MOSQ_denominator = scalar(@C_freqs_MOSQ);
##	my $C_MOSQ_mean_prop = $C_MOSQ_numerator / $C_MOSQ_denominator;
##	
##	
##	
##	my @G_freqs_BIRD = @{$master_frequencies_hh{$position}->{BIRD}->{G_props_arr}};
##	my @G_freqs_MOSQ = @{$master_frequencies_hh{$position}->{MOSQ}->{G_props_arr}};
##	
##	my @G_freqs_BIRD = @{$master_frequencies_hh{$position}->{BIRD}->{G_props_arr}};
##	my @G_freqs_MOSQ = @{$master_frequencies_hh{$position}->{MOSQ}->{G_props_arr}};
##	
##	my $G_BIRD_numerator = sum(@G_freqs_BIRD);
##	my $G_BIRD_denominator = scalar(@G_freqs_BIRD);
##	my $G_BIRD_mean_prop = $G_BIRD_numerator / $G_BIRD_denominator;
##	
##	
##	my $G_MOSQ_numerator = sum(@G_freqs_MOSQ);
##	my $G_MOSQ_denominator = scalar(@G_freqs_MOSQ);
##	my $G_MOSQ_mean_prop = $G_MOSQ_numerator / $G_MOSQ_denominator;
##	
##	
##	
##	
##	my @T_freqs_BIRD = @{$master_frequencies_hh{$position}->{BIRD}->{T_props_arr}};
##	my @T_freqs_MOSQ = @{$master_frequencies_hh{$position}->{MOSQ}->{T_props_arr}};
##	
##	
##	my @T_freqs_BIRD = @{$master_frequencies_hh{$position}->{BIRD}->{T_props_arr}};
##	my @T_freqs_MOSQ = @{$master_frequencies_hh{$position}->{MOSQ}->{T_props_arr}};
##	
##	my $T_BIRD_numerator = sum(@T_freqs_BIRD);
##	my $T_BIRD_denominator = scalar(@T_freqs_BIRD);
##	my $T_BIRD_mean_prop = $T_BIRD_numerator / $T_BIRD_denominator;
##	
##	
##	my $T_MOSQ_numerator = sum(@T_freqs_MOSQ);
##	my $T_MOSQ_denominator = scalar(@T_freqs_MOSQ);
##	my $T_MOSQ_mean_prop = $T_MOSQ_numerator / $T_MOSQ_denominator;
##	
##	
##	if($position == 4543) {
##		print "\nFreq of BIRD A at $position is @A_freqs_BIRD, mean $A_BIRD_mean_prop\n";
##		print "\nFreq of MOSQ A at $position is @A_freqs_MOSQ, mean $A_MOSQ_mean_prop\n";
##		print "\nFreq of BIRD C at $position is @C_freqs_BIRD, mean $C_BIRD_mean_prop\n";
##		print "\nFreq of MOSQ C at $position is @C_freqs_MOSQ, mean $C_MOSQ_mean_prop\n";
##		print "\nFreq of BIRD G at $position is @G_freqs_BIRD, mean $G_BIRD_mean_prop\n";
##		print "\nFreq of MOSQ G at $position is @G_freqs_MOSQ, mean $G_MOSQ_mean_prop\n";
##		print "\nFreq of BIRD T at $position is @T_freqs_BIRD, mean $T_BIRD_mean_prop\n";
##		print "\nFreq of MOSQ T at $position is @T_freqs_MOSQ, mean $T_MOSQ_mean_prop\n";
##	}
	
	
	
### 	# THIS DOES NOT WORK RIGHT, frequencies being attributed to incorrect nucleotides
### 	foreach my $sample (@sorted_groups) {
### 		
### 		my @A_freqs = @{$master_frequencies_hh{$position}->{$sample}->{A_props_arr}};
### 		my $A_numerator = sum(@A_freqs);
### 		my $A_denominator = scalar(@A_freqs);
### 		my $A_mean_prop = $A_numerator / $A_denominator;
### 		#print "\nFreq of $sample A at $position is @A_freqs, mean $A_mean_prop\n";
### 		
### 		if($position == 4543) {
### 			print "\nsite=$position sample=$sample A A_num=$A_numerator A_denom=$A_denominator A_prop=$A_mean_prop\nA_freqs: @A_freqs\n";
### 		}
### 		
### 		if($last_A_prop eq '') {
### 			$last_A_prop = $A_mean_prop;
### 		} else {
### #			if(abs($A_mean_prop - $last_A_prop) >= $prop_diff_threshold) {
### #				$interesting_site = 1;
### #			}
### 			if(abs($A_mean_prop - $last_A_prop) > $curr_max_diff) {
### 				$curr_max_diff = abs($A_mean_prop - $last_A_prop);
### 				$curr_max_diff_nt = 'A';
### 			}
### 		}
### 		
### 		$output_line .= "$A_mean_prop\t";
### 	}
### 	
### 	my $last_C_prop = '';
### 	
### 	foreach my $sample (@sorted_groups) {
### 		
### 		my @C_freqs = @{$master_frequencies_hh{$position}->{$sample}->{C_props_arr}};
### 		my $C_numerator = sum(@C_freqs);
### 		my $C_denominator = scalar(@C_freqs);
### 		my $C_mean_prop = $C_numerator / $C_denominator;
### 		#print "\nFreq of $sample C at $position is @C_freqs, mean $C_mean_prop\n";
### 		
### 		if($position == 4543) {
### 			print "\nsite=$position sample=$sample C C_num=$C_numerator C_denom=$C_denominator C_prop=$C_mean_prop\nC_freqs: @C_freqs\n";
### 		}
### 		
### 		if($last_C_prop eq '') {
### 			$last_C_prop = $C_mean_prop;
### 		} else {
### 			if(abs($C_mean_prop - $last_C_prop) > $curr_max_diff) {
### 				$curr_max_diff = abs($C_mean_prop - $last_C_prop);
### 				$curr_max_diff_nt = 'C';
### 			}
### 		}
### 		
### 		$output_line .= "$C_mean_prop\t";
### 	}
### 	
### 	my $last_G_prop = '';
### 	
### 	foreach my $sample (@sorted_groups) {
### 		
### 		my @G_freqs = @{$master_frequencies_hh{$position}->{$sample}->{G_props_arr}};
### 		my $G_numerator = sum(@G_freqs);
### 		my $G_denominator = scalar(@G_freqs);
### 		my $G_mean_prop = $G_numerator / $G_denominator;
### 		#print "\nFreq of $sample G at $position is @G_freqs, mean $G_mean_prop\n";
### 		
### 		if($position == 4543) {
### 			print "\nsite=$position sample=$sample G G_num=$G_numerator G_denom=$G_denominator G_prop=$G_mean_prop\nG_freqs: @G_freqs\n";
### 		}
### 		
### 		if($last_G_prop eq '') {
### 			$last_G_prop = $G_mean_prop;
### 		} else {
### 			if(abs($G_mean_prop - $last_G_prop) > $curr_max_diff) {
### 				$curr_max_diff = abs($G_mean_prop - $last_G_prop);
### 				$curr_max_diff_nt = 'G';
### 			}
### 		}
### 		
### 		
### 		$output_line .= "$G_mean_prop\t";
### 	}
### 	
### 	my $last_T_prop = '';
### 	
### 	foreach my $sample (@sorted_groups) {
### 	
### 		my @T_freqs = @{$master_frequencies_hh{$position}->{$sample}->{T_props_arr}};
### 		my $T_numerator = sum(@T_freqs);
### 		my $T_denominator = scalar(@T_freqs);
### 		my $T_mean_prop = $T_numerator / $T_denominator;
### 		#print "\nFreq of $sample T at $position is @T_freqs, mean $T_mean_prop\n";
### 		
### 		if($position == 4543) {
### 			print "\nsite=$position sample=$sample T T_num=$T_numerator T_denom=$T_denominator T_prop=$T_mean_prop\nT_freqs: @T_freqs\n";
### 		}
### 		
### 		if($last_T_prop eq '') {
### 			$last_T_prop = $T_mean_prop;
### 		} else {
### 			if(abs($T_mean_prop - $last_T_prop) > $curr_max_diff) {
### 				$curr_max_diff = abs($T_mean_prop - $last_T_prop);
### 				$curr_max_diff_nt = 'T';
### 			}
### 		}
### 		
### 		
### 		$output_line .= "$T_mean_prop\t";
### 	}
### 	
### 	$output_line .= "$curr_max_diff_nt\t$curr_max_diff";
### 	
### 	print OUT_METAPOP "$output_line\n";
### 	
### }
### 
### close OUT_METAPOP;



### 
### # FST ANALYSIS # push(@{$master_frequencies_hh{$position}-> {BIRD} -> {A_props_arr}} , $A_prop);
### #				^ this includes all, not just polymorphic sites
### 
### open(OUT_FST, ">>snpgenie_FST_results.txt");
### 
### # First column is population ID
### my $fst_header = "population\t";
### 
### # Subsequent columns are sites. Could do all, or just polymorphic. Start with all.
### # While doing this, determine meta-consensus nucleotides at each site
### 
### my %meta_consensus_nts; # $meta_consensus_nts{site}->{nt/freq} = 'A'/0.75
### 
### for (my $i = 0; $i < length($seq); $i++) {
### 	my $site = $i+1;
### 
### 	my @this_site_A_freqs;
### 	my @this_site_C_freqs;
### 	my @this_site_G_freqs;
### 	my @this_site_T_freqs;
### 	
### 	foreach my $sample_name (@sorted_samples) { # GROUP names (e.g., BIRD and MOSQ) have been removed
### 		push(@this_site_A_freqs, $master_frequencies_hh{$site}->{$sample_name}->{A_freq});
### 		push(@this_site_C_freqs, $master_frequencies_hh{$site}->{$sample_name}->{C_freq});
### 		push(@this_site_G_freqs, $master_frequencies_hh{$site}->{$sample_name}->{G_freq});
### 		push(@this_site_T_freqs, $master_frequencies_hh{$site}->{$sample_name}->{T_freq});
### 	}
### 	
### 	# Find meta-consensus nt and its mean freq
### 	my $meta_A_freq = sum(@this_site_A_freqs) / scalar(@this_site_A_freqs);
### 	my $meta_C_freq = sum(@this_site_C_freqs) / scalar(@this_site_C_freqs);
### 	my $meta_G_freq = sum(@this_site_G_freqs) / scalar(@this_site_G_freqs);
### 	my $meta_T_freq = sum(@this_site_T_freqs) / scalar(@this_site_T_freqs);
### 	
### 	my $meta_consensus_nt;
### 	my $meta_consensus_nt_freq = 0;
### 	
### 	if($meta_A_freq > $meta_consensus_nt_freq) {
### 		$meta_consensus_nt = 'A';
### 		$meta_consensus_nt_freq = $meta_A_freq;
### 	}
### 	
### 	if($meta_C_freq > $meta_consensus_nt_freq) {
### 		$meta_consensus_nt = 'C';
### 		$meta_consensus_nt_freq = $meta_C_freq;
### 	}
### 	
### 	if($meta_G_freq > $meta_consensus_nt_freq) {
### 		$meta_consensus_nt = 'G';
### 		$meta_consensus_nt_freq = $meta_G_freq;
### 	}
### 	
### 	if($meta_T_freq > $meta_consensus_nt_freq) {
### 		$meta_consensus_nt = 'T';
### 		$meta_consensus_nt_freq = $meta_T_freq;
### 	}
### 	
### 	$meta_consensus_nts{$site}->{nt} = $meta_consensus_nt;
### 	$meta_consensus_nts{$site}->{freq} = $meta_consensus_nt_freq;
### 	
### 	$fst_header .= "site_$site\_$meta_consensus_nt\t";
### 	
### }
### 
### chop($fst_header);
### 
### #my @nucleotides = qw(A C G T); # done above
### 
### print OUT_FST "$fst_header\n";
### 
### 
### 
### 
### 
### 
### foreach my $sample_name (@sorted_samples) { # each LINE of output is a sample
### 	my $this_sample_out_line = "$sample_name\t";
### 	
### 	for (my $i = 0; $i < length($seq); $i++) {
### 		my $site = $i+1;
### 		my $site_meta_consensus_nt = $meta_consensus_nts{$site}->{nt};
### 		
### 		my $freq_key = "$site_meta_consensus_nt\_freq";
### 		my $meta_consensus_freq_here = $master_frequencies_hh{$site}->{$sample_name}->{$freq_key};
### 		
### 		$this_sample_out_line .= "$meta_consensus_freq_here\t";
### 	}
### 	
### 	chop($this_sample_out_line);
### 	print OUT_FST "$this_sample_out_line\n";
### 	
### }



#for (my $i = 0; $i < length($seq); $i++) {
#	my $position = $i+1;
#	
#	my $output_line = "$position\t";
#	
#	foreach my $sample (@sorted_groups) {
#	
#		my $sample_size = scalar(@{$master_frequencies_hh{$position}->{$sample}->{A_props_arr}}); # The A's will do for sample size
#		$output_line .= "$sample_size\t";
#		
#	}
#	
#	# MUST MAINTAIN ORDER OF NUCLEOTIDES HERE
#	#my $interesting_site = '';
#	my $curr_max_diff = 0;
#	my $curr_max_diff_nt;
#	
#	my $last_A_prop = '';	
#	
#	# THIS DOES NOT WORK RIGHT, frequencies being attributed to incorrect nucleotides
#	foreach my $sample (@sorted_groups) {
#		
#		my @A_freqs = @{$master_frequencies_hh{$position}->{$sample}->{A_props_arr}};
#		my $A_numerator = sum(@A_freqs);
#		my $A_denominator = scalar(@A_freqs);
#		my $A_mean_prop = $A_numerator / $A_denominator;
#		#print "\nFreq of $sample A at $position is @A_freqs, mean $A_mean_prop\n";
#		
#		if($position == 4543) {
#			print "\nsite=$position sample=$sample A A_num=$A_numerator A_denom=$A_denominator A_prop=$A_mean_prop\nA_freqs: @A_freqs\n";
#		}
#		
#		if($last_A_prop eq '') {
#			$last_A_prop = $A_mean_prop;
#		} else {
##			if(abs($A_mean_prop - $last_A_prop) >= $prop_diff_threshold) {
##				$interesting_site = 1;
##			}
#			if(abs($A_mean_prop - $last_A_prop) > $curr_max_diff) {
#				$curr_max_diff = abs($A_mean_prop - $last_A_prop);
#				$curr_max_diff_nt = 'A';
#			}
#		}
#		
#		$output_line .= "$A_mean_prop\t";
#	}
#	
#	my $last_C_prop = '';
#	
#	foreach my $sample (@sorted_groups) {
#		
#		my @C_freqs = @{$master_frequencies_hh{$position}->{$sample}->{C_props_arr}};
#		my $C_numerator = sum(@C_freqs);
#		my $C_denominator = scalar(@C_freqs);
#		my $C_mean_prop = $C_numerator / $C_denominator;
#		#print "\nFreq of $sample C at $position is @C_freqs, mean $C_mean_prop\n";
#		
#		if($position == 4543) {
#			print "\nsite=$position sample=$sample C C_num=$C_numerator C_denom=$C_denominator C_prop=$C_mean_prop\nC_freqs: @C_freqs\n";
#		}
#		
#		if($last_C_prop eq '') {
#			$last_C_prop = $C_mean_prop;
#		} else {
#			if(abs($C_mean_prop - $last_C_prop) > $curr_max_diff) {
#				$curr_max_diff = abs($C_mean_prop - $last_C_prop);
#				$curr_max_diff_nt = 'C';
#			}
#		}
#		
#		$output_line .= "$C_mean_prop\t";
#	}
#	
#	my $last_G_prop = '';
#	
#	foreach my $sample (@sorted_groups) {
#		
#		my @G_freqs = @{$master_frequencies_hh{$position}->{$sample}->{G_props_arr}};
#		my $G_numerator = sum(@G_freqs);
#		my $G_denominator = scalar(@G_freqs);
#		my $G_mean_prop = $G_numerator / $G_denominator;
#		#print "\nFreq of $sample G at $position is @G_freqs, mean $G_mean_prop\n";
#		
#		if($position == 4543) {
#			print "\nsite=$position sample=$sample G G_num=$G_numerator G_denom=$G_denominator G_prop=$G_mean_prop\nG_freqs: @G_freqs\n";
#		}
#		
#		if($last_G_prop eq '') {
#			$last_G_prop = $G_mean_prop;
#		} else {
#			if(abs($G_mean_prop - $last_G_prop) > $curr_max_diff) {
#				$curr_max_diff = abs($G_mean_prop - $last_G_prop);
#				$curr_max_diff_nt = 'G';
#			}
#		}
#		
#		
#		$output_line .= "$G_mean_prop\t";
#	}
#	
#	my $last_T_prop = '';
#	
#	foreach my $sample (@sorted_groups) {
#	
#		my @T_freqs = @{$master_frequencies_hh{$position}->{$sample}->{T_props_arr}};
#		my $T_numerator = sum(@T_freqs);
#		my $T_denominator = scalar(@T_freqs);
#		my $T_mean_prop = $T_numerator / $T_denominator;
#		#print "\nFreq of $sample T at $position is @T_freqs, mean $T_mean_prop\n";
#		
#		if($position == 4543) {
#			print "\nsite=$position sample=$sample T T_num=$T_numerator T_denom=$T_denominator T_prop=$T_mean_prop\nT_freqs: @T_freqs\n";
#		}
#		
#		if($last_T_prop eq '') {
#			$last_T_prop = $T_mean_prop;
#		} else {
#			if(abs($T_mean_prop - $last_T_prop) > $curr_max_diff) {
#				$curr_max_diff = abs($T_mean_prop - $last_T_prop);
#				$curr_max_diff_nt = 'T';
#			}
#		}
#		
#		
#		$output_line .= "$T_mean_prop\t";
#	}
#	
#	$output_line .= "$curr_max_diff_nt\t$curr_max_diff";
#	
#	print OUT_FST "$output_line\n";
#	
#}



### close OUT_FST;



# Print a completion message to screen
&end_the_program;

exit;

#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################


#########################################################################################
sub median {
    my @values = sort {$a <=> $b} @_;
    my $length = @values;
    if($length%2) { # odd number of elements: return middle
        return $values[int($length/2)];
    } else { # even number of elements: return mean of middle two
        return ($values[int($length/2)-1] + $values[int($length/2)])/2;
    }
}


#########################################################################################
sub get_header_names {
	# Originally assumed that we've received a tempfile ending in "_snpg9temp.txt"
	# However, we're now calling it at least once before creating the tempfile to 
	# see what kind of processing (e.g., Geneious to CLC) is needed prior to tempfile
	# creation. Must include capability to get headers for .CSV file
	my ($curr_snp_report_filename,$filename) = @_;
	#print "\n$curr_snp_report_filename\n";
	
	#my $newline_char = &detect_newline_char($curr_snp_report_filename);
	#my $old_newline = $/;
	#$/ = $newline_char;
	
	my $seen_tab_delimited = 0;
	my $seen_comma_delimited = 0;
	my $seen_vcf_tab_delimited = 0;
	my @line_arr;
	
	my $line = 0;
	open (CURRINFILE, $curr_snp_report_filename);
	#seek(CURRINFILE,0,0);
	while (<CURRINFILE>) {
		#print "$_";
		if($line == 0) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			if($_ =~/\t\w+\t/) { # it's TAB-delimited
				@line_arr = split(/\t/,$_,-1);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,\w+,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_,-1);
				#print "COMMA!!!!!";
				last;
			}

			$line++;
		} elsif($line > 0 && $_ =~ /^##/) {
			$line++;
		} elsif($line > 0 && ($_ =~ /^#CHROM/)) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			if($_ =~/\t/) { # it's TAB-delimited
				@line_arr = split(/\t/,$_,-1);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_,-1);
				#print "COMMA!!!!!";
				last;
			} else {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | CODON | WARNING
				
				# No change OR error should occur if the file does not, in fact, end
				# with this SUFFIX
				my $file_nm = $curr_snp_report_filename;
				#$file_nm =~ s/_snpg9temp.txt/.txt/;
				$file_nm =~ s/_\w\w\w\w.txt/.txt/;
				
				print ERROR_FILE "$filename\tN/A\tN/A\t".
					"File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
				close ERROR_FILE;
				chdir('..');
				
				#unlink $curr_snp_report_filename;
				
				die "\n\n## WARNING: The SNP Report $filename is ".
					"not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie ".
					"terminated\n\n";
			}
		} else {
			chdir('SNPGenie_Results');
			open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
			# FILE | PRODUCT | SITE | CODON | WARNING
			
			# No change OR error should occur if the file does not, in fact, end
			# with this SUFFIX
			my $file_nm = $curr_snp_report_filename;
			#$file_nm =~ s/_snpg9temp.txt/.txt/;
			$file_nm =~ s/_\w\w\w\w.txt/.txt/;
			
			print ERROR_FILE "$filename\tN/A\tN/A\t".
				"File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
			close ERROR_FILE;
			chdir('..');
			
			#unlink $curr_snp_report_filename;
			
			die "\n\n## WARNING: The SNP Report $filename is ".
				"not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie ".
				"terminated\n\n";
		}
	}
	seek(CURRINFILE,0,0);
	close CURRINFILE;
	#$/ = $old_newline;
	return @line_arr;
}

#########################################################################################
#sub get_product_names {
#	my ($curr_snp_report_filename,$index_over_annot) = @_;
#	#print "\n\n$curr_snp_report_filename\n\n";
#	my $line = 0;
#	my %products_hash;
#	open (CURRINFILE, $curr_snp_report_filename);
#	while (<CURRINFILE>) {
#		if ($line == 0) {
#			$line++;
#		} else {
#			#chomp;
#			
#			# CHOMP for 3 operating systems
#			if($_ =~ /\r\n$/) {
#				$_ =~ s/\r\n//;
#			} elsif($_ =~ /\r$/) {
#				$_ =~ s/\r//;
#			} elsif($_ =~ /\n$/) {
#				$_ =~ s/\n//;
#			}
#			
#			my @line_arr = split(/\t/,$_);
#			my $over_annot = $line_arr[$index_over_annot];
#			
#			if ($over_annot =~/Mature peptide: ([\w\s\.']+)/) {
#				if (! exists $products_hash{$1}) {
#					$products_hash{$1} = 1;
#				}
#			} #elsif ($over_annot =~/CDS: (\w+)/) { # ORIGINAL
#			#	if (! exists $products_hash{$1}) {
#			#		$products_hash{$1} = 1;
#			#	}
#			#}
#			
#			if ($over_annot =~/CDS: ([\w\s\.']+)/) {
#				if (! exists $products_hash{$1}) {
#					$products_hash{$1} = 1;
#				}
#			} #elsif ($over_annot =~/Gene: (\w+)/) { # CHANGED THIS TO IGNORE GENE-ONLY ANNOTATIONS: WANT CDS
#			#	if (! exists $products_hash{$1}) {
#			#		$products_hash{$1} = 1;
#			#	}
#			#} 
#		}
#	}
#	close CURRINFILE;
#	my @product_names = keys %products_hash;
#	#print "\n@product_names\n\n";
#	return @product_names;
#}

#########################################################################################
sub get_product_names_from_gtf {
	my ($cds_file) = @_;
	#print "\n\n$cds_file\n\n";
	my %products_hash;
	open (CURRINFILE, $cds_file);
	while (<CURRINFILE>) {
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			#print "this_line: $_";
			if($_ =~/\s*gene_id\s*\"gene\:([\w\s\.\-\:']+)\"/) { # transcript_id not a problem
				$products_hash{$1} = 1;
				$seen_sense_strand_products = 1;
			} elsif($_ =~ /\s*gene_id\s*\"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
				$seen_sense_strand_products = 1;
			} elsif($_ =~/\s*gene_id\s*\"([\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
				$seen_sense_strand_products = 1;
			} else {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | CODON | WARNING
				
				print ERROR_FILE "$cds_file\tN/A\tN/A\t".
					"CDS annotation(s) does not have a gene_id. SNPGenie terminated.\n";
				close ERROR_FILE;
				chdir('..');
				
				#unlink $curr_snp_report_filename;
				
				die "\n\n## WARNING: CDS annotation(s) in $cds_file does not have a ".
					"gene_id. SNPGenie terminated.\n\n";
			}
		}
	}
	close CURRINFILE;
	
	my @product_names = keys %products_hash;
	#print "\n@product_names\n\n";
	return @product_names;
}

#########################################################################################
#sub get_product_names_vcf {
#	my ($gtf_file_nm) = @_;
#	my %products_hash;
#	open (CURRINFILE, $gtf_file_nm);
#	while (<CURRINFILE>) {
#		chomp;
#		# CHOMP for 3 operating systems
#		#if($_ =~ /\r\n$/) {
#		#	$_ =~ s/\r\n//;
#		#} elsif($_ =~ /\r$/) {
#		#	$_ =~ s/\r//;
#		#} elsif($_ =~ /\n$/) {
#		#	$_ =~ s/\n//;
#		#}
#		
#		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
#			if($_ =~/gene_id \"gene\:([\w\s\.']+)\"/) {
#				my $product = $1;
#				
#				if ((! exists $products_hash{$product}) && ($product ne '')) {
#					$products_hash{$product} = 1;
#				}
#			} elsif($_ =~ /gene_id \"([\w\s\.']+ [\w\s\.']+)\"/) {
#				my $product = $1;
#				
#				if ((! exists $products_hash{$product}) && ($product ne '')) {
#					$products_hash{$product} = 1;
#				}
#			} elsif($_ =~ /gene_id \"([\w\s\.']+)\"/) {
#				my $product = $1;
#				
#				if ((! exists $products_hash{$product}) && ($product ne '')) {
#					$products_hash{$product} = 1;
#				}
#			}
#		}
#		
#	}
#	close CURRINFILE;
#	my @product_names = keys %products_hash;
#	#foreach (@product_names) {
#	#	print "$_\n";
#	#}
#	return @product_names;
#}

#########################################################################################
sub determine_complement_mode {
	my ($gtf_file_nm) = @_;
	my $complement_mode;
	open (CURRINFILE, $gtf_file_nm);
	while (<CURRINFILE>) {
		chomp;
		# CHOMP for 3 operating systems
		if($_ =~ /\r\n$/) {
			$_ =~ s/\r\n//;
		} elsif($_ =~ /\r$/) {
			$_ =~ s/\r//;
		} elsif($_ =~ /\n$/) {
			$_ =~ s/\n//;
		}
		
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t-/) {
			$complement_mode = 1;
			last;
		}
	}
	close CURRINFILE;
	return $complement_mode;
}

#########################################################################################
sub reverse_complement_from_fasta {
	my ($filename) = @_;
	
	# Read in the sequence from the file
	my $seq = '';
	
	open(IN, "$filename") or die "\nCould not open FASTA file $filename\n\n";
	while(<IN>) {
		unless(/>/) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
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
# Obtains all file names in current directory ending in .fa and/or .fasta
sub get_fasta_file_names { 
	my @fasta_file_names = glob "*.fa";
	my @other_fasta_file_names;
	
	#if (scalar(@fasta_file_names) == 0) {
	#	@fasta_file_names = glob "*.fasta";
	#} else {
		@other_fasta_file_names = glob "*.fasta";
		push (@fasta_file_names,@other_fasta_file_names);
	#}
	
	#print "\n\n@fasta_file_names\n\n";
	
	if (scalar(@fasta_file_names) == 0) {
		#chdir('SNPGenie_Results');
		#open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		## FILE | PRODUCT | SITE | CODON | WARNING
		#print ERROR_FILE "N/A\tN/A\tN/A\t".
		#	"No FASTA (.fa or .fasta) files in directory. SNPGenie terminated.\n";
		#close ERROR_FILE;
		#chdir('..');
		
		die "\n\n## WARNING: There are no .fa or .fasta files. SNPGenie terminated.\n\n";
	}
	
	#print "\n\n@fasta_file_names\n\n";
	return 	@fasta_file_names;
}

#########################################################################################
# Obtains all file names in current directory ending in .txt
sub get_txt_file_names { 
	my @txt_file_names = glob "*.txt";
#	if (scalar (@txt_file_names) == 0) {
		#chdir('SNPGenie_Results');
		#open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		## FILE | PRODUCT | SITE | CODON | WARNING
		#print ERROR_FILE "N/A\tN/A\tN/A\t".
		##	"No SNP Reports (.txt) files in directory. SNPGenie terminated.\n";
		#	"No SNP Reports (.txt) files in directory.\n";
		#close ERROR_FILE;
		#chdir('..');
		
		#die "\n\n## WARNING: There are no .txt SNP Reports. SNPGenie terminated.\n\n";
#		print "\n\n## WARNING: There are no .txt SNP Reports.\n\n";
#	}
	return 	@txt_file_names;
}

#########################################################################################
# Obtains all file names in current directory ending in .txt
sub get_csv_file_names { 
	my @csv_file_names = glob "*.csv";
#	if (scalar (@csv_file_names) == 0) {
		#chdir('SNPGenie_Results');
		#open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		## FILE | PRODUCT | SITE | CODON | WARNING
		#print ERROR_FILE "N/A\tN/A\tN/A\t".
		##	"No SNP Reports (.txt) files in directory. SNPGenie terminated.\n";
		#	"No SNP Reports (.csv) files in directory.\n";
		#close ERROR_FILE;
		#chdir('..');
		
		#die "\n\n## WARNING: There are no .csv SNP Reports. SNPGenie terminated.\n\n";
#		print "\n\n## WARNING: There are no .csv SNP Reports.\n\n";
#	}
	return 	@csv_file_names;
}

#########################################################################################
# Obtains all file names in current directory ending in .txt
sub get_vcf_file_names { 
	my @csv_file_names = glob "*.vcf";
#	if (scalar (@csv_file_names) == 0) {
		#chdir('SNPGenie_Results');
		#open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		## FILE | PRODUCT | SITE | CODON | WARNING
		#print ERROR_FILE "N/A\tN/A\tN/A\t".
		##	"No SNP Reports (.txt) files in directory. SNPGenie terminated.\n";
		#	"No SNP Reports (.csv) files in directory.\n";
		#close ERROR_FILE;
		#chdir('..');
		
		#die "\n\n## WARNING: There are no .csv SNP Reports. SNPGenie terminated.\n\n";
#		print "\n\n## WARNING: There are no .csv SNP Reports.\n\n";
#	}
	return 	@csv_file_names;
}

#########################################################################################
# (1) Is passed an array of all files to process, assumed to be in the working directory
# (2) Processes those SNP reports for various anomalies in the Overlapping annotations 
# column, chiefly to get variants that are present in multiple products ON THEIR OWN 
# LINES. An inelegant solution, but it beats re-doing the whole algorithm.
# (3) Obtains and returns all file names in current directory ending in _snpg9temp.txt
#sub create_tempfiles {
#	#my @snp_report_file_names_arr = @_;
#	#print "\nMy SNP Reports are: @snp_report_file_names\n\n";
#	
##	foreach (@snp_report_file_names_arr) {
#		my $snpr_file_name = $_[0];
#		
#		# Generate new file name prefix
#		my $new_file_prefix;
#		if($snpr_file_name =~/\.txt/) { 
#			$new_file_prefix = $`;
#		} else {
#			$new_file_prefix = "inFile";
#		}
#		
#		# New file name
#		#my $new_file_name = $new_file_prefix . "_snpg9temp.txt";
#		my $new_file_name = $new_file_prefix . "_\w\w\w\w.txt";
#		
#		my @header_names_arr = &get_header_names($snpr_file_name,$snpr_file_name);
#		#print "@header_names_arr";
#		
#		# Important lesson here: if a function is called, $_ is updated to its return value,
#		# even when you assign that value to something else.
#		#print "\n\n$_\n\n";
#		
#		my $line = 0;
#		open(SNPR_FILE, $snpr_file_name) or die "Could not open SNP Report file ".
#			"$snpr_file_name\n";
#		while (<SNPR_FILE>) {
#			if($line == 0) {
#				open(NEW_SNPR,">$new_file_name");
#				print NEW_SNPR $_;
#				close NEW_SNPR;
#				$line++;
#			} else {
#				open(NEW_SNPR,">>$new_file_name");
#				# First, replace the damned parentheses
#				if($_ =~ /CDS: ""(\w+ \w+)""/) {
#					my $replacement = $1;
#					#print "\n\nI had a replacement\n\n";
#					$_ =~ s/CDS: ""\w+ \w+""/CDS: $replacement/;
#				}
#				
#				# To account for prime (') symbols in the ORF names
#				if($_ =~ /ORF: ([\w\s\.']+), ORF: ([\w\s\.']+), ORF: ([\w\s\.']+)/) {
#					#print "We got here2\n\n";
#					my $ORF1 = $1;
#					my $ORF2 = $2;
#					my $ORF3 = $3;
#					
#					my $new_line_ORF1 = $_;
#					my $new_line_ORF2 = $_;
#					my $new_line_ORF3 = $_;
#					
#					$new_line_ORF1 =~ s/ORF: [\w\s\.']+, ORF: [\w\s\.']+, ORF: [\w\s\.']+/CDS: $ORF1/;
#					$new_line_ORF2 =~ s/ORF: [\w\s\.']+, ORF: [\w\s\.']+, ORF: [\w\s\.']+/CDS: $ORF2/;
#					$new_line_ORF3 =~ s/ORF: [\w\s\.']+, ORF: [\w\s\.']+, ORF: [\w\s\.']+/CDS: $ORF3/;
#					
#					print NEW_SNPR $new_line_ORF1;
#					print NEW_SNPR $new_line_ORF2;
#					print NEW_SNPR $new_line_ORF3;
#				} elsif($_ =~ /CDS: ([\w\s\.']+), CDS: ([\w\s\.']+), CDS: ([\w\s\.']+)/) {
#					#print "We got here2\n\n";
#					my $ORF1 = $1;
#					my $ORF2 = $2;
#					my $ORF3 = $3;
#					
#					my $new_line_ORF1 = $_;
#					my $new_line_ORF2 = $_;
#					my $new_line_ORF3 = $_;
#					
#					$new_line_ORF1 =~ s/CDS: [\w\s\.']+, CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF1/;
#					$new_line_ORF2 =~ s/CDS: [\w\s\.']+, CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF2/;
#					$new_line_ORF3 =~ s/CDS: [\w\s\.']+, CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF3/;
#					
#					print NEW_SNPR $new_line_ORF1;
#					print NEW_SNPR $new_line_ORF2;
#					print NEW_SNPR $new_line_ORF3;
#				} elsif($_ =~ /ORF: ([\w\s\.']+), CDS: ([\w\s\.']+), CDS: ([\w\s\.']+)/) {
#					#print "We got here2\n\n";
#					my $ORF1 = $1;
#					my $ORF2 = $2;
#					my $ORF3 = $3;
#					
#					my $new_line_ORF1 = $_;
#					my $new_line_ORF2 = $_;
#					my $new_line_ORF3 = $_;
#					
#					$new_line_ORF1 =~ s/ORF: [\w\s\.']+, CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF1/;
#					$new_line_ORF2 =~ s/ORF: [\w\s\.']+, CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF2/;
#					$new_line_ORF3 =~ s/ORF: [\w\s\.']+, CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF3/;
#					
#					print NEW_SNPR $new_line_ORF1;
#					print NEW_SNPR $new_line_ORF2;
#					print NEW_SNPR $new_line_ORF3;
#				} elsif($_ =~ /ORF: ([\w\s\.']+), ORF: ([\w\s\.']+)/) {
#					my $ORF1 = $1;
#					my $ORF2 = $2;
#					
#					my $new_line_ORF1 = $_;
#					my $new_line_ORF2 = $_;
#					
#					$new_line_ORF1 =~ s/ORF: [\w\s\.']+, ORF: [\w\s\.']+/CDS: $ORF1/;
#					$new_line_ORF2 =~ s/ORF: [\w\s\.']+, ORF: [\w\s\.']+/CDS: $ORF2/;
#					
#					print NEW_SNPR $new_line_ORF1;
#					print NEW_SNPR $new_line_ORF2;
#				} elsif($_ =~ /CDS: ([\w\s\.']+), CDS: ([\w\s\.']+)/) {
#					my $ORF1 = $1;
#					my $ORF2 = $2;
#					
#					my $new_line_ORF1 = $_;
#					my $new_line_ORF2 = $_;
#					
#					$new_line_ORF1 =~ s/CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF1/;
#					$new_line_ORF2 =~ s/CDS: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF2/;
#					
#					print NEW_SNPR $new_line_ORF1;
#					print NEW_SNPR $new_line_ORF2;
#				} elsif($_ =~ /ORF: ([\w\s\.']+), CDS: ([\w\s\.']+)/) {
#					my $ORF1 = $1;
#					my $ORF2 = $2;
#					
#					my $new_line_ORF1 = $_;
#					my $new_line_ORF2 = $_;
#					
#					$new_line_ORF1 =~ s/ORF: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF1/;
#					$new_line_ORF2 =~ s/ORF: [\w\s\.']+, CDS: [\w\s\.']+/CDS: $ORF2/;
#					
#					print NEW_SNPR $new_line_ORF1;
#					print NEW_SNPR $new_line_ORF2;
#				} elsif($_ =~ /ORF: ([\w\s\.']+)/) {
#					my $ORF1 = $1;
#					
#					my $new_line_ORF1 = $_;
#					
#					$new_line_ORF1 =~ s/ORF: [\w\s\.']+/CDS: $ORF1/;
#					
#					print NEW_SNPR $new_line_ORF1;
#				} elsif($_ =~ /CDS: ([\w\s\.']+)/) {
#					my $ORF1 = $1;
#					
#					my $new_line_ORF1 = $_;
#					
#					$new_line_ORF1 =~ s/CDS: [\w\s\.']+/CDS: $ORF1/;
#					
#					print NEW_SNPR $new_line_ORF1;
#				} else {
#					print NEW_SNPR $_;
#				}
#				close NEW_SNPR;
#			}
#		}
#		close SNPR_FILE;
##	}
#	
#	#my @new_snp_report_file_names = glob "*_snpg9temp.txt";
##	if (scalar (@new_snp_report_file_names) == 0) {
##		chdir('SNPGenie_Results');
##		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
##		# FILE | PRODUCT | SITE | CODON | WARNING
##		print ERROR_FILE "N/A\tN/A\tN/A\t".
##			"Error processing SNP Reports due to anomalous CDS information. SNPGenie ".
##			"terminated\n";
##		close ERROR_FILE;
##		chdir('..');
##		
##		die "\n\n## WARNING: Error processing SNP Reports due to anomalous CDS ".
##			"information.\n\n##SNPGenie terminated.\n\n";
##	}
#	return $new_file_name;
#	#return @new_snp_report_file_names;
#}


#########################################################################################
# (1) Is passed an array with: [0] $curr_snp_report_name -- the name of the actual SNP 
# Report; and [1] $temp_snp_report_name -- the name of the tempfile to populate
# (2) Processes the SNP report for various anomalies in the Overlapping annotations 
# column, chiefly to get variants that are present in multiple products ON THEIR OWN 
# LINES, and places each processed line into the tempfile
sub populate_tempfile_clc {
	if(scalar @_ != 2) {
		die "\n\n## WARNING: The subroutine populate_tempfile_clc needs exactly 2 ".
			"arguments. SNPgenie terminated.\n\n";
	}
	
	my $curr_snp_report_name = $_[0];
	my $temp_snp_report_name = $_[1];
	
	print "\nConverting $curr_snp_report_name to SNPGenie format...\n";
	
	my @header_names_arr = &get_header_names($temp_snp_report_name,$curr_snp_report_name);
	#print "@header_names_arr";
	#print "\n\n$_\n\n";
	
	my $newline_char = &detect_newline_char($curr_snp_report_name);
	#my $old_newline = $/;
	#$/ = $newline_char;
	my $newline_type;
	
	if($newline_char eq "\r\n") {
		$newline_type = "Windows (CRLF, \\r\\n\)";
	} elsif($newline_char eq "\r") {
		$newline_type = "Mac (CR, \\r\)";
	} elsif($newline_char eq "\n") {
		$newline_type = "Unix (LF, \\n\)";
	}
	
	print "\nIn file $curr_snp_report_name, the newline type is: $newline_type\n";
	
	my $line = 0;
	open(SNPR_FILE, $curr_snp_report_name) or die "\nCould not open SNP Report file ".
		"$curr_snp_report_name\n";
	while (<SNPR_FILE>) {
		if($line == 0) {
			open(TEMP_FILE,">$temp_snp_report_name");
			print TEMP_FILE $_;
			close TEMP_FILE;
			$line++;
		} else {
			open(TEMP_FILE,">>$temp_snp_report_name");
			# First, replace the damned parentheses
			if($_ =~ /CDS: ""([\w\s\.\-\:']+ [\w\s\.\-\:']+)""/) {
				my $replacement = $1;
				#print "\n\nI had a replacement\n\n";
				$_ =~ s/CDS: ""[\w\s\.\-\:']+ [\w\s\.\-\:']+""/CDS: $replacement/;
			}
			
			# To account for prime (') symbols in the ORF names
			if($_ =~ /ORF: ([\w\s\.\-\:']+)[\.\,] ORF: ([\w\s\.\-\:']+)[\.\,] ORF: ([\w\s\.\-\:']+)/) {
				#print "We got here2\n\n";
				my $ORF1 = $1;
				my $ORF2 = $2;
				my $ORF3 = $3;
				
				my $new_line_ORF1 = $_;
				my $new_line_ORF2 = $_;
				my $new_line_ORF3 = $_;
				
				$new_line_ORF1 =~ s/ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+/CDS: $ORF1/;
				$new_line_ORF2 =~ s/ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+/CDS: $ORF2/;
				$new_line_ORF3 =~ s/ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+/CDS: $ORF3/;
				
				print TEMP_FILE $new_line_ORF1;
				print TEMP_FILE $new_line_ORF2;
				print TEMP_FILE $new_line_ORF3;
			} elsif($_ =~ /CDS: ([\w\s\.\-\:']+)[\.\,] CDS: ([\w\s\.\-\:']+)[\.\,] CDS: ([\w\s\.\-\:']+)/) {
				#print "We got here2\n\n";
				my $ORF1 = $1;
				my $ORF2 = $2;
				my $ORF3 = $3;
				
				my $new_line_ORF1 = $_;
				my $new_line_ORF2 = $_;
				my $new_line_ORF3 = $_;
				
				$new_line_ORF1 =~ s/CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF1/;
				$new_line_ORF2 =~ s/CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF2/;
				$new_line_ORF3 =~ s/CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF3/;
				
				print TEMP_FILE $new_line_ORF1;
				print TEMP_FILE $new_line_ORF2;
				print TEMP_FILE $new_line_ORF3;
			} elsif($_ =~ /ORF: ([\w\s\.\-\:']+)[\.\,] CDS: ([\w\s\.\-\:']+)[\.\,] CDS: ([\w\s\.\-\:']+)/) {
				#print "We got here2\n\n";
				my $ORF1 = $1;
				my $ORF2 = $2;
				my $ORF3 = $3;
				
				my $new_line_ORF1 = $_;
				my $new_line_ORF2 = $_;
				my $new_line_ORF3 = $_;
				
				$new_line_ORF1 =~ s/ORF: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF1/;
				$new_line_ORF2 =~ s/ORF: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF2/;
				$new_line_ORF3 =~ s/ORF: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF3/;
				
				print TEMP_FILE $new_line_ORF1;
				print TEMP_FILE $new_line_ORF2;
				print TEMP_FILE $new_line_ORF3;
			} elsif($_ =~ /ORF: ([\w\s\.\-\:']+)[\.\,] ORF: ([\w\s\.\-\:']+)/) {
				my $ORF1 = $1;
				my $ORF2 = $2;
				
				my $new_line_ORF1 = $_;
				my $new_line_ORF2 = $_;
				
				$new_line_ORF1 =~ s/ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+/CDS: $ORF1/;
				$new_line_ORF2 =~ s/ORF: [\w\s\.\-\:']+[\.\,] ORF: [\w\s\.\-\:']+/CDS: $ORF2/;
				
				print TEMP_FILE $new_line_ORF1;
				print TEMP_FILE $new_line_ORF2;
			} elsif($_ =~ /CDS: ([\w\s\.\-\:']+)[\.\,] CDS: ([\w\s\.\-\:']+)/) {
				my $ORF1 = $1;
				my $ORF2 = $2;
				
				my $new_line_ORF1 = $_;
				my $new_line_ORF2 = $_;
				
				$new_line_ORF1 =~ s/CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF1/;
				$new_line_ORF2 =~ s/CDS: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF2/;
				
				print TEMP_FILE $new_line_ORF1;
				print TEMP_FILE $new_line_ORF2;
			} elsif($_ =~ /ORF: ([\w\s\.\-\:']+)[\.\,] CDS: ([\w\s\.\-\:']+)/) {
				my $ORF1 = $1;
				my $ORF2 = $2;
				
				my $new_line_ORF1 = $_;
				my $new_line_ORF2 = $_;
				
				$new_line_ORF1 =~ s/ORF: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF1/;
				$new_line_ORF2 =~ s/ORF: [\w\s\.\-\:']+[\.\,] CDS: [\w\s\.\-\:']+/CDS: $ORF2/;
				
				print TEMP_FILE $new_line_ORF1;
				print TEMP_FILE $new_line_ORF2;
			} elsif($_ =~ /ORF: ([\w\s\.\-\:']+)/) {
				my $ORF1 = $1;
				
				my $new_line_ORF1 = $_;
				
				$new_line_ORF1 =~ s/ORF: [\w\s\.\-\:']+/CDS: $ORF1/;
				
				print TEMP_FILE $new_line_ORF1;
			} elsif($_ =~ /CDS: ([\w\s\.\-\:']+)/) {
				my $ORF1 = $1;
				
				my $new_line_ORF1 = $_;
				
				$new_line_ORF1 =~ s/CDS: [\w\s\.\-\:']+/CDS: $ORF1/;
				
				print TEMP_FILE $new_line_ORF1;
			} else {
				print TEMP_FILE $_;
			}
			close TEMP_FILE;
		}
	}
	close SNPR_FILE;
	#$/ = $old_newline;
}


#########################################################################################
# (1) Is passed an array with: [0] $curr_snp_report_name -- the name of the actual SNP 
# Report; and [1] $temp_snp_report_name -- the name of the tempfile to populate
# (2) Processes the SNP report for various anomalies in the Overlapping annotations 
# column, chiefly to get variants that are present in multiple products ON THEIR OWN 
# LINES, and places each processed line into the tempfile.
# For this GENEIOUS version, we have two missions to accomplish:
# [1] snpgenie_prep_geneious
# [2] snpgenie_geneious_to_clc
sub populate_tempfile_geneious {
	if(scalar @_ != 2) {
		die "\n\n## WARNING: The subroutine populate_tempfile_clc needs exactly 2 ".
			"arguments. SNPgenie terminated.\n\n";
	}
	
	my $curr_snp_report_name = $_[0]; # what we're reading from: the original file
	my $temp_snp_report_name = $_[1]; # what we're populating after processing
	
	# Preparing Geneious' output for SNP Genie by making it tab-delimited and also removing parentheses
	# ONLY CONSIDERS SNPs of LENGTH 1 -- not the multiple-nt linked ones -- but this seem irrelevant

	print "\nConverting $curr_snp_report_name to SNPGenie format...\n";
	
	## PART 1 ## SNPGENIE_PREP_GENEIOUS
	my $temp_processed_snp_report_TEMPLATE = "temp_processed_snp_report_XXXX";
	my ($TEMP_PROCESSED_SNP_REPORT_HANDLE,$temp_processed_snp_report_name) = 
		tempfile($temp_processed_snp_report_TEMPLATE, SUFFIX => ".txt", UNLINK => 1);
#	my ($TEMP_PROCESSED_SNP_REPORT_HANDLE,$temp_processed_snp_report_name) = 
#		tempfile($temp_processed_snp_report_TEMPLATE, SUFFIX => ".txt", UNLINK => 0, OPEN => 0);
	
	#print "My temp processed handle is: $TEMP_PROCESSED_SNP_REPORT_HANDLE\n";
	
	my $newline_char = &detect_newline_char($curr_snp_report_name);
	#my $old_newline = $/;
	#$/ = $newline_char;
	my $newline_type;
	
	if($newline_char eq "\r\n") {
		$newline_type = "Windows (CRLF, \\r\\n\)";
	} elsif($newline_char eq "\r") {
		$newline_type = "Mac (CR, \\r\)";
	} elsif($newline_char eq "\n") {
		$newline_type = "Unix (LF, \\n\)";
	}
	
	my $line = 0;
	open(ORIGINAL_SNP_REPORT,$curr_snp_report_name);
	
	# GO THROUGH EACH LINE OF THE SNP REPORT AND PRINT IT to a tempfile AFTER PROCESSING
	while (<ORIGINAL_SNP_REPORT>) {
		# Get rid of COMMA-DELIMITED VALUES and PARANTHESES
		$_ =~ s/"(\w+),(\w+)"/$1$2/g;
		$_ =~ s/"(\w+),(\w+) -> "/$1$2 -> /g;
		$_ =~ s/"(\w+),(\w+) -> (\w+),(\w+)"/$1\/$2 -> $3\/$4/g;
		$_ =~ s/"(\w+) (\w+), (\(first expressed exon\))"/$1 $2 $3/g;
		
		# REPLACE COMMAS WITH TAB DELIMITERS, which include HEADER requirements
		$_ =~ s/,/\t/g;
		
		print $TEMP_PROCESSED_SNP_REPORT_HANDLE "$_";

	}
	close ORIGINAL_SNP_REPORT;
	#close $TEMP_PROCESSED_SNP_REPORT_HANDLE; # closing it should NOT unlink the file...
	#	CLOSING CAUSES MAJOR PROBLEMS
	seek($TEMP_PROCESSED_SNP_REPORT_HANDLE,0,0);
	
	#$/ = $old_newline;
	
	## PART 2 ## SNPGENIE_GENEIOUS_TO_CLC
	print "\nIn file $curr_snp_report_name, the newline type is: $newline_type\n";
	
	my @header_names_arr = &get_header_names($temp_processed_snp_report_name);
	
	my $index_min;
	my $index_max;
	my $index_type;
	my $index_poly_type;
	my $index_cds_position;
	my $index_change;
	my $index_percent;
	my $index_cov;
	my $index_product;
	
	my $seen_index_min = 0;
	my $seen_index_max = 0;
	my $seen_index_type = 0;
	my $seen_index_poly_type = 0;
	my $seen_index_cds_position = 0;
	my $seen_index_change = 0;
	my $seen_index_percent = 0;
	my $seen_index_cov = 0;
	my $seen_index_product = 0;
	
	# Determine the index of each column
	for (my $i=0; $i<scalar(@header_names_arr); $i++) {
		if ($header_names_arr[$i] =~ /Min \(original sequence\)/) {
			$index_min = $i;
			$seen_index_min = 1;
		} elsif ($header_names_arr[$i] =~ /Minimum/) { 
			$index_min = $i;
			$seen_index_min = 1;
		} elsif ($header_names_arr[$i] =~ /Max \(original sequence\)/) {
			$index_max = $i;
			$seen_index_max = 1;
		} elsif ($header_names_arr[$i] =~ /Maximum/) {
			$index_max = $i;
			$seen_index_max = 1;
		} elsif ($header_names_arr[$i] =~ /Polymorphism Type/) { # Will check the value BEGINS with 'SNP'
			$index_poly_type = $i;
			$seen_index_poly_type = 1;
		} elsif ($header_names_arr[$i] =~ /Type/) { # Will check the value BEGINS with 'SNP'
			$index_type = $i;
			$seen_index_type = 1;
		} elsif ($header_names_arr[$i] =~ /CDS Position/) {
			$index_cds_position = $i;
			$seen_index_cds_position = 1;
		} elsif ($header_names_arr[$i] eq "Change") {
			$index_change = $i;
			$seen_index_change = 1;
		} elsif ($header_names_arr[$i] =~ /Variant Frequency/) {
			$index_percent = $i;
			$seen_index_percent = 1;
		} elsif ($header_names_arr[$i] =~ /Coverage/) {
			$index_cov = $i;
			$seen_index_cov = 1;
		} elsif ($header_names_arr[$i] =~ /product/) {
			$index_product = $i;
			$seen_index_product = 1;
		}
	}
	
	# DIE if one of the headers have not been seen
	if ($seen_index_min == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Minimum\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Minimum\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_max == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Maximum\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Maximum\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_poly_type == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Polymorphism Type\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Polymorphism Type\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_type == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Type\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Type\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_cds_position == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"CDS Position\". Proceed with extreme caution\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		print "\n\n## WARNING: $curr_snp_report_name does not contain the column \"CDS Position\". It is HIGHLY\n".
			"## RECOMMENDED that you include this column, because the \"Minimum\" column is prone to\n".
			"## SUBSTANTIAL ERROR. Proceed with extreme caution.\n\n";
	} elsif ($seen_index_change == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Change\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Change\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_percent == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Variant Frequency\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Variant Frequency\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_cov == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"Coverage\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"Coverage\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_product == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"product\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"product\". SNPGenie terminated.\n\n";	
	}
	
#HEREE	my @product_names_arr = &get_product_names_geneious($temp_processed_snp_report_name,$index_product,$index_type);
	
	#foreach my $prod (@product_names_arr) {
	#	print "Product: $prod\n";
	#}
	
	# BUILD A CLC VERSION
	# Now we want to cycle through the file again in order to BUILD A CLC VERSION file
	# for output
	$line = 0;
	
	$newline_char = &detect_newline_char($temp_processed_snp_report_name);
	#$old_newline = $/;
	#$/ = $newline_char;
	
	open(OUTFILE,">>$temp_snp_report_name");
	#open(CURRINFILE, $temp_processed_snp_report_name) or die "Could not open $TEMP_PROCESSED_SNP_REPORT_HANDLE\n";
	#while(<CURRINFILE>) {
	#open($TEMP_PROCESSED_SNP_REPORT_HANDLE) or die "Could not open $TEMP_PROCESSED_SNP_REPORT_HANDLE\n";
	while(<$TEMP_PROCESSED_SNP_REPORT_HANDLE>) {
		if($line == 0) {
			#print "$_";
			if(!($_ =~/\t/)) {
				die "\n\n## WARNING:\n## The processed SNP Report $curr_snp_report_name is ".
					"not TAB-delimited (\\t), or there is only one column.\n\n";
			}
			
			my $clc_format_header = "File\tReference Position\tCDS Position\tType\tReference\t".
				"Allele\tCount\tCoverage\tFrequency\tOverlapping annotations\n";
			print OUTFILE "$clc_format_header";
			#print "$clc_format_header";
			$line++;
			
		} else {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			my @line_arr = split(/\t/,$_,-1);
			
			# ONLY DO THE THING IF THE GENEIOUS TYPE IS "Polymorphism"; not CDS.
			my $type = $line_arr[$index_type];
			#print "$type\n";
			if($type =~ 'Polymorphism') {
				#print "Type is Polymorphism\n";
				# Save the pure records; the GENEIOUS ORDER is:
				# min => 350
				# max => 350
				# type => Polymorphism OR CDS
				# cds_position => 480
				# change => G -> A
				# percent => 20.56%
				# cov => 1313
				# poly_type => SNP (transition) or SNP (transversion) or Substitution or ""
				# product => gag
				
				my $poly_type = $line_arr[$index_poly_type];
				my $percent = $line_arr[$index_percent];
				
				if(($poly_type =~ /SNP/ || $poly_type =~ /Substitution/) && ($percent > 0)) {
					#print "Poly Type is SNP\n";
					my $min = $line_arr[$index_min];					
					my $max = $line_arr[$index_max];
					my $cds_position = $line_arr[$index_cds_position]; # 1-based
					my $change = $line_arr[$index_change];
					my $coverage = $line_arr[$index_cov];
					my $product_name = $line_arr[$index_product];
					
					if($min eq '...') {
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						# FILE | PRODUCT | SITE | WARNING
						print ERROR_FILE "$curr_snp_report_name\t$product_name\t$min\t".
							"In this Geneious SNP report, a site in the Minimum column is reported as '...'. ".
							"SNPGenie will try to infer the stie from the CDS Position. However, if possible, ".
							"please try to produce a SNP report without errors.\n";
						close ERROR_FILE;
						chdir('..');
						
						print "\n\n## WARNING: In file $curr_snp_report_name, product $product_name, site $min,\n".
							"## a site in the Minimum column is reported as '...'. SNPGenie will try to infer the stie\n".
							"## from the CDS Position. However, if possible, please try to produce a SNP report without errors.\n";
					}
					
					# New segments approach
					# We want to USE CDS POSITION to determine the "Minimum" value,
					# because the given Minimum value is so often incorrect. As in:
					# SWITCHED BACK FROM THIS FOR THE SIV DATA BECAUSE OF MULTI-SEGMENT PRODUCTS
					#my $min = ($line_arr[$index_cds_position] + $cds_coordinates{$product_name}->{start} - 1);
					# SO we need to know the starting position of the relevant product. 
					# Do we have it yet? We DO have the product name.
					my @cds_coords_arr; # product start site-based
					# IF there is CDS data AND WE CAN, change the Minimum to match
					if(($product_name ne '') && ($cds_position > 0)) { # STRINGS are 0
						@cds_coords_arr = @{$product_coordinates_harr{$product_name}->{product_coord_arr}};
						
						#print "\n$product_name @cds_coords_arr\n";
						
						my $num_segments = (@cds_coords_arr / 2);
						my %this_product_starts;
						my %this_product_stops;
						for(my $i=1; $i<=scalar(@cds_coords_arr); $i++) {
							$this_product_starts{$i} = $cds_coords_arr[2*$i-2];
							$this_product_stops{$i} = $cds_coords_arr[2*$i-1];
						}
						
						my $cumul_length_prev_product = 0;
						# For each segment, step through each site, add to length and sequence
						for(my $i = 1; $i <= $num_segments; $i++) {
							
							my $segment_start_site = $this_product_starts{$i};
							my $segment_stop_site = $this_product_stops{$i};
							my $current_length = $segment_stop_site - $segment_start_site + 1;
							
							if($cds_position > $cumul_length_prev_product && $cds_position <= ($cumul_length_prev_product + $current_length)) { # within this segment's range of sites
								my $actual = ($segment_start_site + ($cds_position - $cumul_length_prev_product) - 1);
								
								if($min != $actual) {
									
									chdir('SNPGenie_Results');
									open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
									# FILE | PRODUCT | SITE | WARNING
									print ERROR_FILE "$curr_snp_report_name\t$product_name\t$min\t".
										"There is a conflict between the Minimum ($min) and the actual site ($actual) implied by the CDS Position. The ".
										"latter has been used to determine the correct site; please verify these results\n";
									close ERROR_FILE;
									chdir('..');
									
									print "\n\n## WARNING: In file $curr_snp_report_name, product $product_name, site $min,\n".
									"## there is a conflict between the Minimum ($min) and the actual site ($actual) implied by the CDS Position.\n".
									"## The latter has been used to determine the correct site; please verify these results.\n";
									
									$min = $actual;
								}
							}
							
							$cumul_length_prev_product += $current_length;
							
						}			
					}
					
					# Extract the needed values; the CLC ORDER will be:
					# file_nm
					# ref_pos => 350
					# cds_pos => 28 # WE ARE ADDING THIS TO ENABLE ERROR CHECKING LATER, *IF*
						# WE ARE IN GENEIOUS MODE
					# type => SNV OR MNV OR Insertion OR etc.
					# ref => G
					# allele => A
					# count => 270
					# cov => 1313
					# freq => 20.56
					# over_annot => CDS: gag
				
					# REF_POS
					# In case there's a less-than sign
					if($line_arr[$index_min] =~ /\<(\d+)/) {
						$min = $1;
						#print "We found the min $min\n";
					}
					my $ref_pos = $min;
					
					# TYPE
					# Find length, nucleotides involved, etc.
					my $reference_nts;
					my $variant_nts;
					my $reference_length;
					my $variant_length;
					my $is_change = 0;
					if($change =~/(\w+) -> (\w+)/) {
						#print "We found a change\n";
						$reference_nts = $1;
						$variant_nts = $2;
						$reference_length = length($reference_nts);
						$variant_length = length($variant_nts);
						if($reference_nts ne $variant_nts) {
							$is_change = 1;
						}
					}
					
					if($is_change == 1 && $reference_length == $variant_length) {
						#print "We verified a change\n";
						my $clc_type;
						if($reference_length == 1) {
							$clc_type = 'SNV';
						} elsif($reference_length > 1) {
							$clc_type = 'MNV';
						}
						# what others must we account for?
						
						my $confirm_max = ($min + $variant_length - 1);
						
						#if($max == 0) {
						#	print "Yes, the string '...' does == 0.\n";
						#}
						
						if(($max != 0) && ($max != $confirm_max)) {
							#die "\n\n## WARNING: In the SNP Report $curr_snp_report_name, site $min, variant $change,\n".
							#	"## the Maximum should be $confirm_max but is $max. This error often accompanies\n".
							#	"## mistakes in the site coordinate (Minimum) column, so you may wish to verify\n".
							#	"## that the values in the Minimum column match those expected given the CDS Position.\n".
							#	"## Please correct such Geneious errors before proceeding. SNPGENIE TERMINATED.\n\n";
							
							print "\n\n## WARNING: In the SNP Report $curr_snp_report_name, site $min, variant $change,\n".
								"## the Maximum should be $confirm_max but is $max. This error often accompanies\n".
								"## mistakes in the site coordinate (Minimum) column, so you may wish to verify\n".
								"## that the values in the Minimum column match those expected given the CDS Position.\n".
								"## The latter has been used to determine the correct site; please verify these results.\n\n";
						}
						
						my $product_print;
						if($product_name ne '') { 
							#print "Product name isn't blank\n";
							$product_print = "CDS: $product_name";
						} else {
							#print "Product name is blank\n";
							$product_print = $product_name;
						}
						
						#print "Product name is: $product_name\nBlah";
						
						# FREQ and COUNT
						my $freq = $line_arr[$index_percent] / 100; # this automatically trims the "%"
						my $count = ($line_arr[$index_cov] * $freq);
						
						my $this_line = "$curr_snp_report_name\t$ref_pos\t$cds_position\t$clc_type\t$reference_nts\t".
							"$variant_nts\t$count\t$coverage\t$percent\t$product_print\n";
						
						print OUTFILE "$this_line";
						#print "$this_line";
					} # If it's a change
				} # If poly type contains 'SNP' or 'Substitution'
			} # If type contains 'Polymorphism'
			
			$line ++;
		} # $line does not equal 0
	} # END WHILE LOOP
	#close CURRINFILE;
	close $TEMP_PROCESSED_SNP_REPORT_HANDLE;
	close OUTFILE;
	#$/ = $old_newline;
	unlink $temp_processed_snp_report_name;
}


## DIFFERENT FORMATS DIVEGRE WITHIN
#########################################################################################
sub populate_tempfile_vcf {
	if(scalar @_ != 3) {
		die "\n\n## WARNING: The subroutine populate_tempfile_clc needs exactly 3 ".
			"arguments. SNPgenie terminated.\n\n";
	}
	
	my $curr_snp_report_name = $_[0]; # what we're reading from: the original file
	my $temp_snp_report_name = $_[1]; # what we're populating after processing
	my $gtf_file_nm = $_[2];
	
	print "\nConverting $curr_snp_report_name to SNPGenie format...\n";

	my $newline_char = &detect_newline_char($curr_snp_report_name);
	#$/ = $newline_char;
	my $newline_type;
	
	if($newline_char eq "\r\n") {
		$newline_type = "Windows (CRLF, \\r\\n\)";
	} elsif($newline_char eq "\r") {
		$newline_type = "Mac (CR, \\r\)";
	} elsif($newline_char eq "\n") {
		$newline_type = "Unix (LF, \\n\)";
	}
	
	print "\nIn file $curr_snp_report_name, the newline type is: $newline_type\n";
	
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
	#my $index_sample1;
	
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
	#my $seen_index_sample1;
	
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
	my $header_line;
	open (ORIGINAL_SNP_REPORT, $curr_snp_report_name);
	while (<ORIGINAL_SNP_REPORT>) {	
		chomp;
		# CHOMP for 3 operating systems
		if($_ =~ /\r\n$/) {
			$_ =~ s/\r\n//;
		} elsif($_ =~ /\r$/) {
			$_ =~ s/\r//;
		} elsif($_ =~ /\n$/) {
			$_ =~ s/\n//;
		}
			
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
	my @header_arr = split("\t",$header_line,-1);
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
		} #elsif ($header_arr[$i] =~ /sample1/) {
#			$index_sample1 = $i;
#			$seen_index_sample1 = 1;
#		}
	}
	
#	# Count the number of samples
#	my $num_samples = 0;
#	my @sample_names;
#	for(my $samp_index = $index_sample; $samp_index < scalar(@header_arr); $samp_index++) {
#		if($header_arr[$samp_index] =~ /[\w\d\_\.]+/) {
#			$num_samples++;
#			push(@sample_names, $&);
#		}
#	}
#	
#	print "\nThis VCF contains $num_samples samples: @sample_names\n\n";
	
	# DIE if one of the headers have not been seen
	if ($seen_index_chrom == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"CHROM\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"CHROM\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_pos == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"POS\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"POS\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_id == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"ID\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"ID\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_ref == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"REF\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"REF\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_alt == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"ALT\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"ALT\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_qual == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"QUAL\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"QUAL\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_filter == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"FILTER\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FILTER\". SNPGenie terminated.\n\n";	
	} elsif ($seen_index_info == 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
			"Does not contain the column header \"INFO\". SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"INFO\". SNPGenie terminated.\n\n";	
	} #elsif ($seen_index_format == 0) {
#		chdir('SNPGenie_Results');
#		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
#		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
#			"Does not contain the column header \"FORMAT\". SNPGenie terminated.\n";
#		close ERROR_FILE;
#		chdir('..');
#		
#		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FORMAT\". SNPGenie terminated.\n\n";	
#	} #elsif ($seen_index_sample1 == 0) {
#		chdir('SNPGenie_Results');
#		open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
#		print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
#			"Does not contain the column header \"sample1\". SNPGenie terminated.\n";
#		close ERROR_FILE;
#		chdir('..');
#		
#		#unlink $curr_snp_report_name;
#		
#		die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"sample1\". SNPGenie terminated.\n\n";	
#	}
	
	# NEED TO BUILD A HASH WITH keys as ALL PRODUCT POSITIONS IN THE GENOME, and values being an array
	# of all PRODUCTS overlapping that position. Do we want all the products just on this strand,
	# or all product on both strands (sense and antisense)? Well, we are running JUST ONE VCF file,
	# while running SNPGenie with a single GTF file for the forward, then again another for the 
	# reverse strands. ON THIS RUN, however, we're just concerned with assigning THIS STRAND'S
	# PRODUCTS to the right sites. We can go in later, after creating the tempfile, to get
	# information about how many products actually overlap each site. In other words, we're only
	# going to worry about this strand's GTF file for the purposes of this subroutine.
	
	# Loop through GTF and store (1) positions having products in %positions_with_product_hash
	# and (2) product names in keys of %products_hash
	# For this, the NUMBER of segments for each gene don't matter.
	my %positions_with_product_hash; # $positions_with_product_hash{SITE}->{@PRODUCTS}
	my %products_hash; # $products_hash{PRODUCT_NAME}->{1}; this is to check if it yet exists
	open (GTF_INFILE, $gtf_file_nm);
	while (<GTF_INFILE>) {
		chomp;
		# CHOMP for 3 operating systems
		if($_ =~ /\r\n$/) {
			$_ =~ s/\r\n//;
		} elsif($_ =~ /\r$/) {
			$_ =~ s/\r//;
		} elsif($_ =~ /\n$/) {
			$_ =~ s/\n//;
		}
		
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Make sure it's on the + strand
			my $product;
			if($_ =~ /\s*gene_id\s*\"gene\:([\w\s\.\-\:']+)\"/) {
				$product = $1;
				
				if ((! exists $products_hash{$product}) && ($product ne '')) {
					$products_hash{$product} = 1;
				}
			} elsif($_ =~ /\s*gene_id\s*\"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
				$product = $1;
				
				if ((! exists $products_hash{$product}) && ($product ne '')) {
					$products_hash{$product} = 1;
				}
			} elsif($_ =~ /\s*gene_id\s*\"([\w\s\.\-\:']+)\"/) {
				$product = $1;
				
				if ((! exists $products_hash{$product}) && ($product ne '')) {
					$products_hash{$product} = 1;
				}
			}
			
			# This is not redundant, because we're getting the names now
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
	
	my @product_names_arr = sort(keys %products_hash);
	
	#foreach (@product_names_arr) {
	#	print "$_\n";
	#}
	
	# INCLUDES "hypothetical protein" -- this not a specific tag?
	#foreach my $prod (@product_names_arr) {
	#	print "Product: $prod\n";
	#}
	
	# Now we want to cycle through the file again in order to build a CLC version file
	# for output
	#open(OUTFILE,">>snpgenie_tempfile.txt");
	open(TEMP_FILE,">>$temp_snp_report_name");
	my $clc_format_header = "File\tReference Position\tType\tReference\t".
			"Allele\tCount\tCoverage\tFrequency\tOverlapping annotations\n";
	print TEMP_FILE "$clc_format_header";
	
	open (ORIGINAL_SNP_REPORT, $curr_snp_report_name);
	while (<ORIGINAL_SNP_REPORT>) {
		unless(/^#/) { # lines that begins with "##" or "#" are metadata
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			my @line_arr = split(/\t/,$_,-1);
			
			# ONLY DO THE THING IF THE GENEIOUS TYPE IS "Polymorphism"; not CDS.
			my $id = $line_arr[$index_id];
			#print "$type\n";
			#if($id eq '.') { # It's a SNV
			
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
			
			my $variant1 = 0;
			my $variant2 = 0;
			my $variant3 = 0;
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
			# Let's also add to this the EXCLUSION of lines which have 0 variants;
			# this should save a considerable amount of processing time later.
			my $product_entry = '';
			if($equal_lengths_variants) {
				if($variant3) { # THERE ARE THREE VARIANTS -- ADD A FLAG!
					
					##SAMVCF VCF FORMAT #1
					if($vcfformat == 1) {
					#if($info_value =~ /NS=(\d+)/) { # We've got a VCF SUMMARIZING INDIVIDUALS
						if($warn_file_type_not_supported == 0) {
							print "\n### WARNING: VCF FORMAT TYPE 1 IS NOT FULLY SUPPORTED ###\n";
							$warn_file_type_not_supported ++;
						}
						
						my $num_samples;
						if($info_value =~ /NS=(\d+)/) {
							$num_samples = $1;
						}
						
						my $variant_freq1;
						my $variant_freq2;
						my $variant_freq3;
						if($info_value =~ /AF=([\d\.e\-]+),([\d\.e\-]+),([\d\.e\-]+)/) {
							$variant_freq1 = $1;
							$variant_freq2 = $2;
							$variant_freq3 = $3;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated.\n\n";	
						}
						# elsif($info_value =~ /AF=([\d\.\e\-]+),([\d\.\e\-]+)/) {
						#	$variant_freq1 = $1;
						#	$variant_freq2 = $2;
						#} elsif($info_value =~ /AF=([\d\.\e\-]+)/) {
						#	$variant_freq1 = $1;
						#} 
						
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
										"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\t$product_entry\n";
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
								
								if($variant_freq3 > 0) {
									my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant3\t$variant_count3\t$num_samples\t$variant_pct3\t$product_entry\n";
									
									print TEMP_FILE "$this_line3";
								}
							}
						} else { # no products at this site, so no product entries (blank)
							# PRINT 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$num_samples\t$variant_pct3\t$product_entry\n";
								
								print TEMP_FILE "$this_line3";
							}
						}
						
					##SAMVCF VCF FORMAT #2	
					} elsif($vcfformat == 2) { # We've got a VCF of POOL ##SAMVCF
						
						my $coverage;
						my $variant_freq1;
						my $variant_freq2;
						my $variant_freq3;
						
						if($info_value =~ /DP=(\d+)/) { # We've got a VCF of POOL
							$coverage = $1;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated.\n\n";	
						}
						
						if($info_value =~ /AF=([\d\.e\-]+),([\d\.e\-]+),([\d\.e\-]+)/) { # We've got a VCF of POOL
							$variant_freq1 = $1;
							$variant_freq2 = $2;
							$variant_freq3 = $3;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated.\n\n";	
						}
						#elsif($info_value =~ /AF=([\d\.\e\-]+),([\d\.\e\-]+)/) {
						#	$variant_freq1 = $1;
						#	$variant_freq2 = $2;
						#} elsif($info_value =~ /AF=([\d\.\e\-]+)/) {
						#	$variant_freq1 = $1;
						#} 
						
						# COUNTS and FREQS
						my $variant_count1 = $variant_freq1 * $coverage;
						my $variant_count2 = $variant_freq2 * $coverage;
						my $variant_count3 = $variant_freq3 * $coverage;
						
						my $alt_count = ($variant_count1 + $variant_count2 + $variant_count3);
						my $ref_count = $coverage - $alt_count;
						
						my $variant_pct1 = (100 * $variant_freq1);
						my $variant_pct2 = (100 * $variant_freq2);
						my $variant_pct3 = (100 * $variant_freq3);
						
						$product_entry = '';
						if(@this_site_products) {
							foreach my $product (@this_site_products) {
								$product_entry = 'CDS: ' . $product;
								
								# PRINT UP TO 3 LINES TO FILE
								if($variant_freq1 > 0) {
									my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
								
								if($variant_freq3 > 0) {
									my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
									
									print TEMP_FILE "$this_line3";
								}
							}
						} else {
							# PRINT UP TO 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
								
								print TEMP_FILE "$this_line3";
							}
						}
						
					##SAMVCF VCF FORMAT #3
					} elsif($vcfformat == 3) {
					#} elsif($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) { # We've got a VCF of POOL
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						print ERROR_FILE "$curr_snp_report_name\t". $this_site_products[0] .
							"\t$ref_pos\t".
							"Site has three variants in a pooled VCF file. Variant frequencies".
							" have been approximated\n";
						close ERROR_FILE;
						chdir('..');
						
						my $fwd_ref_reads;
						my $rev_ref_reads;
						my $fwd_alt_reads;
						my $rev_alt_reads;
						
						if($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) {
							# These are high-quality reads, so may be less that the actual coverage
							$fwd_ref_reads = $1;
							$rev_ref_reads = $2;
							$fwd_alt_reads = $3;
							$rev_alt_reads = $4;
						}
						
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
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
								
								if($variant_freq3 > 0) {
									my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
									
									print TEMP_FILE "$this_line3";
								}
							}
						} else {
							# PRINT 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
								
								print TEMP_FILE "$this_line3";
							}
						}
						
					##SAMVCF VCF FORMAT #4
					} elsif($vcfformat == 4) {	
					#} elsif($format_value =~ /AD/) { # We've got a VCF of POOL; we need AD and DP
						# Die if there wasn't a FORMAT column, necessary for this format
						if ($seen_index_format == 0) {
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"Does not contain the column header \"FORMAT\". SNPGenie terminated.\n";
							close ERROR_FILE;
							chdir('..');
		
							die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FORMAT\". SNPGenie terminated.\n\n";	
						}
						
						# Find out how many ":" appear before AD, if any
						my $prior_to_AD;
						
						if($format_value =~ /AD/) {
							$prior_to_AD = $`;
						}
						
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
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"VCF file $curr_snp_report_name contains AD but not DP data. SNPGenie terminated.\n";
							close ERROR_FILE;
							chdir('..');
							
							die "\n\n## WARNING: $curr_snp_report_name contains AD but not DP data. SNPGenie terminated.\n\n";	 
						}
						
						# EXTRACT the VALUE of AD
						my @sample_value_arr = split(/\:/,$sample_value,-1);
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
						if ($DP < $coverage) {
							warn "\n## WARNING: In $curr_snp_report_name site $ref_pos".
									",\n## the reads total should ".
									"equal the coverage ($DP) but is instead: $coverage.".
									"\n## The reads total has been used. Please verify your data.\n";
							
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							# FILE | PRODUCT | SITE | CODON | WARNING
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"The reads total should equal the coverage ($DP) but is instead: ".
								"$coverage. The reads total has been used. Please verify your data.\n";
							close ERROR_FILE;
							chdir('..');
						} # it seems common that DP > coverage, because some reads get filtered
								
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
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
								
								if($variant_freq3 > 0) {
									my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
									
									print TEMP_FILE "$this_line3";
								}
							}
						} else {
							# PRINT 3 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
							
							if($variant_freq3 > 0) {
								my $this_line3 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant3\t$variant_count3\t$coverage\t$variant_pct3\t$product_entry\n";
								
								print TEMP_FILE "$this_line3";
							}
						}
					}
					
				} elsif($variant2) { # THERE ARE TWO VARIANTS -- ADD A FLAG!
					
					##SAMVCF VCF FORMAT #1
					if($vcfformat == 1) {
					#if($info_value =~ /NS=(\d+)/) { # We've got a VCF summarizing INDIVIDUALS
						if($warn_file_type_not_supported == 0) {
							print "\n### WARNING: SNP REPORT FILE TYPE NOT FULLY SUPPORTED ###\n";
							$warn_file_type_not_supported ++;
						}
						
						my $num_samples;
						if($info_value =~ /NS=(\d+)/) {
							$num_samples = $1;
						}
						
						my $variant_freq1;
						my $variant_freq2;
						if($info_value =~ /AF=([\d\.e\-]+),([\d\.e\-]+)/) {
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
										"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\t$product_entry\n";
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
							}
						} else {
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$num_samples\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
						}
						
					##SAMVCF VCF FORMAT #2	
					} elsif($vcfformat == 2) { # We've got a VCF of POOL ##SAMVCF
						
						my $coverage;
						my $variant_freq1;
						my $variant_freq2;
						
						if($info_value =~ /DP=(\d+)/) { # We've got a VCF of POOL
							$coverage = $1;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated.\n\n";	
						}
						
						if($info_value =~ /AF=([\d\.e\-]+),([\d\.e\-]+)/) { # We've got a VCF of POOL
							$variant_freq1 = $1;
							$variant_freq2 = $2;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated.\n\n";	
						}
						
						# COUNTS and FREQS
						my $variant_count1 = $variant_freq1 * $coverage;
						my $variant_count2 = $variant_freq2 * $coverage;
						
						my $alt_count = ($variant_count1 + $variant_count2);
						my $ref_count = $coverage - $alt_count;
						
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
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
							}
						} else {
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
								
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
						}
					
					##SAMVCF VCF FORMAT #3
					} elsif($vcfformat == 3) { # We've got a VCF of POOL ##SAMVCF
					#} elsif($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) { # We've got a VCF of POOL
						chdir('SNPGenie_Results');
						open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
						print ERROR_FILE "$curr_snp_report_name\t". $this_site_products[0] .
							"\t$ref_pos\t".
							"Site has two variants in a pooled VCF file. Variant frequencies".
							" have been approximated\n";
						close ERROR_FILE;
						chdir('..');
						
						my $fwd_ref_reads;
						my $rev_ref_reads;
						my $fwd_alt_reads;
						my $rev_alt_reads;
						
						if($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) {
							# These are high-quality reads, so may be less that the actual coverage
							$fwd_ref_reads = $1;
							$rev_ref_reads = $2;
							$fwd_alt_reads = $3;
							$rev_alt_reads = $4;
						}
						
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
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
							}
						} else {
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
								
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
						}
					
					##SAMVCF VCF FORMAT #4, variant length 2
					} elsif($vcfformat == 4) {
					#} elsif($format_value =~ /AD/) { # We've got a VCF of POOL; we need AD and DP
						# Die if there wasn't a FORMAT column, necessary for this format
						if ($seen_index_format == 0) {
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"Does not contain the column header \"FORMAT\". SNPGenie terminated.\n";
							close ERROR_FILE;
							chdir('..');
		
							die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FORMAT\". SNPGenie terminated.\n\n";	
						}
						
						# Find out how many ":" appear before AD, if any
						my $prior_to_AD;
						
						if($format_value =~ /AD/) {
							$prior_to_AD = $`;
						}	
					
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
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"VCF file $curr_snp_report_name contains AD but not DP data. SNPGenie terminated.\n";
							close ERROR_FILE;
							chdir('..');
							
							die "\n\n## WARNING: $curr_snp_report_name contains AD but not DP data. SNPGenie terminated.\n\n";	 
						}
						
						# GENERATE a regex for what comes before the VALUE of AD
						#my $prior_to_AD_regex;
						#for (my $i=0; $i<$colon_count_before_AD; $i++) {
						#	$prior_to_AD_regex .= "[\d\w\/\,]+\:";
						#}
						
						# EXTRACT the VALUE of AD
						my @sample_value_arr = split(/\:/,$sample_value,-1);
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
						if ($DP < $coverage) {
							warn "\n## WARNING: In $curr_snp_report_name site $ref_pos".
									",\n## the reads total should ".
									"equal the coverage ($DP) but is instead: $coverage.".
									"\n## The reads total has been used. Please verify your data.\n";
							
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							# FILE | PRODUCT | SITE | CODON | WARNING
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"The reads total should equal the coverage ($DP) but is instead: ".
								"$coverage. The reads total has been used. Please verify your data.\n";
							close ERROR_FILE;
							chdir('..');
						} # it seems common that DP > coverage, because some reads get filtered
								
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
									
									print TEMP_FILE "$this_line1";
								}
								
								if($variant_freq2 > 0) {
									my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
									
									print TEMP_FILE "$this_line2";
								}
							}
						} else {
							# PRINT 2 LINES TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
							
							if($variant_freq2 > 0) {
								my $this_line2 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant2\t$variant_count2\t$coverage\t$variant_pct2\t$product_entry\n";
								
								print TEMP_FILE "$this_line2";
							}
						}
					} # variant length 2, VCF FORMAT 4
				
				} elsif($variant1) { # THERE IS ONE VARIANT -- no flag needed
					
					##SAMVCF VCF FORMAT #1
					if($vcfformat == 1) {
					#if($info_value =~ /NS=(\d+)/) { # We've got a VCF summarizing INDIVIDUALS
						if($warn_file_type_not_supported == 0) {
							print "\n### WARNING: SNP REPORT FILE TYPE NOT FULLY SUPPORTED ###\n";
							$warn_file_type_not_supported ++;
						}
						
						my $num_samples;
						if($info_value =~ /NS=(\d+)/) {
							$num_samples = $1;
						}
						
						my $variant_freq1;
						if($info_value =~ /AF=([\d\.e\-]+)/) {
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
										"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\t$product_entry\n";
							
									print TEMP_FILE "$this_line1";
								}
							}
						} else {
							# PRINT 1 LINE TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$num_samples\t$variant_pct1\t$product_entry\n";
							
								print TEMP_FILE "$this_line1";
							}
						}
					
					##SAMVCF VCF FORMAT #2	
					} elsif($vcfformat == 2) { # We've got a VCF of POOL ##SAMVCF
						
						my $coverage;
						my $variant_freq1;
						
						if($info_value =~ /DP=(\d+)/) { # We've got a VCF of POOL
							$coverage = $1;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated 1.\n\n";	
						}
						
						if($info_value =~ /AF=([\d\.e\-]+)/) { # We've got a VCF of POOL
							$variant_freq1 = $1;
						} else {
							die "\n\n## WARNING: $curr_snp_report_name does not conform to ".
								"VCF format $vcfformat. SNPGenie terminated 2.\n\n";	
						}
						
						# COUNTS and FREQS
						my $variant_count1 = $variant_freq1 * $coverage;
						
						my $alt_count = ($variant_count1);
						my $ref_count = $coverage - $alt_count;
						
						my $variant_pct1 = (100 * $variant_freq1);
						
						$product_entry = '';
						if(@this_site_products) {
							foreach my $product (@this_site_products) {
								$product_entry = 'CDS: ' . $product;
								
								# PRINT 1 LINE TO FILE
								if($variant_freq1 > 0) {
									my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
										"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
									print TEMP_FILE "$this_line1";
								}
							}
						} else {
							# PRINT 1 LINE TO FILE; $product_entry remains BLANK
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
						}
					
					##SAMVCF VCF FORMAT #3
					} elsif($vcfformat == 3) { # We've got a VCF of POOL ##SAMVCF
					#} elsif($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) { # We've got a VCF of POOL
						
						my $fwd_ref_reads;
						my $rev_ref_reads;
						my $fwd_alt_reads;
						my $rev_alt_reads;
						
						if($info_value =~ /DP4=(\d+),(\d+),(\d+),(\d+)/) {
							# These are high-quality reads, so may be less that the actual coverage
							$fwd_ref_reads = $1;
							$rev_ref_reads = $2;
							$fwd_alt_reads = $3;
							$rev_alt_reads = $4;
						}
						
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
								
									print TEMP_FILE "$this_line1";
								}
							}
						} else {
							# PRINT 1 LINE TO FILE; $product_entry remains BLANK
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
						}
						
					##SAMVCF VCF FORMAT #4 variant length 1
					} elsif($vcfformat == 4) {
					#} elsif($format_value =~ /AD/) { # We've got a VCF of POOL; we need AD and DP
						# Die if there wasn't a FORMAT column, necessary for this format
						if ($seen_index_format == 0) {
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"Does not contain the column header \"FORMAT\". SNPGenie terminated.\n";
							close ERROR_FILE;
							chdir('..');
		
							die "\n\n## WARNING: $curr_snp_report_name does not contain the column header \"FORMAT\". SNPGenie terminated.\n\n";	
						}
						
						# Find out how many ":" appear before AD, if any
						my $prior_to_AD;
						
						if($format_value =~ /AD/) {
							$prior_to_AD = $`;
						}		

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
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"VCF file $curr_snp_report_name contains AD but not DP data. SNPGenie terminated.\n";
							close ERROR_FILE;
							chdir('..');
							
							die "\n\n## WARNING: $curr_snp_report_name contains AD but not DP data. SNPGenie terminated.\n\n";	 
						}
						
						# EXTRACT the VALUE of AD
						my @sample_value_arr = split(/\:/,$sample_value,-1);
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
						if ($DP < $coverage) {
							warn "\n## WARNING: In $curr_snp_report_name site $ref_pos".
									",\n## the reads total should ".
									"equal the coverage ($DP) but is instead: $coverage.".
									"\n## The reads total has been used. Please verify your data.\n";
							
							chdir('SNPGenie_Results');
							open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
							# FILE | PRODUCT | SITE | CODON | WARNING
							print ERROR_FILE "$curr_snp_report_name\tNA\tNA\t".
								"The reads total should equal the coverage ($DP) but is instead: ".
								"$coverage. The reads total has been used. Please verify your data.\n";
							close ERROR_FILE;
							chdir('..');
						} # it seems common that DP > coverage, because some reads get filtered
						
						####################
						if ($coverage == 0) {
							warn "\n### WARNING: coverage is 0 in $curr_snp_report_name site $ref_pos.\n\n";
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
									
									print TEMP_FILE "$this_line1";
								}
							}
						} else {
							# PRINT 1 LINE TO FILE
							if($variant_freq1 > 0) {
								my $this_line1 = "$curr_snp_report_name\t$ref_pos\t$clc_type\t$reference_nts\t".
									"$variant1\t$variant_count1\t$coverage\t$variant_pct1\t$product_entry\n";
								
								print TEMP_FILE "$this_line1";
							}
						}
					} # end of AD format
				} # end of length 1 variant
			} # equal length variants
		} # don't begin with # (or ##)
	} # done with SNP report
	close ORIGINAL_SNP_REPORT;
	close TEMP_FILE;
}

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
#sub get_product_names_geneious {
#	my ($curr_snp_report_filename,$index_product,$index_type) = @_;
#	my $line = 0;
#	my %products_hash;
#	open (CURRINFILE, $curr_snp_report_filename);
#	while (<CURRINFILE>) {
#		if ($line == 0) {
#			$line++;
#		} else {
#			#chomp;
#			# CHOMP for 3 operating systems
#			if($_ =~ /\r\n$/) {
#				$_ =~ s/\r\n//;
#			} elsif($_ =~ /\r$/) {
#				$_ =~ s/\r//;
#			} elsif($_ =~ /\n$/) {
#				$_ =~ s/\n//;
#			}
#			
#			if($_ =~ /1.40E-19/) {
#				print "\n$_\n\n";
#			}
#			
#			my @line_arr = split(/\t/,$_);
#			my $product = $line_arr[$index_product];
#			my $type = $line_arr[$index_type];
#			
#			if (! exists $products_hash{$product} && 
#				(($type eq 'Polymorphism') || ($type eq 'CDS')) &&
#				($product ne '')) {
#				
#				$products_hash{$product} = 1;
#			}
#		}
#	}
#	close CURRINFILE;
#	my @product_names = keys %products_hash;
#	return @product_names;
#}

#########################################################################################
# (1) Is passed an array of all files to process, assumed to be in the working directory
# (2) Processes those SNP reports for various anomalies in the Overlapping annotations 
# column, chiefly to get variants that are present in multiple products ON THEIR OWN 
# LINES. An inelegant solution, but it beats re-doing the whole algorithm.
# (3) Obtains and returns all file names in current directory ending in _snpg9temp.txt
#sub remove_tempfiles {
#	my @tempfile_names_arr = glob "*_snpg9temp.txt";
#	unlink @tempfile_names_arr;
#}

#########################################################################################
# Obtains the file in the current directory ending in .gtf
sub get_cds_file_name { 
	my $cds_file_name = glob "*\.gtf";
	if ($cds_file_name eq '') {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		# FILE | PRODUCT | SITE | CODON | WARNING
		print ERROR_FILE "N/A\tN/A\tN/A\t".
			"No .gtf file in the working directory for product information. SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		die "\n\n## WARNING: There is no .gtf file with CDS information. SNPGenie terminated.\n".
		"\n";
	}
	return 	$cds_file_name;
}

#########################################################################################
# Get an array with the product's start and stop sites, found in the .gtf file.
# Cause the program to DIE if the coordinates do not yield a multiple of three (i.e., 
# there is an incomplete codon)
# Returns an array with:
#	returned[0] = starting site
#	returned[1] = stop site
#	IF MORE SEGMENTS:
#	returned[2] = starting site 2
#	returned[3] = stop site 2
#	etc.
sub get_product_coordinates {
	my $product = $_[0];
	
	#print "\nThis time called for product: $product\n";
	
	my %start_site_h; # {segment #}->{start coordinate}
	my %stop_site_h;
	my %segment_lengths_h;
	
	open (CURRINFILE, $cds_file);
	while (<CURRINFILE>) { # go through the GTF file
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand #COMEBACK to \d+]
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			my $this_line = $_;
			my $product_present = 0;
			my $this_start;
			my $this_stop;
			
			#print "\nPRODUCT: $product\n";
			
			# incomplete segments
			if ($this_line =~/\s*gene_id\s+"$product\"/) { 
				if ($this_line =~/CDS\t(\d+)\t(\d+)/) {
					$product_present = 1;
					$this_start = $1;
					$this_stop = $2;
					
					#print "\n## $product\: $this_start\-$this_stop\n";
				}
			} elsif($this_line =~/\s*gene_id\s+"gene\:$product\"/) {
				if ($this_line =~/CDS\t(\d+)\t(\d+)/) {
					$product_present = 1;
					$this_start = $1;
					$this_stop = $2;
					#print "\n## $product\: $this_start\-$this_stop\n";
				}
			}
			
			# ADD SOME ERROR IF THE GENE NAME WASN'T RECOGNIZED BUT OUGHT TO HAVE BEEN?
			
			if($product_present == 1) {
				my $curr_max_key = max(keys %start_site_h); # this is the number of segments so far
				#print "curr_max_key is: $curr_max_key\n";
				my $curr_key = ($curr_max_key + 1);
				#print "curr_key is: $curr_key\n";
				$start_site_h{$curr_key} = $this_start;
				$stop_site_h{$curr_key} = $this_stop;
			}
		}
	}
	close CURRINFILE;
	
	#new incomplete segments
	my $num_segments = max(keys %start_site_h); # this is the number of segments total
	#print "\nAfter while, $product num segments: $num_segments\n";
	my $sum_of_lengths;
	
	for(my $i=1; $i<=$num_segments; $i++) {
		$segment_lengths_h{$i} = ($stop_site_h{$i} - $start_site_h{$i} + 1);
		$sum_of_lengths += $segment_lengths_h{$i};
	}
	
	# Make sure the sum of all nucleotides for this coding product add to a multiple of 3
	# incomplete segments
	if (($sum_of_lengths % 3) != 0) {
		chdir('SNPGenie_Results');
		open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
		# FILE | PRODUCT | SITE | CODON | WARNING
		print ERROR_FILE "N/A\t$product\tN/A\t".
			"Coordinates do not specify complete codon set (nucleotides a multiple of 3). SNPGenie terminated.\n";
		close ERROR_FILE;
		chdir('..');
		
		#unlink $curr_snp_report_name;
		
		die "\n\n## WARNING: The CDS coordinates for gene $product in the gtf file ".
			"do not yield a set of complete codons,\n".
			"## or are absent from the file. The number of nucleotides must ".
			"be a multiple of 3.\n## SNPGenie terminated.\n\n";
	}
	
	my @coord_arr;
	
	# Sort the segments
	my @sorted_segment_numbers;
	
	@sorted_segment_numbers = sort {($start_site_h{$a} <=> $start_site_h{$b}) || ($a <=> $b) } keys %start_site_h;
	
	foreach(@sorted_segment_numbers) {
		push(@coord_arr,$start_site_h{$_});
		push(@coord_arr,$stop_site_h{$_});
	}
	
#	for(my $i=1; $i<=$num_segments; $i++) {
#		push(@coord_arr,$start_site_h{$i});
#		push(@coord_arr,$stop_site_h{$i});
#	}
	
	return @coord_arr;
}

#########################################################################################
# Get the amino acid (single-letter code) encoded by a given DNA or RNA codon
# Returns an array with:
#	returned[0] = number of nonsynonymous sites
#	returned[1] = number of synonymous sites
sub get_amino_acid {
	#my ($codon) = @_;
	my ($codon) = uc($_[0]); # uc returns uppercase
	$codon =~ tr/U/T/;
	
	my $amino_acid;
	
	# Establish genetic code for use with synonymous sites; DNA or RNA
	my %code = (
		"AAA"=>"K","AAC"=>"N","AAG"=>"K","AAT"=>"N","ACA"=>"T","ACC"=>"T","ACG"=>"T",
		"ACT"=>"T","AGA"=>"R","AGC"=>"S","AGG"=>"R","AGT"=>"S","ATA"=>"I","ATC"=>"I",
		"ATG"=>"M","ATT"=>"I","CAA"=>"Q","CAC"=>"H","CAG"=>"Q","CAT"=>"H","CCA"=>"P",
		"CCC"=>"P","CCG"=>"P","CCT"=>"P","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R",
		"CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","GAA"=>"E","GAC"=>"D","GAG"=>"E",
		"GAT"=>"D","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GGA"=>"G","GGC"=>"G",
		"GGG"=>"G","GGT"=>"G","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","TAA"=>"*",
		"TAC"=>"Y","TAG"=>"*","TAT"=>"Y","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S",
		"TGA"=>"*","TGC"=>"C","TGG"=>"W","TGT"=>"C","TTA"=>"L","TTC"=>"F","TTG"=>"L",
		"TTT"=>"F"
	);
	
	$amino_acid = $code{$codon};
}


#########################################################################################
# Get the number of synonymous and nonsynonymous sites in a particular codon position
sub get_number_of_sites {
	my ($codon,$position) = @_;
	
	$codon = uc($codon);
	$codon =~ tr/U/T/;
	
	my $codon_N_position = "$codon"."_N_"."$position";
	my $codon_S_position = "$codon"."_S_"."$position";
	my $num_N;
	my $num_S;
	my @num_sites;
	
	# Numbers of nonsynonymous and synonymous sites for all three sites of every codon.
	# Per Nei-Gojobori, STOP codons have 0 sites. Results in arrays codon name in the format:
	# SITE1N SITE1S SITE2N SITE2S SITE3N SITE3S
	my %codon_type_position_hash = (
		"AAA_N_1" => 1,"AAA_S_1" => 0,"AAA_N_2" => 1,"AAA_S_2" => 0,"AAA_N_3" => 2/3,"AAA_S_3" => 1/3,
		"AAC_N_1" => 1,"AAC_S_1" => 0,"AAC_N_2" => 1,"AAC_S_2" => 0,"AAC_N_3" => 2/3,"AAC_S_3" => 1/3,
		"AAG_N_1" => 1,"AAG_S_1" => 0,"AAG_N_2" => 1,"AAG_S_2" => 0,"AAG_N_3" => 2/3,"AAG_S_3" => 1/3,
		"AAT_N_1" => 1,"AAT_S_1" => 0,"AAT_N_2" => 1,"AAT_S_2" => 0,"AAT_N_3" => 2/3,"AAT_S_3" => 1/3,
		"ACA_N_1" => 1,"ACA_S_1" => 0,"ACA_N_2" => 1,"ACA_S_2" => 0,"ACA_N_3" => 0,"ACA_S_3" => 1,
		"ACC_N_1" => 1,"ACC_S_1" => 0,"ACC_N_2" => 1,"ACC_S_2" => 0,"ACC_N_3" => 0,"ACC_S_3" => 1,
		"ACG_N_1" => 1,"ACG_S_1" => 0,"ACG_N_2" => 1,"ACG_S_2" => 0,"ACG_N_3" => 0,"ACG_S_3" => 1,
		"ACT_N_1" => 1,"ACT_S_1" => 0,"ACT_N_2" => 1,"ACT_S_2" => 0,"ACT_N_3" => 0,"ACT_S_3" => 1,
		"AGA_N_1" => 1/2,"AGA_S_1" => 1/2,"AGA_N_2" => 1,"AGA_S_2" => 0,"AGA_N_3" => 2/3,"AGA_S_3" => 1/3,
		"AGC_N_1" => 1,"AGC_S_1" => 0,"AGC_N_2" => 1,"AGC_S_2" => 0,"AGC_N_3" => 2/3,"AGC_S_3" => 1/3,
		"AGG_N_1" => 2/3,"AGG_S_1" => 1/3,"AGG_N_2" => 1,"AGG_S_2" => 0,"AGG_N_3" => 2/3,"AGG_S_3" => 1/3,
		"AGT_N_1" => 1,"AGT_S_1" => 0,"AGT_N_2" => 1,"AGT_S_2" => 0,"AGT_N_3" => 2/3,"AGT_S_3" => 1/3,
		"ATA_N_1" => 1,"ATA_S_1" => 0,"ATA_N_2" => 1,"ATA_S_2" => 0,"ATA_N_3" => 1/3,"ATA_S_3" => 2/3,
		"ATC_N_1" => 1,"ATC_S_1" => 0,"ATC_N_2" => 1,"ATC_S_2" => 0,"ATC_N_3" => 1/3,"ATC_S_3" => 2/3,
		"ATG_N_1" => 1,"ATG_S_1" => 0,"ATG_N_2" => 1,"ATG_S_2" => 0,"ATG_N_3" => 1,"ATG_S_3" => 0,
		"ATT_N_1" => 1,"ATT_S_1" => 0,"ATT_N_2" => 1,"ATT_S_2" => 0,"ATT_N_3" => 1/3,"ATT_S_3" => 2/3,
	
		"CAA_N_1" => 1,"CAA_S_1" => 0,"CAA_N_2" => 1,"CAA_S_2" => 0,"CAA_N_3" => 2/3,"CAA_S_3" => 1/3,
		"CAC_N_1" => 1,"CAC_S_1" => 0,"CAC_N_2" => 1,"CAC_S_2" => 0,"CAC_N_3" => 2/3,"CAC_S_3" => 1/3,
		"CAG_N_1" => 1,"CAG_S_1" => 0,"CAG_N_2" => 1,"CAG_S_2" => 0,"CAG_N_3" => 2/3,"CAG_S_3" => 1/3,
		"CAT_N_1" => 1,"CAT_S_1" => 0,"CAT_N_2" => 1,"CAT_S_2" => 0,"CAT_N_3" => 2/3,"CAT_S_3" => 1/3,
		"CCA_N_1" => 1,"CCA_S_1" => 0,"CCA_N_2" => 1,"CCA_S_2" => 0,"CCA_N_3" => 0,"CCA_S_3" => 1,
		"CCC_N_1" => 1,"CCC_S_1" => 0,"CCC_N_2" => 1,"CCC_S_2" => 0,"CCC_N_3" => 0,"CCC_S_3" => 1,
		"CCG_N_1" => 1,"CCG_S_1" => 0,"CCG_N_2" => 1,"CCG_S_2" => 0,"CCG_N_3" => 0,"CCG_S_3" => 1,
		"CCT_N_1" => 1,"CCT_S_1" => 0,"CCT_N_2" => 1,"CCT_S_2" => 0,"CCT_N_3" => 0,"CCT_S_3" => 1,
		"CGA_N_1" => 1/2,"CGA_S_1" => 1/2,"CGA_N_2" => 1,"CGA_S_2" => 0,"CGA_N_3" => 0,"CGA_S_3" => 1,
		"CGC_N_1" => 1,"CGC_S_1" => 0,"CGC_N_2" => 1,"CGC_S_2" => 0,"CGC_N_3" => 0,"CGC_S_3" => 1,
		"CGG_N_1" => 2/3,"CGG_S_1" => 1/3,"CGG_N_2" => 1,"CGG_S_2" => 0,"CGG_N_3" => 0,"CGG_S_3" => 1,
		"CGT_N_1" => 1,"CGT_S_1" => 0,"CGT_N_2" => 1,"CGT_S_2" => 0,"CGT_N_3" => 0,"CGT_S_3" => 1,
		"CTA_N_1" => 2/3,"CTA_S_1" => 1/3,"CTA_N_2" => 1,"CTA_S_2" => 0,"CTA_N_3" => 0,"CTA_S_3" => 1,
		"CTC_N_1" => 1,"CTC_S_1" => 0,"CTC_N_2" => 1,"CTC_S_2" => 0,"CTC_N_3" => 0,"CTC_S_3" => 1,
		"CTG_N_1" => 2/3,"CTG_S_1" => 1/3,"CTG_N_2" => 1,"CTG_S_2" => 0,"CTG_N_3" => 0,"CTG_S_3" => 1,
		"CTT_N_1" => 1,"CTT_S_1" => 0,"CTT_N_2" => 1,"CTT_S_2" => 0,"CTT_N_3" => 0,"CTT_S_3" => 1,
	
		"GAA_N_1" => 1,"GAA_S_1" => 0,"GAA_N_2" => 1,"GAA_S_2" => 0,"GAA_N_3" => 2/3,
		"GAA_S_3" => 1/3,"GAC_N_1" => 1,"GAC_S_1" => 0,"GAC_N_2" => 1,"GAC_S_2" => 0,
		"GAC_N_3" => 2/3,"GAC_S_3" => 1/3,"GAG_N_1" => 1,"GAG_S_1" => 0,"GAG_N_2" => 1,
		"GAG_S_2" => 0,"GAG_N_3" => 2/3,"GAG_S_3" => 1/3,"GAT_N_1" => 1,"GAT_S_1" => 0,
		"GAT_N_2" => 1,"GAT_S_2" => 0,"GAT_N_3" => 2/3,"GAT_S_3" => 1/3,"GCA_N_1" => 1,
		"GCA_S_1" => 0,"GCA_N_2" => 1,"GCA_S_2" => 0,"GCA_N_3" => 0,"GCA_S_3" => 1,
		"GCC_N_1" => 1,"GCC_S_1" => 0,"GCC_N_2" => 1,"GCC_S_2" => 0,"GCC_N_3" => 0,
		"GCC_S_3" => 1,"GCG_N_1" => 1,"GCG_S_1" => 0,"GCG_N_2" => 1,"GCG_S_2" => 0,
		"GCG_N_3" => 0,"GCG_S_3" => 1,"GCT_N_1" => 1,"GCT_S_1" => 0,"GCT_N_2" => 1,
		"GCT_S_2" => 0,"GCT_N_3" => 0,"GCT_S_3" => 1,"GGA_N_1" => 1,"GGA_S_1" => 0,
		"GGA_N_2" => 1,"GGA_S_2" => 0,"GGA_N_3" => 0,"GGA_S_3" => 1,"GGC_N_1" => 1,
		"GGC_S_1" => 0,"GGC_N_2" => 1,"GGC_S_2" => 0,"GGC_N_3" => 0,"GGC_S_3" => 1,
		"GGG_N_1" => 1,"GGG_S_1" => 0,"GGG_N_2" => 1,"GGG_S_2" => 0,"GGG_N_3" => 0,
		"GGG_S_3" => 1,"GGT_N_1" => 1,"GGT_S_1" => 0,"GGT_N_2" => 1,"GGT_S_2" => 0,
		"GGT_N_3" => 0,"GGT_S_3" => 1,"GTA_N_1" => 1,"GTA_S_1" => 0,"GTA_N_2" => 1,
		"GTA_S_2" => 0,"GTA_N_3" => 0,"GTA_S_3" => 1,"GTC_N_1" => 1,"GTC_S_1" => 0,
		"GTC_N_2" => 1,"GTC_S_2" => 0,"GTC_N_3" => 0,"GTC_S_3" => 1,"GTG_N_1" => 1,
		"GTG_S_1" => 0,"GTG_N_2" => 1,"GTG_S_2" => 0,"GTG_N_3" => 0,"GTG_S_3" => 1,
		"GTT_N_1" => 1,"GTT_S_1" => 0,"GTT_N_2" => 1,"GTT_S_2" => 0,"GTT_N_3" => 0,
		"GTT_S_3" => 1,
	
		"TAA_N_1" => 0,"TAA_S_1" => 0,"TAA_N_2" => 0,"TAA_S_2" => 0,"TAA_N_3" => 0,"TAA_S_3" => 0, # STOP
		"TAC_N_1" => 1,"TAC_S_1" => 0,"TAC_N_2" => 1,"TAC_S_2" => 0,"TAC_N_3" => 0,"TAC_S_3" => 1,
		"TAG_N_1" => 0,"TAG_S_1" => 0,"TAG_N_2" => 0,"TAG_S_2" => 0,"TAG_N_3" => 0,"TAG_S_3" => 0, # STOP
		"TAT_N_1" => 1,"TAT_S_1" => 0,"TAT_N_2" => 1,"TAT_S_2" => 0,"TAT_N_3" => 0,
		"TAT_S_3" => 1,"TCA_N_1" => 1,"TCA_S_1" => 0,"TCA_N_2" => 1,"TCA_S_2" => 0,
		"TCA_N_3" => 0,"TCA_S_3" => 1,"TCC_N_1" => 1,"TCC_S_1" => 0,"TCC_N_2" => 1,
		"TCC_S_2" => 0,"TCC_N_3" => 0,"TCC_S_3" => 1,"TCG_N_1" => 1,"TCG_S_1" => 0,
		"TCG_N_2" => 1,"TCG_S_2" => 0,"TCG_N_3" => 0,"TCG_S_3" => 1,"TCT_N_1" => 1,
		"TCT_S_1" => 0,"TCT_N_2" => 1,"TCT_S_2" => 0,"TCT_N_3" => 0,"TCT_S_3" => 1,
		"TGA_N_1" => 0,"TGA_S_1" => 0,"TGA_N_2" => 0,"TGA_S_2" => 0,"TGA_N_3" => 0,"TGA_S_3" => 0, # STOP
		"TGC_N_1" => 1,"TGC_S_1" => 0,"TGC_N_2" => 1,"TGC_S_2" => 0,"TGC_N_3" => 1/2,"TGC_S_3" => 1/2,
		"TGG_N_1" => 1,"TGG_S_1" => 0,"TGG_N_2" => 1,"TGG_S_2" => 0,"TGG_N_3" => 1,"TGG_S_3" => 0,
		"TGT_N_1" => 1,"TGT_S_1" => 0,"TGT_N_2" => 1,"TGT_S_2" => 0,"TGT_N_3" => 1/2,"TGT_S_3" => 1/2,
		"TTA_N_1" => 2/3,"TTA_S_1" => 1/3,"TTA_N_2" => 1,"TTA_S_2" => 0,"TTA_N_3" => 2/3,"TTA_S_3" => 1/3,
		"TTC_N_1" => 1,"TTC_S_1" => 0,"TTC_N_2" => 1,"TTC_S_2" => 0,"TTC_N_3" => 2/3,"TTC_S_3" => 1/3,
		"TTG_N_1" => 2/3,"TTG_S_1" => 1/3,"TTG_N_2" => 1,"TTG_S_2" => 0,"TTG_N_3" => 2/3,"TTG_S_3" => 1/3,
		"TTT_N_1" => 1,"TTT_S_1" => 0,"TTT_N_2" => 1,"TTT_S_2" => 0,"TTT_N_3" => 2/3,"TTT_S_3" => 1/3,
	);
	
	$num_N = $codon_type_position_hash{$codon_N_position};
	$num_S = $codon_type_position_hash{$codon_S_position};
	
	@num_sites = ($num_N,$num_S);
	
	return @num_sites;
}


#########################################################################################
# Returns the number of differences between two codons, averaged over all minimal-evolution paths
# Per Nei-Gojobori, we return 0 differences for STOP codons, because they have 0 sites.
sub return_avg_diffs {
	my ($codon1,$codon2) = @_;
	
	my %all_diffs_hh = (
		'AAA' => {
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 0, 'S' => 1},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 0},
			'ACC' => {'N' => 1.5, 'S' => 0.5},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1.5, 'S' => 0.5},
			'AGA' => {'N' => 1, 'S' => 0},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 1, 'S' => 1},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1.5, 'S' => 0.5},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1.5, 'S' => 0.5},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 1, 'S' => 1},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 0},
			'CCC' => {'N' => 2.5, 'S' => 0.5},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2.5, 'S' => 0.5},
			'CGA' => {'N' => 1.5, 'S' => 0.5},
			'CGC' => {'N' => 2.5, 'S' => 0.5},
			'CGG' => {'N' => 1.5, 'S' => 1.5},
			'CGT' => {'N' => 2.5, 'S' => 0.5},
			'CTA' => {'N' => 2, 'S' => 0},
			'CTC' => {'N' => 2.5, 'S' => 0.5},
			'CTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTT' => {'N' => 2.5, 'S' => 0.5},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 1, 'S' => 1},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2.5, 'S' => 0.5},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2.5, 'S' => 0.5},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2.5, 'S' => 0.5},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 3, 'S' => 0},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 3, 'S' => 0},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 2.75, 'S' => 0.25},
			'TTG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTT' => {'N' => 2.75, 'S' => 0.25},
		},
		'AAC' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 0, 'S' => 1},
			'ACA' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 1, 'S' => 0},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 1},
			'CCA' => {'N' => 2.5, 'S' => 0.5},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.5, 'S' => 0.5},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2.5, 'S' => 0.5},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 1},
			'GCA' => {'N' => 2.5, 'S' => 0.5},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.5, 'S' => 0.5},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2.5, 'S' => 0.5},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.25, 'S' => 0.75},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 3, 'S' => 0},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2.75, 'S' => 0.25},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 3, 'S' => 0},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'AAG' => {
			'AAA' => {'N' => 0, 'S' => 1},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1.5, 'S' => 0.5},
			'ACG' => {'N' => 1, 'S' => 0},
			'ACT' => {'N' => 1.5, 'S' => 0.5},
			'AGA' => {'N' => 1, 'S' => 1},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 1, 'S' => 1},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2.5, 'S' => 0.5},
			'CCG' => {'N' => 2, 'S' => 0},
			'CCT' => {'N' => 2.5, 'S' => 0.5},
			'CGA' => {'N' => 1.5, 'S' => 1.5},
			'CGC' => {'N' => 2.5, 'S' => 0.5},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 2.5, 'S' => 0.5},
			'CTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CTG' => {'N' => 2, 'S' => 0},
			'CTT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GAA' => {'N' => 1, 'S' => 1},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2.5, 'S' => 0.5},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2.5, 'S' => 0.5},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 3, 'S' => 0},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 3, 'S' => 0},
			'TTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTC' => {'N' => 3, 'S' => 0},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 3, 'S' => 0},
		},
		'AAT' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAC' => {'N' => 0, 'S' => 1},
			'AAG' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 0},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 1, 'S' => 1},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 2.5, 'S' => 0.5},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.5, 'S' => 0.5},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 2.5, 'S' => 0.5},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 1, 'S' => 1},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 2.5, 'S' => 0.5},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.5, 'S' => 0.5},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 2.5, 'S' => 0.5},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 0},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.25, 'S' => 0.75},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 3, 'S' => 0},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2.75, 'S' => 0.25},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 3, 'S' => 0},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'ACA' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAC' => {'N' => 1.5, 'S' => 0.5},
			'AAG' => {'N' => 1, 'S' => 1},
			'AAT' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 0, 'S' => 1},
			'ACG' => {'N' => 0, 'S' => 1},
			'ACT' => {'N' => 0, 'S' => 1},
			'AGA' => {'N' => 1, 'S' => 0},
			'AGC' => {'N' => 1.5, 'S' => 0.5},
			'AGG' => {'N' => 1, 'S' => 1},
			'AGT' => {'N' => 1.5, 'S' => 0.5},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 1},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 0},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1.5, 'S' => 0.5},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 1.5, 'S' => 1.5},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 0},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 0},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1, 'S' => 0},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.25, 'S' => 0.75},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 2.25, 'S' => 0.75},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
		},
		'ACC' => {
			'AAA' => {'N' => 1.5, 'S' => 0.5},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 1.5, 'S' => 0.5},
			'AAT' => {'N' => 1, 'S' => 1},
			'ACA' => {'N' => 0, 'S' => 1},
			'ACG' => {'N' => 0, 'S' => 1},
			'ACT' => {'N' => 0, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 1, 'S' => 1},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 0},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 0},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 0},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.5, 'S' => 0.5},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2.5, 'S' => 0.5},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'ACG' => {
			'AAA' => {'N' => 1, 'S' => 1},
			'AAC' => {'N' => 1.5, 'S' => 0.5},
			'AAG' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 1.5, 'S' => 0.5},
			'ACA' => {'N' => 0, 'S' => 1},
			'ACC' => {'N' => 0, 'S' => 1},
			'ACT' => {'N' => 0, 'S' => 1},
			'AGA' => {'N' => 1, 'S' => 1},
			'AGC' => {'N' => 1.5, 'S' => 0.5},
			'AGG' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 1.5, 'S' => 0.5},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1.5, 'S' => 0.5},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 1.5, 'S' => 0.5},
			'CAA' => {'N' => 2, 'S' => 1},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 0},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1.5, 'S' => 1.5},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTG' => {'N' => 2, 'S' => 0},
			'CTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 0},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 0},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.5, 'S' => 0.5},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 2.5, 'S' => 0.5},
			'TTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTC' => {'N' => 2.5, 'S' => 0.5},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 2.5, 'S' => 0.5},
		},
		'ACT' => {
			'AAA' => {'N' => 1.5, 'S' => 0.5},
			'AAC' => {'N' => 1, 'S' => 1},
			'AAG' => {'N' => 1.5, 'S' => 0.5},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 0, 'S' => 1},
			'ACC' => {'N' => 0, 'S' => 1},
			'ACG' => {'N' => 0, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 1, 'S' => 1},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 0},
			'CGA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 0},
			'GGA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.5, 'S' => 0.5},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2.5, 'S' => 0.5},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'AGA' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 1},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 0},
			'ACC' => {'N' => 1.5, 'S' => 0.5},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 0, 'S' => 1},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1.5, 'S' => 0.5},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1.5, 'S' => 0.5},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 1.5, 'S' => 1.5},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 1.5, 'S' => 0.5},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 1.5, 'S' => 1.5},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 0, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 0, 'S' => 2},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 1.5, 'S' => 0.5},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 1.66666666666667, 'S' => 1.33333333333333},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGA' => {'N' => 1, 'S' => 0},
			'GGC' => {'N' => 1.5, 'S' => 0.5},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1.5, 'S' => 0.5},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 3, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 3, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 1, 'S' => 1},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 2.75, 'S' => 0.25},
			'TTG' => {'N' => 2.25, 'S' => 0.75},
			'TTT' => {'N' => 2.75, 'S' => 0.25},
		},
		'AGC' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 1},
			'ACA' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 1, 'S' => 0},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 0, 'S' => 1},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 0},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 1.5, 'S' => 0.5},
			'GGC' => {'N' => 1, 'S' => 0},
			'GGG' => {'N' => 1.5, 'S' => 0.5},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.5, 'S' => 0.5},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 2.75, 'S' => 0.25},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 3, 'S' => 0},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'AGG' => {
			'AAA' => {'N' => 1, 'S' => 1},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1.5, 'S' => 0.5},
			'ACG' => {'N' => 1, 'S' => 0},
			'ACT' => {'N' => 1.5, 'S' => 0.5},
			'AGA' => {'N' => 0, 'S' => 1},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 1.5, 'S' => 1.5},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 1.5, 'S' => 1.5},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 1.5, 'S' => 0.5},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 0, 'S' => 2},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 0, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 1.66666666666667, 'S' => 1.33333333333333},
			'CTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1.5, 'S' => 0.5},
			'GGG' => {'N' => 1, 'S' => 0},
			'GGT' => {'N' => 1.5, 'S' => 0.5},
			'GTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTC' => {'N' => 2.5, 'S' => 0.5},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 3, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 3, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2.25, 'S' => 0.75},
			'TTC' => {'N' => 3, 'S' => 0},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 3, 'S' => 0},
		},
		'AGT' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 1},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 0},
			'AGA' => {'N' => 1, 'S' => 0},
			'AGC' => {'N' => 0, 'S' => 1},
			'AGG' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 0},
			'CTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 1.5, 'S' => 0.5},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1.5, 'S' => 0.5},
			'GGT' => {'N' => 1, 'S' => 0},
			'GTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.5, 'S' => 0.5},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 2.75, 'S' => 0.25},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 3, 'S' => 0},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'ATA' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAC' => {'N' => 1.5, 'S' => 0.5},
			'AAG' => {'N' => 1.5, 'S' => 0.5},
			'AAT' => {'N' => 1.5, 'S' => 0.5},
			'ACA' => {'N' => 1, 'S' => 0},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 1, 'S' => 0},
			'AGC' => {'N' => 1.5, 'S' => 0.5},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 0, 'S' => 1},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 0, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 2, 'S' => 0},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1.5, 'S' => 0.5},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 0},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTA' => {'N' => 1, 'S' => 0},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.5, 'S' => 0.5},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.5, 'S' => 0.5},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TCT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.5, 'S' => 0.5},
			'TGG' => {'N' => 2.5, 'S' => 0.5},
			'TGT' => {'N' => 2.5, 'S' => 0.5},
			'TTA' => {'N' => 1, 'S' => 0},
			'TTC' => {'N' => 1.5, 'S' => 0.5},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 1.5, 'S' => 0.5},
		},
		'ATC' => {
			'AAA' => {'N' => 1.5, 'S' => 0.5},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 1},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 0},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 0, 'S' => 1},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 0, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 0},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 0},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.5, 'S' => 0.5},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 3, 'S' => 0},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'ATG' => {
			'AAA' => {'N' => 1.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 1.5, 'S' => 0.5},
			'ACG' => {'N' => 1, 'S' => 0},
			'ACT' => {'N' => 1.5, 'S' => 0.5},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CAC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'CCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCG' => {'N' => 2, 'S' => 0},
			'CCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CTA' => {'N' => 1.5, 'S' => 0.5},
			'CTC' => {'N' => 1.5, 'S' => 0.5},
			'CTG' => {'N' => 1, 'S' => 0},
			'CTT' => {'N' => 1.5, 'S' => 0.5},
			'GAA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GAC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGC' => {'N' => 2.5, 'S' => 0.5},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2.5, 'S' => 0.5},
			'GTA' => {'N' => 1.5, 'S' => 0.5},
			'GTC' => {'N' => 1.5, 'S' => 0.5},
			'GTG' => {'N' => 1, 'S' => 0},
			'GTT' => {'N' => 1.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 3, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 3, 'S' => 0},
			'TCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 3, 'S' => 0},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 3, 'S' => 0},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 1, 'S' => 0},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'ATT' => {
			'AAA' => {'N' => 1.5, 'S' => 0.5},
			'AAC' => {'N' => 1, 'S' => 1},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1, 'S' => 0},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 0, 'S' => 1},
			'ATC' => {'N' => 0, 'S' => 1},
			'ATG' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 1, 'S' => 0},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.5, 'S' => 0.5},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 3, 'S' => 0},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'CAA' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 1},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 0},
			'ACC' => {'N' => 2.5, 'S' => 0.5},
			'ACG' => {'N' => 2, 'S' => 1},
			'ACT' => {'N' => 2.5, 'S' => 0.5},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AGG' => {'N' => 1.5, 'S' => 1.5},
			'AGT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'ATA' => {'N' => 2, 'S' => 0},
			'ATC' => {'N' => 2.5, 'S' => 0.5},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 0, 'S' => 1},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 0},
			'CCC' => {'N' => 1.5, 'S' => 0.5},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1.5, 'S' => 0.5},
			'CGA' => {'N' => 1, 'S' => 0},
			'CGC' => {'N' => 1.5, 'S' => 0.5},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1.5, 'S' => 0.5},
			'CTA' => {'N' => 1, 'S' => 0},
			'CTC' => {'N' => 1.5, 'S' => 0.5},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1.5, 'S' => 0.5},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 1, 'S' => 1},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2.5, 'S' => 0.5},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2.5, 'S' => 0.5},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2.5, 'S' => 0.5},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2.5, 'S' => 0.5},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2.5, 'S' => 0.5},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 2.5, 'S' => 0.5},
			'TTG' => {'N' => 1, 'S' => 2},
			'TTT' => {'N' => 2.5, 'S' => 0.5},
		},
		'CAC' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 1},
			'ACA' => {'N' => 2.5, 'S' => 0.5},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2.5, 'S' => 0.5},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2.5, 'S' => 0.5},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 0, 'S' => 1},
			'CCA' => {'N' => 1.5, 'S' => 0.5},
			'CCC' => {'N' => 1, 'S' => 0},
			'CCG' => {'N' => 1.5, 'S' => 0.5},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1.5, 'S' => 0.5},
			'CGC' => {'N' => 1, 'S' => 0},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 1.5, 'S' => 0.5},
			'CTC' => {'N' => 1, 'S' => 0},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 1},
			'GCA' => {'N' => 2.5, 'S' => 0.5},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.5, 'S' => 0.5},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2.5, 'S' => 0.5},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2.5, 'S' => 0.5},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.25, 'S' => 0.75},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.75, 'S' => 0.25},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2.25, 'S' => 0.75},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2.25, 'S' => 0.75},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'CAG' => {
			'AAA' => {'N' => 1, 'S' => 1},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2.5, 'S' => 0.5},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.5, 'S' => 0.5},
			'AGA' => {'N' => 1.5, 'S' => 1.5},
			'AGC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'ATA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'CAA' => {'N' => 0, 'S' => 1},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1.5, 'S' => 0.5},
			'CCG' => {'N' => 1, 'S' => 0},
			'CCT' => {'N' => 1.5, 'S' => 0.5},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1.5, 'S' => 0.5},
			'CGG' => {'N' => 1, 'S' => 0},
			'CGT' => {'N' => 1.5, 'S' => 0.5},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1.5, 'S' => 0.5},
			'CTG' => {'N' => 1, 'S' => 0},
			'CTT' => {'N' => 1.5, 'S' => 0.5},
			'GAA' => {'N' => 1, 'S' => 1},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2.5, 'S' => 0.5},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2.5, 'S' => 0.5},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2.5, 'S' => 0.5},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2.5, 'S' => 0.5},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2.5, 'S' => 0.5},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.75, 'S' => 0.25},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 2.75, 'S' => 0.25},
			'TTA' => {'N' => 1, 'S' => 2},
			'TTC' => {'N' => 2.5, 'S' => 0.5},
			'TTG' => {'N' => 1, 'S' => 1},
			'TTT' => {'N' => 2.5, 'S' => 0.5},
		},
		'CAT' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 1},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 2.5, 'S' => 0.5},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 2.5, 'S' => 0.5},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2.5, 'S' => 0.5},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAC' => {'N' => 0, 'S' => 1},
			'CAG' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 1.5, 'S' => 0.5},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1.5, 'S' => 0.5},
			'CCT' => {'N' => 1, 'S' => 0},
			'CGA' => {'N' => 1.5, 'S' => 0.5},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 1, 'S' => 0},
			'CTA' => {'N' => 1.5, 'S' => 0.5},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 1, 'S' => 0},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 1, 'S' => 1},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 2.5, 'S' => 0.5},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.5, 'S' => 0.5},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 2.5, 'S' => 0.5},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 2.5, 'S' => 0.5},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 0},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.25, 'S' => 0.75},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.75, 'S' => 0.25},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2.25, 'S' => 0.75},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2.25, 'S' => 0.75},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'CCA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.5, 'S' => 0.5},
			'AAG' => {'N' => 2, 'S' => 1},
			'AAT' => {'N' => 2.5, 'S' => 0.5},
			'ACA' => {'N' => 1, 'S' => 0},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGG' => {'N' => 1.5, 'S' => 1.5},
			'AGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATA' => {'N' => 2, 'S' => 0},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAC' => {'N' => 1.5, 'S' => 0.5},
			'CAG' => {'N' => 1, 'S' => 1},
			'CAT' => {'N' => 1.5, 'S' => 0.5},
			'CCC' => {'N' => 0, 'S' => 1},
			'CCG' => {'N' => 0, 'S' => 1},
			'CCT' => {'N' => 0, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 0},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 0},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 0},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1, 'S' => 0},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTG' => {'N' => 1.5, 'S' => 1.5},
			'TTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
		},
		'CCC' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 0},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2, 'S' => 1},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 2, 'S' => 1},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 1, 'S' => 1},
			'CCA' => {'N' => 0, 'S' => 1},
			'CCG' => {'N' => 0, 'S' => 1},
			'CCT' => {'N' => 0, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 0},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 0},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 0},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 0},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2, 'S' => 1},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'CCG' => {
			'AAA' => {'N' => 2, 'S' => 1},
			'AAC' => {'N' => 2.5, 'S' => 0.5},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 2.5, 'S' => 0.5},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 0},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 1.5},
			'AGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CAA' => {'N' => 1, 'S' => 1},
			'CAC' => {'N' => 1.5, 'S' => 0.5},
			'CAG' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 1.5, 'S' => 0.5},
			'CCA' => {'N' => 0, 'S' => 1},
			'CCC' => {'N' => 0, 'S' => 1},
			'CCT' => {'N' => 0, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 0},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 0},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 0},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 0},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTA' => {'N' => 1.5, 'S' => 1.5},
			'TTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
		},
		'CCT' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 0},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2, 'S' => 1},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2, 'S' => 1},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 1},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 0, 'S' => 1},
			'CCC' => {'N' => 0, 'S' => 1},
			'CCG' => {'N' => 0, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 0},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 0},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 0},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2, 'S' => 1},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'CGA' => {
			'AAA' => {'N' => 1.5, 'S' => 0.5},
			'AAC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AAG' => {'N' => 1.5, 'S' => 1.5},
			'AAT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ACA' => {'N' => 1.5, 'S' => 0.5},
			'ACC' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'ACG' => {'N' => 1.5, 'S' => 1.5},
			'ACT' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'AGA' => {'N' => 0, 'S' => 1},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 0, 'S' => 2},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'ATG' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'ATT' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAC' => {'N' => 1.5, 'S' => 0.5},
			'CAG' => {'N' => 1, 'S' => 1},
			'CAT' => {'N' => 1.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 0},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 0, 'S' => 1},
			'CGG' => {'N' => 0, 'S' => 1},
			'CGT' => {'N' => 0, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 0},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 0},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 1, 'S' => 1},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 1.25, 'S' => 1.75},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'CGC' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2, 'S' => 1},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 1, 'S' => 1},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 1, 'S' => 1},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 2, 'S' => 1},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 1, 'S' => 1},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 0},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 0, 'S' => 1},
			'CGG' => {'N' => 0, 'S' => 1},
			'CGT' => {'N' => 0, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 0},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 0},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'CGG' => {
			'AAA' => {'N' => 1.5, 'S' => 1.5},
			'AAC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AAG' => {'N' => 1.5, 'S' => 0.5},
			'AAT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ACA' => {'N' => 1.5, 'S' => 1.5},
			'ACC' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'ACG' => {'N' => 1.5, 'S' => 0.5},
			'ACT' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'AGA' => {'N' => 0, 'S' => 2},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 0, 'S' => 1},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'ATC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CAA' => {'N' => 1, 'S' => 1},
			'CAC' => {'N' => 1.5, 'S' => 0.5},
			'CAG' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 1.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 0},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 0, 'S' => 1},
			'CGC' => {'N' => 0, 'S' => 1},
			'CGT' => {'N' => 0, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 0},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 0},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.5, 'S' => 0.5},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.5, 'S' => 0.5},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1.5, 'S' => 0.5},
			'TGG' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 1.5, 'S' => 0.5},
			'TTA' => {'N' => 1.25, 'S' => 1.75},
			'TTC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
		},
		'CGT' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2, 'S' => 1},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 1, 'S' => 1},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 1, 'S' => 1},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 2, 'S' => 1},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 1},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 0},
			'CGA' => {'N' => 0, 'S' => 1},
			'CGC' => {'N' => 0, 'S' => 1},
			'CGG' => {'N' => 0, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 0},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 0},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'CTA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.5, 'S' => 0.5},
			'AAG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AAT' => {'N' => 2.5, 'S' => 0.5},
			'ACA' => {'N' => 2, 'S' => 0},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGG' => {'N' => 1.66666666666667, 'S' => 1.33333333333333},
			'AGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAC' => {'N' => 1.5, 'S' => 0.5},
			'CAG' => {'N' => 1, 'S' => 1},
			'CAT' => {'N' => 1.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 0},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 0},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 0, 'S' => 1},
			'CTG' => {'N' => 0, 'S' => 1},
			'CTT' => {'N' => 0, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 0},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1.5, 'S' => 0.5},
			'TCC' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'TCG' => {'N' => 1.5, 'S' => 1.5},
			'TCT' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 1.5, 'S' => 1.5},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 0, 'S' => 1},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 0, 'S' => 2},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'CTC' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 1, 'S' => 1},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 1, 'S' => 1},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 0},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 0},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 0, 'S' => 1},
			'CTG' => {'N' => 0, 'S' => 1},
			'CTT' => {'N' => 0, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 0},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 1, 'S' => 1},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'CTG' => {
			'AAA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AAC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'ACA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGA' => {'N' => 1.66666666666667, 'S' => 1.33333333333333},
			'AGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1.5, 'S' => 0.5},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 1.5, 'S' => 0.5},
			'CAA' => {'N' => 1, 'S' => 1},
			'CAC' => {'N' => 1.5, 'S' => 0.5},
			'CAG' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 1.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 0},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 0},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 0, 'S' => 1},
			'CTC' => {'N' => 0, 'S' => 1},
			'CTT' => {'N' => 0, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.5, 'S' => 0.5},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.5, 'S' => 0.5},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 0},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1.5, 'S' => 1.5},
			'TCC' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'TCG' => {'N' => 1.5, 'S' => 0.5},
			'TCT' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TTA' => {'N' => 0, 'S' => 2},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 0, 'S' => 1},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'CTT' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 1, 'S' => 1},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 1.5, 'S' => 0.5},
			'CAC' => {'N' => 1, 'S' => 1},
			'CAG' => {'N' => 1.5, 'S' => 0.5},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 0},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 0},
			'CTA' => {'N' => 0, 'S' => 1},
			'CTC' => {'N' => 0, 'S' => 1},
			'CTG' => {'N' => 0, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 1, 'S' => 1},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'GAA' => {
			'AAA' => {'N' => 1, 'S' => 0},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 1},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 0},
			'ACC' => {'N' => 2.5, 'S' => 0.5},
			'ACG' => {'N' => 2, 'S' => 1},
			'ACT' => {'N' => 2.5, 'S' => 0.5},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'AGG' => {'N' => 2, 'S' => 1},
			'AGT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATA' => {'N' => 2, 'S' => 0},
			'ATC' => {'N' => 2.5, 'S' => 0.5},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2.5, 'S' => 0.5},
			'CAA' => {'N' => 1, 'S' => 0},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 1, 'S' => 1},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 0},
			'CCC' => {'N' => 2.5, 'S' => 0.5},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2.5, 'S' => 0.5},
			'CGA' => {'N' => 2, 'S' => 0},
			'CGC' => {'N' => 2.5, 'S' => 0.5},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2.5, 'S' => 0.5},
			'CTA' => {'N' => 2, 'S' => 0},
			'CTC' => {'N' => 2.5, 'S' => 0.5},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 0, 'S' => 1},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 0},
			'GCC' => {'N' => 1.5, 'S' => 0.5},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1.5, 'S' => 0.5},
			'GGA' => {'N' => 1, 'S' => 0},
			'GGC' => {'N' => 1.5, 'S' => 0.5},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1.5, 'S' => 0.5},
			'GTA' => {'N' => 1, 'S' => 0},
			'GTC' => {'N' => 1.5, 'S' => 0.5},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 2.75, 'S' => 0.25},
			'TTG' => {'N' => 2, 'S' => 1},
			'TTT' => {'N' => 2.75, 'S' => 0.25},
		},
		'GAC' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 1},
			'ACA' => {'N' => 2.5, 'S' => 0.5},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 1},
			'CCA' => {'N' => 2.5, 'S' => 0.5},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.5, 'S' => 0.5},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2.5, 'S' => 0.5},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2.5, 'S' => 0.5},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2.5, 'S' => 0.5},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2.5, 'S' => 0.5},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 0, 'S' => 1},
			'GCA' => {'N' => 1.5, 'S' => 0.5},
			'GCC' => {'N' => 1, 'S' => 0},
			'GCG' => {'N' => 1.5, 'S' => 0.5},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 1.5, 'S' => 0.5},
			'GGC' => {'N' => 1, 'S' => 0},
			'GGG' => {'N' => 1.5, 'S' => 0.5},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 1.5, 'S' => 0.5},
			'GTC' => {'N' => 1, 'S' => 0},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.25, 'S' => 0.75},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.75, 'S' => 0.25},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2.75, 'S' => 0.25},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2.75, 'S' => 0.25},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'GAG' => {
			'AAA' => {'N' => 1, 'S' => 1},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 1, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2.5, 'S' => 0.5},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.5, 'S' => 0.5},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATC' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'CAA' => {'N' => 1, 'S' => 1},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 1, 'S' => 0},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2.5, 'S' => 0.5},
			'CCG' => {'N' => 2, 'S' => 0},
			'CCT' => {'N' => 2.5, 'S' => 0.5},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2.5, 'S' => 0.5},
			'CGG' => {'N' => 2, 'S' => 0},
			'CGT' => {'N' => 2.5, 'S' => 0.5},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2.5, 'S' => 0.5},
			'CTG' => {'N' => 2, 'S' => 0},
			'CTT' => {'N' => 2.5, 'S' => 0.5},
			'GAA' => {'N' => 0, 'S' => 1},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1.5, 'S' => 0.5},
			'GCG' => {'N' => 1, 'S' => 0},
			'GCT' => {'N' => 1.5, 'S' => 0.5},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1.5, 'S' => 0.5},
			'GGG' => {'N' => 1, 'S' => 0},
			'GGT' => {'N' => 1.5, 'S' => 0.5},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1.5, 'S' => 0.5},
			'GTG' => {'N' => 1, 'S' => 0},
			'GTT' => {'N' => 1.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.5, 'S' => 0.5},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.75, 'S' => 0.25},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 2.75, 'S' => 0.25},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2.75, 'S' => 0.25},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 2.75, 'S' => 0.25},
		},
		'GAT' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 1},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 2.5, 'S' => 0.5},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.83333333333333, 'S' => 0.166666666666667},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 1, 'S' => 1},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 2.5, 'S' => 0.5},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.5, 'S' => 0.5},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 2.5, 'S' => 0.5},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2.5, 'S' => 0.5},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 2.5, 'S' => 0.5},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.5, 'S' => 0.5},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAC' => {'N' => 0, 'S' => 1},
			'GAG' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 1.5, 'S' => 0.5},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1.5, 'S' => 0.5},
			'GCT' => {'N' => 1, 'S' => 0},
			'GGA' => {'N' => 1.5, 'S' => 0.5},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1.5, 'S' => 0.5},
			'GGT' => {'N' => 1, 'S' => 0},
			'GTA' => {'N' => 1.5, 'S' => 0.5},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 0},
			'TCA' => {'N' => 2.25, 'S' => 0.75},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.25, 'S' => 0.75},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.75, 'S' => 0.25},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2.75, 'S' => 0.25},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2.75, 'S' => 0.25},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'GCA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.5, 'S' => 0.5},
			'AAG' => {'N' => 2, 'S' => 1},
			'AAT' => {'N' => 2.5, 'S' => 0.5},
			'ACA' => {'N' => 1, 'S' => 0},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGG' => {'N' => 2, 'S' => 1},
			'AGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATA' => {'N' => 2, 'S' => 0},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 1},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 0},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 0},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 0},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAC' => {'N' => 1.5, 'S' => 0.5},
			'GAG' => {'N' => 1, 'S' => 1},
			'GAT' => {'N' => 1.5, 'S' => 0.5},
			'GCC' => {'N' => 0, 'S' => 1},
			'GCG' => {'N' => 0, 'S' => 1},
			'GCT' => {'N' => 0, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 0},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 0},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1, 'S' => 0},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTG' => {'N' => 2, 'S' => 1},
			'TTT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
		},
		'GCC' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 0},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 2, 'S' => 1},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 0},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 1.5, 'S' => 0.5},
			'GAT' => {'N' => 1, 'S' => 1},
			'GCA' => {'N' => 0, 'S' => 1},
			'GCG' => {'N' => 0, 'S' => 1},
			'GCT' => {'N' => 0, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 0},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 0},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 0},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'GCG' => {
			'AAA' => {'N' => 2, 'S' => 1},
			'AAC' => {'N' => 2.5, 'S' => 0.5},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 2.5, 'S' => 0.5},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 0},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CAA' => {'N' => 2, 'S' => 1},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 0},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 0},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2, 'S' => 0},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 1},
			'GAC' => {'N' => 1.5, 'S' => 0.5},
			'GAG' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 1.5, 'S' => 0.5},
			'GCA' => {'N' => 0, 'S' => 1},
			'GCC' => {'N' => 0, 'S' => 1},
			'GCT' => {'N' => 0, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 0},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 0},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.25, 'S' => 0.75},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.25, 'S' => 0.75},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 0},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
		},
		'GCT' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 0},
			'AGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2, 'S' => 1},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 0},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 1.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 1},
			'GAG' => {'N' => 1.5, 'S' => 0.5},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 0, 'S' => 1},
			'GCC' => {'N' => 0, 'S' => 1},
			'GCG' => {'N' => 0, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 0},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'GGA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAG' => {'N' => 2, 'S' => 1},
			'AAT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'ACA' => {'N' => 2, 'S' => 0},
			'ACC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACG' => {'N' => 2, 'S' => 1},
			'ACT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGA' => {'N' => 1, 'S' => 0},
			'AGC' => {'N' => 1.5, 'S' => 0.5},
			'AGG' => {'N' => 1, 'S' => 1},
			'AGT' => {'N' => 1.5, 'S' => 0.5},
			'ATA' => {'N' => 2, 'S' => 0},
			'ATC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 1},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 2, 'S' => 0},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 0},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 0},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAC' => {'N' => 1.5, 'S' => 0.5},
			'GAG' => {'N' => 1, 'S' => 1},
			'GAT' => {'N' => 1.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 0},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 0, 'S' => 1},
			'GGG' => {'N' => 0, 'S' => 1},
			'GGT' => {'N' => 0, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 0},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 1, 'S' => 1},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 2.25, 'S' => 0.75},
			'TTG' => {'N' => 2, 'S' => 1},
			'TTT' => {'N' => 2.25, 'S' => 0.75},
		},
		'GGC' => {
			'AAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.5, 'S' => 0.5},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 0},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 1.5, 'S' => 0.5},
			'GAT' => {'N' => 1, 'S' => 1},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 0},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 0, 'S' => 1},
			'GGG' => {'N' => 0, 'S' => 1},
			'GGT' => {'N' => 0, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 0},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 2.5, 'S' => 0.5},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 2.5, 'S' => 0.5},
			'TTT' => {'N' => 2, 'S' => 1},
		},
		'GGG' => {
			'AAA' => {'N' => 2, 'S' => 1},
			'AAC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGA' => {'N' => 1, 'S' => 1},
			'AGC' => {'N' => 1.5, 'S' => 0.5},
			'AGG' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 1.5, 'S' => 0.5},
			'ATA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATC' => {'N' => 2.5, 'S' => 0.5},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 2.5, 'S' => 0.5},
			'CAA' => {'N' => 2, 'S' => 1},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2, 'S' => 0},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 0},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2, 'S' => 0},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 1},
			'GAC' => {'N' => 1.5, 'S' => 0.5},
			'GAG' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 1.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 0},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 0, 'S' => 1},
			'GGC' => {'N' => 0, 'S' => 1},
			'GGT' => {'N' => 0, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 0},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.5, 'S' => 0.5},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.5, 'S' => 0.5},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1.5, 'S' => 0.5},
			'TGG' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 1.5, 'S' => 0.5},
			'TTA' => {'N' => 2, 'S' => 1},
			'TTC' => {'N' => 2.5, 'S' => 0.5},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 2.5, 'S' => 0.5},
		},
		'GGT' => {
			'AAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 1.5, 'S' => 0.5},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 1.5, 'S' => 0.5},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.5, 'S' => 0.5},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1, 'S' => 1},
			'CGT' => {'N' => 1, 'S' => 0},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2, 'S' => 1},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 1.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 1},
			'GAG' => {'N' => 1.5, 'S' => 0.5},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 0},
			'GGA' => {'N' => 0, 'S' => 1},
			'GGC' => {'N' => 0, 'S' => 1},
			'GGG' => {'N' => 0, 'S' => 1},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 2.5, 'S' => 0.5},
			'TTC' => {'N' => 2, 'S' => 1},
			'TTG' => {'N' => 2.5, 'S' => 0.5},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'GTA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.5, 'S' => 0.5},
			'AAG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AAT' => {'N' => 2.5, 'S' => 0.5},
			'ACA' => {'N' => 2, 'S' => 0},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 1},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 2, 'S' => 0},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 0},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 0},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 0},
			'GAC' => {'N' => 1.5, 'S' => 0.5},
			'GAG' => {'N' => 1, 'S' => 1},
			'GAT' => {'N' => 1.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 0},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 0},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 0, 'S' => 1},
			'GTG' => {'N' => 0, 'S' => 1},
			'GTT' => {'N' => 0, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.5, 'S' => 0.5},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.5, 'S' => 0.5},
			'TCA' => {'N' => 2, 'S' => 0},
			'TCC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCG' => {'N' => 2, 'S' => 1},
			'TCT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.25, 'S' => 0.75},
			'TGG' => {'N' => 2, 'S' => 1},
			'TGT' => {'N' => 2.25, 'S' => 0.75},
			'TTA' => {'N' => 1, 'S' => 0},
			'TTC' => {'N' => 1.5, 'S' => 0.5},
			'TTG' => {'N' => 1, 'S' => 1},
			'TTT' => {'N' => 1.5, 'S' => 0.5},
		},
		'GTC' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2.5, 'S' => 0.5},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 1, 'S' => 1},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 0},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 1.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 1.5, 'S' => 0.5},
			'GAT' => {'N' => 1, 'S' => 1},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 0},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 0},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 0, 'S' => 1},
			'GTG' => {'N' => 0, 'S' => 1},
			'GTT' => {'N' => 0, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 1},
			'TCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCC' => {'N' => 2, 'S' => 0},
			'TCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCT' => {'N' => 2, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 2.5, 'S' => 0.5},
			'TGT' => {'N' => 2, 'S' => 1},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'GTG' => {
			'AAA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AAC' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'ACA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'AGC' => {'N' => 2.5, 'S' => 0.5},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 2.5, 'S' => 0.5},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1.5, 'S' => 0.5},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 1.5, 'S' => 0.5},
			'CAA' => {'N' => 2, 'S' => 1},
			'CAC' => {'N' => 2.5, 'S' => 0.5},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.5, 'S' => 0.5},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2, 'S' => 0},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 0},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 0},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 1, 'S' => 1},
			'GAC' => {'N' => 1.5, 'S' => 0.5},
			'GAG' => {'N' => 1, 'S' => 0},
			'GAT' => {'N' => 1.5, 'S' => 0.5},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 0},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 0},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 0, 'S' => 1},
			'GTC' => {'N' => 0, 'S' => 1},
			'GTT' => {'N' => 0, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2.5, 'S' => 0.5},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2.5, 'S' => 0.5},
			'TCA' => {'N' => 2, 'S' => 1},
			'TCC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCG' => {'N' => 2, 'S' => 0},
			'TCT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2.5, 'S' => 0.5},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 2.5, 'S' => 0.5},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 1.5, 'S' => 0.5},
			'TTG' => {'N' => 1, 'S' => 0},
			'TTT' => {'N' => 1.5, 'S' => 0.5},
		},
		'GTT' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2.5, 'S' => 0.5},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 1, 'S' => 1},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2, 'S' => 1},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 0},
			'GAA' => {'N' => 1.5, 'S' => 0.5},
			'GAC' => {'N' => 1, 'S' => 1},
			'GAG' => {'N' => 1.5, 'S' => 0.5},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 0},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1, 'S' => 1},
			'GGT' => {'N' => 1, 'S' => 0},
			'GTA' => {'N' => 0, 'S' => 1},
			'GTC' => {'N' => 0, 'S' => 1},
			'GTG' => {'N' => 0, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCC' => {'N' => 2, 'S' => 1},
			'TCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TCT' => {'N' => 2, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 1},
			'TGG' => {'N' => 2.5, 'S' => 0.5},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'TAA' => { # STOP
			'AAA' => {'N' => 0, 'S' => 0}, # STOP
			'AAC' => {'N' => 0, 'S' => 0}, # STOP
			'AAG' => {'N' => 0, 'S' => 0}, # STOP
			'AAT' => {'N' => 0, 'S' => 0}, # STOP
			'ACA' => {'N' => 0, 'S' => 0}, # STOP
			'ACC' => {'N' => 0, 'S' => 0}, # STOP
			'ACG' => {'N' => 0, 'S' => 0}, # STOP
			'ACT' => {'N' => 0, 'S' => 0}, # STOP
			'AGA' => {'N' => 0, 'S' => 0}, # STOP
			'AGC' => {'N' => 0, 'S' => 0}, # STOP
			'AGG' => {'N' => 0, 'S' => 0}, # STOP
			'AGT' => {'N' => 0, 'S' => 0}, # STOP
			'ATA' => {'N' => 0, 'S' => 0}, # STOP
			'ATC' => {'N' => 0, 'S' => 0}, # STOP
			'ATG' => {'N' => 0, 'S' => 0}, # STOP
			'ATT' => {'N' => 0, 'S' => 0}, # STOP
			'CAA' => {'N' => 0, 'S' => 0}, # STOP
			'CAC' => {'N' => 0, 'S' => 0}, # STOP
			'CAG' => {'N' => 0, 'S' => 0}, # STOP
			'CAT' => {'N' => 0, 'S' => 0}, # STOP
			'CCA' => {'N' => 0, 'S' => 0}, # STOP
			'CCC' => {'N' => 0, 'S' => 0}, # STOP
			'CCG' => {'N' => 0, 'S' => 0}, # STOP
			'CCT' => {'N' => 0, 'S' => 0}, # STOP
			'CGA' => {'N' => 0, 'S' => 0}, # STOP
			'CGC' => {'N' => 0, 'S' => 0}, # STOP
			'CGG' => {'N' => 0, 'S' => 0}, # STOP
			'CGT' => {'N' => 0, 'S' => 0}, # STOP
			'CTA' => {'N' => 0, 'S' => 0}, # STOP
			'CTC' => {'N' => 0, 'S' => 0}, # STOP
			'CTG' => {'N' => 0, 'S' => 0}, # STOP
			'CTT' => {'N' => 0, 'S' => 0}, # STOP
			'GAA' => {'N' => 0, 'S' => 0}, # STOP
			'GAC' => {'N' => 0, 'S' => 0}, # STOP
			'GAG' => {'N' => 0, 'S' => 0}, # STOP
			'GAT' => {'N' => 0, 'S' => 0}, # STOP
			'GCA' => {'N' => 0, 'S' => 0}, # STOP
			'GCC' => {'N' => 0, 'S' => 0}, # STOP
			'GCG' => {'N' => 0, 'S' => 0}, # STOP
			'GCT' => {'N' => 0, 'S' => 0}, # STOP
			'GGA' => {'N' => 0, 'S' => 0}, # STOP
			'GGC' => {'N' => 0, 'S' => 0}, # STOP
			'GGG' => {'N' => 0, 'S' => 0}, # STOP
			'GGT' => {'N' => 0, 'S' => 0}, # STOP
			'GTA' => {'N' => 0, 'S' => 0}, # STOP
			'GTC' => {'N' => 0, 'S' => 0}, # STOP
			'GTG' => {'N' => 0, 'S' => 0}, # STOP
			'GTT' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 0, 'S' => 0}, # STOP
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 0, 'S' => 0}, # STOP
			'TCA' => {'N' => 0, 'S' => 0}, # STOP
			'TCC' => {'N' => 0, 'S' => 0}, # STOP
			'TCG' => {'N' => 0, 'S' => 0}, # STOP
			'TCT' => {'N' => 0, 'S' => 0}, # STOP
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 0, 'S' => 0}, # STOP
			'TGG' => {'N' => '*', 'S' => '*'}, # STOP
			'TGT' => {'N' => 0, 'S' => 0}, # STOP
			'TTA' => {'N' => 0, 'S' => 0}, # STOP
			'TTC' => {'N' => 0, 'S' => 0}, # STOP
			'TTG' => {'N' => 0, 'S' => 0}, # STOP
			'TTT' => {'N' => 0, 'S' => 0}, # STOP
		},
		'TAC' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 1},
			'ACA' => {'N' => 2.25, 'S' => 0.75},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.25, 'S' => 0.75},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 3, 'S' => 0},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 3, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 3, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 1, 'S' => 0},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 1},
			'CCA' => {'N' => 2.25, 'S' => 0.75},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.25, 'S' => 0.75},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2.5, 'S' => 0.5},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 2.25, 'S' => 0.75},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2.25, 'S' => 0.75},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 1, 'S' => 0},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 1},
			'GCA' => {'N' => 2.25, 'S' => 0.75},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.25, 'S' => 0.75},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2.5, 'S' => 0.5},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 0, 'S' => 1},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 0},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'TAG' => { # STOP
			'AAA' => {'N' => 0, 'S' => 0}, # STOP
			'AAC' => {'N' => 0, 'S' => 0}, # STOP
			'AAG' => {'N' => 0, 'S' => 0}, # STOP
			'AAT' => {'N' => 0, 'S' => 0}, # STOP
			'ACA' => {'N' => 0, 'S' => 0}, # STOP
			'ACC' => {'N' => 0, 'S' => 0}, # STOP
			'ACG' => {'N' => 0, 'S' => 0}, # STOP
			'ACT' => {'N' => 0, 'S' => 0}, # STOP
			'AGA' => {'N' => 0, 'S' => 0}, # STOP
			'AGC' => {'N' => 0, 'S' => 0}, # STOP
			'AGG' => {'N' => 0, 'S' => 0}, # STOP
			'AGT' => {'N' => 0, 'S' => 0}, # STOP
			'ATA' => {'N' => 0, 'S' => 0}, # STOP
			'ATC' => {'N' => 0, 'S' => 0}, # STOP
			'ATG' => {'N' => 0, 'S' => 0}, # STOP
			'ATT' => {'N' => 0, 'S' => 0}, # STOP
			'CAA' => {'N' => 0, 'S' => 0}, # STOP
			'CAC' => {'N' => 0, 'S' => 0}, # STOP
			'CAG' => {'N' => 0, 'S' => 0}, # STOP
			'CAT' => {'N' => 0, 'S' => 0}, # STOP
			'CCA' => {'N' => 0, 'S' => 0}, # STOP
			'CCC' => {'N' => 0, 'S' => 0}, # STOP
			'CCG' => {'N' => 0, 'S' => 0}, # STOP
			'CCT' => {'N' => 0, 'S' => 0}, # STOP
			'CGA' => {'N' => 0, 'S' => 0}, # STOP
			'CGC' => {'N' => 0, 'S' => 0}, # STOP
			'CGG' => {'N' => 0, 'S' => 0}, # STOP
			'CGT' => {'N' => 0, 'S' => 0}, # STOP
			'CTA' => {'N' => 0, 'S' => 0}, # STOP
			'CTC' => {'N' => 0, 'S' => 0}, # STOP
			'CTG' => {'N' => 0, 'S' => 0}, # STOP
			'CTT' => {'N' => 0, 'S' => 0}, # STOP
			'GAA' => {'N' => 0, 'S' => 0}, # STOP
			'GAC' => {'N' => 0, 'S' => 0}, # STOP
			'GAG' => {'N' => 0, 'S' => 0}, # STOP
			'GAT' => {'N' => 0, 'S' => 0}, # STOP
			'GCA' => {'N' => 0, 'S' => 0}, # STOP
			'GCC' => {'N' => 0, 'S' => 0}, # STOP
			'GCG' => {'N' => 0, 'S' => 0}, # STOP
			'GCT' => {'N' => 0, 'S' => 0}, # STOP
			'GGA' => {'N' => 0, 'S' => 0}, # STOP
			'GGC' => {'N' => 0, 'S' => 0}, # STOP
			'GGG' => {'N' => 0, 'S' => 0}, # STOP
			'GGT' => {'N' => 0, 'S' => 0}, # STOP
			'GTA' => {'N' => 0, 'S' => 0}, # STOP
			'GTC' => {'N' => 0, 'S' => 0}, # STOP
			'GTG' => {'N' => 0, 'S' => 0}, # STOP
			'GTT' => {'N' => 0, 'S' => 0}, # STOP
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 0, 'S' => 0}, # STOP
			'TCA' => {'N' => 0, 'S' => 0}, # STOP
			'TCC' => {'N' => 0, 'S' => 0}, # STOP
			'TCG' => {'N' => 0, 'S' => 0}, # STOP
			'TCT' => {'N' => 0, 'S' => 0}, # STOP
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 0, 'S' => 0}, # STOP
			'TGG' => {'N' => 0, 'S' => 0}, # STOP
			'TGT' => {'N' => 0, 'S' => 0}, # STOP
			'TTA' => {'N' => 0, 'S' => 0}, # STOP
			'TTC' => {'N' => 0, 'S' => 0}, # STOP
			'TTG' => {'N' => 0, 'S' => 0}, # STOP
			'TTT' => {'N' => 0, 'S' => 0}, # STOP
		},
		'TAT' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 1, 'S' => 1},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 1, 'S' => 0},
			'ACA' => {'N' => 2.25, 'S' => 0.75},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.25, 'S' => 0.75},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 3, 'S' => 0},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 3, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 3, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 1, 'S' => 1},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 1, 'S' => 0},
			'CCA' => {'N' => 2.25, 'S' => 0.75},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.25, 'S' => 0.75},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2.5, 'S' => 0.5},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 2.25, 'S' => 0.75},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.25, 'S' => 0.75},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 1, 'S' => 1},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 1, 'S' => 0},
			'GCA' => {'N' => 2.25, 'S' => 0.75},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.25, 'S' => 0.75},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 2.5, 'S' => 0.5},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 0, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'TCA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.25, 'S' => 0.75},
			'AAG' => {'N' => 2, 'S' => 1},
			'AAT' => {'N' => 2.25, 'S' => 0.75},
			'ACA' => {'N' => 1, 'S' => 0},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 2.25, 'S' => 0.75},
			'AGG' => {'N' => 2, 'S' => 1},
			'AGT' => {'N' => 2.25, 'S' => 0.75},
			'ATA' => {'N' => 2, 'S' => 0},
			'ATC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CAA' => {'N' => 2, 'S' => 0},
			'CAC' => {'N' => 2.25, 'S' => 0.75},
			'CAG' => {'N' => 2, 'S' => 1},
			'CAT' => {'N' => 2.25, 'S' => 0.75},
			'CCA' => {'N' => 1, 'S' => 0},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 0},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2, 'S' => 1},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1.5, 'S' => 0.5},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 1.5, 'S' => 1.5},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.25, 'S' => 0.75},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.25, 'S' => 0.75},
			'GCA' => {'N' => 1, 'S' => 0},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2, 'S' => 0},
			'GTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTG' => {'N' => 2, 'S' => 1},
			'GTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 0, 'S' => 1},
			'TCG' => {'N' => 0, 'S' => 1},
			'TCT' => {'N' => 0, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 1, 'S' => 1},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 1, 'S' => 0},
			'TTC' => {'N' => 1.5, 'S' => 0.5},
			'TTG' => {'N' => 1, 'S' => 1},
			'TTT' => {'N' => 1.5, 'S' => 0.5},
		},
		'TCC' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 0},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2.5, 'S' => 0.5},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 2.5, 'S' => 0.5},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 2.5, 'S' => 0.5},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 0},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 0},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 0, 'S' => 1},
			'TCG' => {'N' => 0, 'S' => 1},
			'TCT' => {'N' => 0, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'TCG' => {
			'AAA' => {'N' => 2, 'S' => 1},
			'AAC' => {'N' => 2.25, 'S' => 0.75},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 2.25, 'S' => 0.75},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 0},
			'ACT' => {'N' => 1, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 1},
			'AGC' => {'N' => 2.5, 'S' => 0.5},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 2.5, 'S' => 0.5},
			'ATA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ATC' => {'N' => 2.5, 'S' => 0.5},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 2.5, 'S' => 0.5},
			'CAA' => {'N' => 2, 'S' => 1},
			'CAC' => {'N' => 2.25, 'S' => 0.75},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.25, 'S' => 0.75},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 0},
			'CCT' => {'N' => 1, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CGG' => {'N' => 2, 'S' => 0},
			'CGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTA' => {'N' => 1.5, 'S' => 1.5},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.25, 'S' => 0.75},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.25, 'S' => 0.75},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 0},
			'GCT' => {'N' => 1, 'S' => 1},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 0, 'S' => 1},
			'TCC' => {'N' => 0, 'S' => 1},
			'TCT' => {'N' => 0, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1.5, 'S' => 0.5},
			'TGG' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 1.5, 'S' => 0.5},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 1.5, 'S' => 0.5},
			'TTG' => {'N' => 1, 'S' => 0},
			'TTT' => {'N' => 1.5, 'S' => 0.5},
		},
		'TCT' => {
			'AAA' => {'N' => 2.5, 'S' => 0.5},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 2.5, 'S' => 0.5},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 1, 'S' => 1},
			'ACC' => {'N' => 1, 'S' => 1},
			'ACG' => {'N' => 1, 'S' => 1},
			'ACT' => {'N' => 1, 'S' => 0},
			'AGA' => {'N' => 2.5, 'S' => 0.5},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 2.5, 'S' => 0.5},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 2.5, 'S' => 0.5},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 1, 'S' => 1},
			'CCC' => {'N' => 1, 'S' => 1},
			'CCG' => {'N' => 1, 'S' => 1},
			'CCT' => {'N' => 1, 'S' => 0},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 1.83333333333333, 'S' => 1.16666666666667},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 2.5, 'S' => 0.5},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.5, 'S' => 0.5},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 1, 'S' => 1},
			'GCC' => {'N' => 1, 'S' => 1},
			'GCG' => {'N' => 1, 'S' => 1},
			'GCT' => {'N' => 1, 'S' => 0},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 0},
			'TCA' => {'N' => 0, 'S' => 1},
			'TCC' => {'N' => 0, 'S' => 1},
			'TCG' => {'N' => 0, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 1.5, 'S' => 0.5},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 1.5, 'S' => 0.5},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 1.5, 'S' => 0.5},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'TGA' => { # STOP
			'AAA' => {'N' => 0, 'S' => 0}, # STOP
			'AAC' => {'N' => 0, 'S' => 0}, # STOP
			'AAG' => {'N' => 0, 'S' => 0}, # STOP
			'AAT' => {'N' => 0, 'S' => 0}, # STOP
			'ACA' => {'N' => 0, 'S' => 0}, # STOP
			'ACC' => {'N' => 0, 'S' => 0}, # STOP
			'ACG' => {'N' => 0, 'S' => 0}, # STOP
			'ACT' => {'N' => 0, 'S' => 0}, # STOP
			'AGA' => {'N' => 0, 'S' => 0}, # STOP
			'AGC' => {'N' => 0, 'S' => 0}, # STOP
			'AGG' => {'N' => 0, 'S' => 0}, # STOP
			'AGT' => {'N' => 0, 'S' => 0}, # STOP
			'ATA' => {'N' => 0, 'S' => 0}, # STOP
			'ATC' => {'N' => 0, 'S' => 0}, # STOP
			'ATG' => {'N' => 0, 'S' => 0}, # STOP
			'ATT' => {'N' => 0, 'S' => 0}, # STOP
			'CAA' => {'N' => 0, 'S' => 0}, # STOP
			'CAC' => {'N' => 0, 'S' => 0}, # STOP
			'CAG' => {'N' => 0, 'S' => 0}, # STOP
			'CAT' => {'N' => 0, 'S' => 0}, # STOP
			'CCA' => {'N' => 0, 'S' => 0}, # STOP
			'CCC' => {'N' => 0, 'S' => 0}, # STOP
			'CCG' => {'N' => 0, 'S' => 0}, # STOP
			'CCT' => {'N' => 0, 'S' => 0}, # STOP
			'CGA' => {'N' => 0, 'S' => 0}, # STOP
			'CGC' => {'N' => 0, 'S' => 0}, # STOP
			'CGG' => {'N' => 0, 'S' => 0}, # STOP
			'CGT' => {'N' => 0, 'S' => 0}, # STOP
			'CTA' => {'N' => 0, 'S' => 0}, # STOP
			'CTC' => {'N' => 0, 'S' => 0}, # STOP
			'CTG' => {'N' => 0, 'S' => 0}, # STOP
			'CTT' => {'N' => 0, 'S' => 0}, # STOP
			'GAA' => {'N' => 0, 'S' => 0}, # STOP
			'GAC' => {'N' => 0, 'S' => 0}, # STOP
			'GAG' => {'N' => 0, 'S' => 0}, # STOP
			'GAT' => {'N' => 0, 'S' => 0}, # STOP
			'GCA' => {'N' => 0, 'S' => 0}, # STOP
			'GCC' => {'N' => 0, 'S' => 0}, # STOP
			'GCG' => {'N' => 0, 'S' => 0}, # STOP
			'GCT' => {'N' => 0, 'S' => 0}, # STOP
			'GGA' => {'N' => 0, 'S' => 0}, # STOP
			'GGC' => {'N' => 0, 'S' => 0}, # STOP
			'GGG' => {'N' => 0, 'S' => 0}, # STOP
			'GGT' => {'N' => 0, 'S' => 0}, # STOP
			'GTA' => {'N' => 0, 'S' => 0}, # STOP
			'GTC' => {'N' => 0, 'S' => 0}, # STOP
			'GTG' => {'N' => 0, 'S' => 0}, # STOP
			'GTT' => {'N' => 0, 'S' => 0}, # STOP
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 0, 'S' => 0}, # STOP
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 0, 'S' => 0}, # STOP
			'TCA' => {'N' => 0, 'S' => 0}, # STOP
			'TCC' => {'N' => 0, 'S' => 0}, # STOP
			'TCG' => {'N' => 0, 'S' => 0}, # STOP
			'TCT' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 0, 'S' => 0}, # STOP
			'TGG' => {'N' => 0, 'S' => 0}, # STOP
			'TGT' => {'N' => 0, 'S' => 0}, # STOP
			'TTA' => {'N' => 0, 'S' => 0}, # STOP
			'TTC' => {'N' => 0, 'S' => 0}, # STOP
			'TTG' => {'N' => 0, 'S' => 0}, # STOP
			'TTT' => {'N' => 0, 'S' => 0}, # STOP
		},
		'TGC' => {
			'AAA' => {'N' => 3, 'S' => 0},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 3, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 2.25, 'S' => 0.75},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 1, 'S' => 0},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 1},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 3, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 1},
			'CAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.75, 'S' => 0.25},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 0},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 1, 'S' => 1},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 0},
			'CTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTT' => {'N' => 2, 'S' => 1},
			'GAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.75, 'S' => 0.25},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 0},
			'GGG' => {'N' => 1.5, 'S' => 0.5},
			'GGT' => {'N' => 1, 'S' => 1},
			'GTA' => {'N' => 2.25, 'S' => 0.75},
			'GTC' => {'N' => 2, 'S' => 0},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 0},
			'TCG' => {'N' => 1.5, 'S' => 0.5},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGG' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 0, 'S' => 1},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 1},
		},
		'TGG' => {
			'AAA' => {'N' => 2, 'S' => 1},
			'AAC' => {'N' => 3, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 3, 'S' => 0},
			'ACA' => {'N' => 2, 'S' => 1},
			'ACC' => {'N' => 2.5, 'S' => 0.5},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.5, 'S' => 0.5},
			'AGA' => {'N' => 1, 'S' => 1},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 1, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 3, 'S' => 0},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 3, 'S' => 0},
			'CAA' => {'N' => 2, 'S' => 1},
			'CAC' => {'N' => 2.75, 'S' => 0.25},
			'CAG' => {'N' => 2, 'S' => 0},
			'CAT' => {'N' => 2.75, 'S' => 0.25},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCG' => {'N' => 2, 'S' => 0},
			'CCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1.5, 'S' => 0.5},
			'CGG' => {'N' => 1, 'S' => 0},
			'CGT' => {'N' => 1.5, 'S' => 0.5},
			'CTA' => {'N' => 1.5, 'S' => 1.5},
			'CTC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CTG' => {'N' => 1.5, 'S' => 0.5},
			'CTT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.75, 'S' => 0.25},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.75, 'S' => 0.25},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1.5, 'S' => 0.5},
			'GGG' => {'N' => 1, 'S' => 0},
			'GGT' => {'N' => 1.5, 'S' => 0.5},
			'GTA' => {'N' => 2, 'S' => 1},
			'GTC' => {'N' => 2.5, 'S' => 0.5},
			'GTG' => {'N' => 2, 'S' => 0},
			'GTT' => {'N' => 2.5, 'S' => 0.5},
			'TAA' => {'N' => '*', 'S' => '*'}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1.5, 'S' => 0.5},
			'TCG' => {'N' => 1, 'S' => 0},
			'TCT' => {'N' => 1.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 1, 'S' => 1},
			'TTC' => {'N' => 2, 'S' => 0},
			'TTG' => {'N' => 1, 'S' => 0},
			'TTT' => {'N' => 2, 'S' => 0},
		},
		'TGT' => {
			'AAA' => {'N' => 3, 'S' => 0},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 3, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2.25, 'S' => 0.75},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 1, 'S' => 1},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 1, 'S' => 0},
			'ATA' => {'N' => 2.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 1},
			'ATG' => {'N' => 3, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.75, 'S' => 0.25},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2, 'S' => 1},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 1, 'S' => 1},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 1, 'S' => 0},
			'CTA' => {'N' => 2, 'S' => 1},
			'CTC' => {'N' => 2, 'S' => 1},
			'CTG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTT' => {'N' => 2, 'S' => 0},
			'GAA' => {'N' => 2.66666666666667, 'S' => 0.333333333333333},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.75, 'S' => 0.25},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 1, 'S' => 1},
			'GGC' => {'N' => 1, 'S' => 1},
			'GGG' => {'N' => 1.5, 'S' => 0.5},
			'GGT' => {'N' => 1, 'S' => 0},
			'GTA' => {'N' => 2.25, 'S' => 0.75},
			'GTC' => {'N' => 2, 'S' => 1},
			'GTG' => {'N' => 2.5, 'S' => 0.5},
			'GTT' => {'N' => 2, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1.5, 'S' => 0.5},
			'TCT' => {'N' => 1, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 0, 'S' => 1},
			'TGG' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 1, 'S' => 1},
			'TTG' => {'N' => 2, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'TTA' => {
			'AAA' => {'N' => 2, 'S' => 0},
			'AAC' => {'N' => 2.75, 'S' => 0.25},
			'AAG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AAT' => {'N' => 2.75, 'S' => 0.25},
			'ACA' => {'N' => 2, 'S' => 0},
			'ACC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ACG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AGA' => {'N' => 2, 'S' => 0},
			'AGC' => {'N' => 2.75, 'S' => 0.25},
			'AGG' => {'N' => 2.25, 'S' => 0.75},
			'AGT' => {'N' => 2.75, 'S' => 0.25},
			'ATA' => {'N' => 1, 'S' => 0},
			'ATC' => {'N' => 1.5, 'S' => 0.5},
			'ATG' => {'N' => 1.5, 'S' => 0.5},
			'ATT' => {'N' => 1.5, 'S' => 0.5},
			'CAA' => {'N' => 1, 'S' => 1},
			'CAC' => {'N' => 2.25, 'S' => 0.75},
			'CAG' => {'N' => 1, 'S' => 2},
			'CAT' => {'N' => 2.25, 'S' => 0.75},
			'CCA' => {'N' => 1.5, 'S' => 0.5},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 1.5, 'S' => 1.5},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 1.25, 'S' => 1.75},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 0, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 0, 'S' => 2},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 0},
			'GAC' => {'N' => 2.75, 'S' => 0.25},
			'GAG' => {'N' => 2, 'S' => 1},
			'GAT' => {'N' => 2.75, 'S' => 0.25},
			'GCA' => {'N' => 2, 'S' => 0},
			'GCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCG' => {'N' => 2, 'S' => 1},
			'GCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGA' => {'N' => 2, 'S' => 0},
			'GGC' => {'N' => 2.5, 'S' => 0.5},
			'GGG' => {'N' => 2, 'S' => 1},
			'GGT' => {'N' => 2.5, 'S' => 0.5},
			'GTA' => {'N' => 1, 'S' => 0},
			'GTC' => {'N' => 1.5, 'S' => 0.5},
			'GTG' => {'N' => 1, 'S' => 1},
			'GTT' => {'N' => 1.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 0},
			'TCC' => {'N' => 1.5, 'S' => 0.5},
			'TCG' => {'N' => 1, 'S' => 1},
			'TCT' => {'N' => 1.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 1, 'S' => 1},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 0, 'S' => 1},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'TTC' => {
			'AAA' => {'N' => 2.75, 'S' => 0.25},
			'AAC' => {'N' => 2, 'S' => 0},
			'AAG' => {'N' => 3, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 1},
			'ACA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ACC' => {'N' => 2, 'S' => 0},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 1},
			'AGA' => {'N' => 2.75, 'S' => 0.25},
			'AGC' => {'N' => 2, 'S' => 0},
			'AGG' => {'N' => 3, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 1},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1, 'S' => 0},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 1},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 0},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 1},
			'CCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCC' => {'N' => 2, 'S' => 0},
			'CCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 0},
			'CGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGT' => {'N' => 2, 'S' => 1},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 0},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2.75, 'S' => 0.25},
			'GAC' => {'N' => 2, 'S' => 0},
			'GAG' => {'N' => 2.75, 'S' => 0.25},
			'GAT' => {'N' => 2, 'S' => 1},
			'GCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCC' => {'N' => 2, 'S' => 0},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 1},
			'GGA' => {'N' => 2.25, 'S' => 0.75},
			'GGC' => {'N' => 2, 'S' => 0},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 1},
			'GTA' => {'N' => 1.5, 'S' => 0.5},
			'GTC' => {'N' => 1, 'S' => 0},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 1},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 1},
			'TCA' => {'N' => 1.5, 'S' => 0.5},
			'TCC' => {'N' => 1, 'S' => 0},
			'TCG' => {'N' => 1.5, 'S' => 0.5},
			'TCT' => {'N' => 1, 'S' => 1},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 0},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 1},
			'TTA' => {'N' => 1, 'S' => 0},
			'TTG' => {'N' => 1, 'S' => 0},
			'TTT' => {'N' => 0, 'S' => 1},
		},
		'TTG' => {
			'AAA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'AAC' => {'N' => 3, 'S' => 0},
			'AAG' => {'N' => 2, 'S' => 0},
			'AAT' => {'N' => 3, 'S' => 0},
			'ACA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'ACC' => {'N' => 2.5, 'S' => 0.5},
			'ACG' => {'N' => 2, 'S' => 0},
			'ACT' => {'N' => 2.5, 'S' => 0.5},
			'AGA' => {'N' => 2.25, 'S' => 0.75},
			'AGC' => {'N' => 3, 'S' => 0},
			'AGG' => {'N' => 2, 'S' => 0},
			'AGT' => {'N' => 3, 'S' => 0},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 2, 'S' => 0},
			'ATG' => {'N' => 1, 'S' => 0},
			'ATT' => {'N' => 2, 'S' => 0},
			'CAA' => {'N' => 1, 'S' => 2},
			'CAC' => {'N' => 2.25, 'S' => 0.75},
			'CAG' => {'N' => 1, 'S' => 1},
			'CAT' => {'N' => 2.25, 'S' => 0.75},
			'CCA' => {'N' => 1.5, 'S' => 1.5},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 1.5, 'S' => 0.5},
			'CCT' => {'N' => 2, 'S' => 1},
			'CGA' => {'N' => 1.25, 'S' => 1.75},
			'CGC' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CGG' => {'N' => 1.5, 'S' => 0.5},
			'CGT' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CTA' => {'N' => 0, 'S' => 2},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 0, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 1},
			'GAA' => {'N' => 2, 'S' => 1},
			'GAC' => {'N' => 2.75, 'S' => 0.25},
			'GAG' => {'N' => 2, 'S' => 0},
			'GAT' => {'N' => 2.75, 'S' => 0.25},
			'GCA' => {'N' => 2, 'S' => 1},
			'GCC' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCG' => {'N' => 2, 'S' => 0},
			'GCT' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GGA' => {'N' => 2, 'S' => 1},
			'GGC' => {'N' => 2.5, 'S' => 0.5},
			'GGG' => {'N' => 2, 'S' => 0},
			'GGT' => {'N' => 2.5, 'S' => 0.5},
			'GTA' => {'N' => 1, 'S' => 1},
			'GTC' => {'N' => 1.5, 'S' => 0.5},
			'GTG' => {'N' => 1, 'S' => 0},
			'GTT' => {'N' => 1.5, 'S' => 0.5},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 2, 'S' => 0},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 2, 'S' => 0},
			'TCA' => {'N' => 1, 'S' => 1},
			'TCC' => {'N' => 1.5, 'S' => 0.5},
			'TCG' => {'N' => 1, 'S' => 0},
			'TCT' => {'N' => 1.5, 'S' => 0.5},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 2, 'S' => 0},
			'TGG' => {'N' => 1, 'S' => 0},
			'TGT' => {'N' => 2, 'S' => 0},
			'TTA' => {'N' => 0, 'S' => 1},
			'TTC' => {'N' => 1, 'S' => 0},
			'TTT' => {'N' => 1, 'S' => 0},
		},
		'TTT' => {
			'AAA' => {'N' => 2.75, 'S' => 0.25},
			'AAC' => {'N' => 2, 'S' => 1},
			'AAG' => {'N' => 3, 'S' => 0},
			'AAT' => {'N' => 2, 'S' => 0},
			'ACA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'ACC' => {'N' => 2, 'S' => 1},
			'ACG' => {'N' => 2.5, 'S' => 0.5},
			'ACT' => {'N' => 2, 'S' => 0},
			'AGA' => {'N' => 2.75, 'S' => 0.25},
			'AGC' => {'N' => 2, 'S' => 1},
			'AGG' => {'N' => 3, 'S' => 0},
			'AGT' => {'N' => 2, 'S' => 0},
			'ATA' => {'N' => 1.5, 'S' => 0.5},
			'ATC' => {'N' => 1, 'S' => 1},
			'ATG' => {'N' => 2, 'S' => 0},
			'ATT' => {'N' => 1, 'S' => 0},
			'CAA' => {'N' => 2.5, 'S' => 0.5},
			'CAC' => {'N' => 2, 'S' => 1},
			'CAG' => {'N' => 2.5, 'S' => 0.5},
			'CAT' => {'N' => 2, 'S' => 0},
			'CCA' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCC' => {'N' => 2, 'S' => 1},
			'CCG' => {'N' => 2.16666666666667, 'S' => 0.833333333333333},
			'CCT' => {'N' => 2, 'S' => 0},
			'CGA' => {'N' => 2, 'S' => 1},
			'CGC' => {'N' => 2, 'S' => 1},
			'CGG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'CGT' => {'N' => 2, 'S' => 0},
			'CTA' => {'N' => 1, 'S' => 1},
			'CTC' => {'N' => 1, 'S' => 1},
			'CTG' => {'N' => 1, 'S' => 1},
			'CTT' => {'N' => 1, 'S' => 0},
			'GAA' => {'N' => 2.75, 'S' => 0.25},
			'GAC' => {'N' => 2, 'S' => 1},
			'GAG' => {'N' => 2.75, 'S' => 0.25},
			'GAT' => {'N' => 2, 'S' => 0},
			'GCA' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCC' => {'N' => 2, 'S' => 1},
			'GCG' => {'N' => 2.33333333333333, 'S' => 0.666666666666667},
			'GCT' => {'N' => 2, 'S' => 0},
			'GGA' => {'N' => 2.25, 'S' => 0.75},
			'GGC' => {'N' => 2, 'S' => 1},
			'GGG' => {'N' => 2.5, 'S' => 0.5},
			'GGT' => {'N' => 2, 'S' => 0},
			'GTA' => {'N' => 1.5, 'S' => 0.5},
			'GTC' => {'N' => 1, 'S' => 1},
			'GTG' => {'N' => 1.5, 'S' => 0.5},
			'GTT' => {'N' => 1, 'S' => 0},
			'TAA' => {'N' => 0, 'S' => 0}, # STOP
			'TAC' => {'N' => 1, 'S' => 1},
			'TAG' => {'N' => 0, 'S' => 0}, # STOP
			'TAT' => {'N' => 1, 'S' => 0},
			'TCA' => {'N' => 1.5, 'S' => 0.5},
			'TCC' => {'N' => 1, 'S' => 1},
			'TCG' => {'N' => 1.5, 'S' => 0.5},
			'TCT' => {'N' => 1, 'S' => 0},
			'TGA' => {'N' => 0, 'S' => 0}, # STOP
			'TGC' => {'N' => 1, 'S' => 1},
			'TGG' => {'N' => 2, 'S' => 0},
			'TGT' => {'N' => 1, 'S' => 0},
			'TTA' => {'N' => 1, 'S' => 0},
			'TTC' => {'N' => 0, 'S' => 1},
			'TTG' => {'N' => 1, 'S' => 0}
		}
	);
	
	my $num_N_diffs = $all_diffs_hh{$codon1}->{$codon2}->{N};
	my $num_S_diffs = $all_diffs_hh{$codon1}->{$codon2}->{S};
	
	my @diffs_arr = ($num_N_diffs,$num_S_diffs);
	
	return @diffs_arr;
}


#########################################################################################
# Perform the sliding window codon analysis using the finished codon_results.txt file
sub sliding_window {
	chdir('SNPGenie_Results');
	my ($sliding_window_size) = @_;
	my $results_file = 'codon_results.txt';
	
	# Generate new file name prefix
	my $new_file_name = "sliding_window_length$slidingwindow\_results.txt";
	if (-e "$new_file_name") { # Can also use "./SNPGenie_Results"; use "-d" to check directory
		die "\n\n## WARNING:\n## The file $new_file_name already exists.\n## Please ".
			"rename or move this file so a new one ".
			"can be created.\n\n";
	}
	
	my @header_names_arr = &get_header_names($results_file,$results_file);
	#print "\n@header_names_arr\n\n";
	
	my $index_file;
	my $index_product;
	my $index_site;
	my $index_codon;
	my $index_N_diffs_codon;
	my $index_S_diffs_codon;
#	my $index_N_diffs_site;
#	my $index_S_diffs_site;
	my $index_N_sites_codon;
	my $index_S_sites_codon;
#	my $index_N_sites_site;
#	my $index_S_sites_site;
	my $index_N_sites_ref;
	my $index_S_sites_ref;
	my $index_N_diffs_ref;
	my $index_S_diffs_ref;
#	my $index_piN;
#	my $index_piS;
#	my $index_piN_over_piS;
	
#	my $index_num_overlap_ORFs;
#	my $index_possible_stops;
#	my $index_mean_dN_ref;
#	my $index_mean_dS_ref;
	
	# Determine the index of each column
	for (my $i=0; $i<scalar(@header_names_arr); $i++) {
		if ($header_names_arr[$i] eq 'file') {
			$index_file = $i;
		} elsif ($header_names_arr[$i] eq 'product') {
			$index_product = $i;
		} elsif ($header_names_arr[$i] eq 'site') {
			$index_site = $i;
		} elsif ($header_names_arr[$i] eq 'codon') {
			$index_codon = $i;
		} elsif ($header_names_arr[$i] eq 'mean_nonsyn_diffs') {
			$index_N_diffs_codon = $i;
		} elsif ($header_names_arr[$i] eq 'mean_syn_diffs') {
			$index_S_diffs_codon = $i;
		} elsif ($header_names_arr[$i] eq 'nonsyn_sites') {
			$index_N_sites_codon = $i;
		} elsif ($header_names_arr[$i] eq 'syn_sites') {
			$index_S_sites_codon = $i;
		} elsif ($header_names_arr[$i] eq 'nonsyn_sites_ref') {
			$index_N_sites_ref = $i;
		} elsif ($header_names_arr[$i] eq 'syn_sites_ref') {
			$index_S_sites_ref = $i;
		} elsif ($header_names_arr[$i] eq "mean_nonsyn_diffs_vs_ref") {
			$index_N_diffs_ref = $i;
		} elsif ($header_names_arr[$i] eq "mean_syn_diffs_vs_ref") {
			$index_S_diffs_ref = $i;
		}
	}
	
	my %results_file_data_hh;
	# Will store file->product->site->codon/N_diffs_codon/S_diffs_codon/...
	
	#my %codon_sites_hh;
	# Will store file->product->
	
	my $line = 0;
	open (INFILE, $results_file);
	while (<INFILE>) {
		if($line == 0) {
			$line++;
		} else {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			my @line_arr = split(/\t/,$_,-1);
			
			my $file = $line_arr[$index_file];
			my $product = $line_arr[$index_product];
			my $site = $line_arr[$index_site];
			my $codon = $line_arr[$index_codon];
			my $N_diffs_codon = $line_arr[$index_N_diffs_codon];
			my $S_diffs_codon = $line_arr[$index_S_diffs_codon];
#			my $N_diffs_site = $line_arr[$index_N_diffs_site];
#			my $S_diffs_site = $line_arr[$index_S_diffs_site];
			my $N_sites_codon = $line_arr[$index_N_sites_codon];
			my $S_sites_codon = $line_arr[$index_S_sites_codon];
#			my $N_sites_site = $line_arr[$index_N_sites_site];
#			my $S_sites_site = $line_arr[$index_S_sites_site];
			my $N_sites_ref = $line_arr[$index_N_sites_ref];
			my $S_sites_ref = $line_arr[$index_S_sites_ref];
			my $N_diffs_ref = $line_arr[$index_N_diffs_ref];
			my $S_diffs_ref = $line_arr[$index_S_diffs_ref];
#			my $piN = $line_arr[$index_piN];
#			my $piS = $line_arr[$index_piS];
#			my $piN_over_piS = $line_arr[$index_piN_over_piS];
			
#			my $num_overlap_ORFs = $line_arr[$index_num_overlap_ORFs];
#			my $possible_stops = $line_arr[$index_possible_stops];
#			my $mean_dN_ref = $line_arr[$index_mean_dN_ref];
#			my $mean_dS_ref = $line_arr[$index_mean_dS_ref];
			
			#push(@codon_site_arr,$site);
			
			$results_file_data_hh{$file}->{$product}->{$site}->{codon} = $codon;
			$results_file_data_hh{$file}->{$product}->{$site}->{N_diffs_codon} = $N_diffs_codon;
			$results_file_data_hh{$file}->{$product}->{$site}->{S_diffs_codon} = $S_diffs_codon;
#			$results_file_data_hh{$file}->{$product}->{$site}->{N_diffs_site} = $N_diffs_site;
#			$results_file_data_hh{$file}->{$product}->{$site}->{S_diffs_site} = $S_diffs_site;
			$results_file_data_hh{$file}->{$product}->{$site}->{N_sites_codon} = $N_sites_codon;
			$results_file_data_hh{$file}->{$product}->{$site}->{S_sites_codon} = $S_sites_codon;
#			$results_file_data_hh{$file}->{$product}->{$site}->{N_sites_site} = $N_sites_site;
#			$results_file_data_hh{$file}->{$product}->{$site}->{S_sites_site} = $S_sites_site;
			$results_file_data_hh{$file}->{$product}->{$site}->{N_sites_ref} = $N_sites_ref;
			$results_file_data_hh{$file}->{$product}->{$site}->{S_sites_ref} = $S_sites_ref;
			$results_file_data_hh{$file}->{$product}->{$site}->{N_diffs_ref} = $N_diffs_ref;
			$results_file_data_hh{$file}->{$product}->{$site}->{S_diffs_ref} = $S_diffs_ref;
#			$results_file_data_hh{$file}->{$product}->{$site}->{piN} = $piN;
#			$results_file_data_hh{$file}->{$product}->{$site}->{piS} = $piS;
#			$results_file_data_hh{$file}->{$product}->{$site}->{piN_over_piS} = $piN_over_piS;
			
#			$results_file_data_hh{$file}->{$product}->{$site}->{num_overlap_ORFs} = $num_overlap_ORFs;
#			$results_file_data_hh{$file}->{$product}->{$site}->{possible_stops} = $possible_stops;
#			$results_file_data_hh{$file}->{$product}->{$site}->{mean_dN_ref} = $mean_dN_ref;
#			$results_file_data_hh{$file}->{$product}->{$site}->{mean_dS_ref} = $mean_dS_ref;
			
		}
	}
	close INFILE;
	
	my %sliding_window_results_hh;
	
	my @sorted_files = sort (keys %results_file_data_hh);
	#print "\n@sorted_files\n\n";
	
	my $first_print = 0;
	
	foreach my $file (@sorted_files) {
		my @products = keys (%{$results_file_data_hh{$file}});
		my @sorted_products = sort (@products);
		#print "\n@sorted_products\n\n";
		
		foreach my $product (@sorted_products) {
			my @sites = keys (%{$results_file_data_hh{$file}->{$product}});
			my @sorted_sites = sort {$a <=> $b} (@sites);
			#print "\n@sorted_sites\n\n";
			my $num_sorted_sites = scalar @sorted_sites;
			my $contiguous_zero_sites = 0;
			
			for (my $site_index = 0; $site_index <= ($num_sorted_sites - $sliding_window_size); $site_index++) {
				if ($num_sorted_sites >= $sliding_window_size) {
						my $site = $sorted_sites[$site_index];
						#print "Site: $site\n";
						
						my @site_indices_in_window_arr;
						
						for (my $site_to_add_index = $site_index; $site_to_add_index < ($site_index + $sliding_window_size); $site_to_add_index++) {
							push(@site_indices_in_window_arr,$site_to_add_index);
						}
						
						#print "\n\n@site_indices_in_window_arr\n\n";
						
						my $last_site = ($sorted_sites[$site_indices_in_window_arr[-1]]+2);
						my $last_codon = $results_file_data_hh{$file}->{$product}->{($last_site-2)}->{codon};
						
						my $first_codon = $results_file_data_hh{$file}->{$product}->{$site}->{codon};
						my $sum_N_diffs_codon;
						my $sum_S_diffs_codon;
#						my $sum_N_diffs_site;
#						my $sum_S_diffs_site;
						my $sum_N_sites_codon;
						my $sum_S_sites_codon;
#						my $sum_N_sites_site;
#						my $sum_S_sites_site;
						my $sum_N_sites_ref;
						my $sum_S_sites_ref;
						my $sum_N_diffs_ref;
						my $sum_S_diffs_ref;
						
						foreach (@site_indices_in_window_arr) { # Indices refer to the index OF the site NUMBER in @sorted_sites
							my $this_site = $sorted_sites[$_];
							$sum_N_diffs_codon += $results_file_data_hh{$file}->{$product}->{$this_site}->{N_diffs_codon};
							$sum_S_diffs_codon += $results_file_data_hh{$file}->{$product}->{$this_site}->{S_diffs_codon};
#							$sum_N_diffs_site += $results_file_data_hh{$file}->{$product}->{$this_site}->{N_diffs_site};
#							$sum_S_diffs_site += $results_file_data_hh{$file}->{$product}->{$this_site}->{S_diffs_site};
							$sum_N_sites_codon += $results_file_data_hh{$file}->{$product}->{$this_site}->{N_sites_codon};
							$sum_S_sites_codon += $results_file_data_hh{$file}->{$product}->{$this_site}->{S_sites_codon};
#							$sum_N_sites_site += $results_file_data_hh{$file}->{$product}->{$this_site}->{N_sites_site};
#							$sum_S_sites_site += $results_file_data_hh{$file}->{$product}->{$this_site}->{S_sites_site};
							$sum_N_sites_ref += $results_file_data_hh{$file}->{$product}->{$this_site}->{N_sites_ref};
							$sum_S_sites_ref += $results_file_data_hh{$file}->{$product}->{$this_site}->{S_sites_ref};
							$sum_N_diffs_ref += $results_file_data_hh{$file}->{$product}->{$this_site}->{N_diffs_ref};
							$sum_S_diffs_ref += $results_file_data_hh{$file}->{$product}->{$this_site}->{S_diffs_ref};
							#print "Site $this_site\nsum_N_diffs_codon $sum_N_diffs_codon\nsum_N_sites_codon $sum_N_sites_codon\n\n";
						}
						
						#if ($sorted_sites[$site_index] == 156 && $product eq 'HA') {
						#	print "Sum N Diff Codon: $sum_N_diffs_codon";
						#}
						
						#if () {
						#	
						#}
						
						my $window_codon_piN;
						if ($sum_N_sites_codon > 0) {
							$window_codon_piN = ($sum_N_diffs_codon / $sum_N_sites_codon);
						} else {
							$window_codon_piN = '*';
						}
						
						my $window_codon_piS;
						if ($sum_S_sites_codon > 0) {
							$window_codon_piS = ($sum_S_diffs_codon / $sum_S_sites_codon);
						} else {
							$window_codon_piS = '*';
						}
						
						my $window_ref_piN;
						if ($sum_N_sites_codon > 0) {
							$window_ref_piN = ($sum_N_diffs_ref / $sum_N_sites_codon);
						} else {
							$window_ref_piN = '*';
						}
						
						my $window_ref_piS;
						if ($sum_S_sites_codon > 0) {
							$window_ref_piS = ($sum_S_diffs_ref / $sum_S_sites_codon);
						} else {
							$window_ref_piS = '*';
						}
						
#						my $window_site_piN;
#						if ($sum_N_sites_site > 0) {
#							$window_site_piN = ($sum_N_diffs_site / $sum_N_sites_site);
#						} else {
#							$window_site_piN = '*';
#						}
#						
#						my $window_site_piS;
#						if ($sum_S_sites_site > 0) {
#							$window_site_piS = ($sum_S_diffs_site / $sum_S_sites_site);
#						} else {
#							$window_site_piS = '*';
#						}
						
						my $window_codon_based_ratio;
						if ($window_codon_piS > 0 && $window_codon_piS ne '*') {
							$window_codon_based_ratio = ($window_codon_piN / $window_codon_piS);
						} else {
							$window_codon_based_ratio = '*';
						}
	
#						my $window_site_based_ratio;
#						if ($window_site_piS > 0 && $window_site_piS ne '*') {
#							$window_site_based_ratio = ($window_site_piN / $window_site_piS);
#						} else {
#							$window_site_based_ratio = '*';
#						} 

						my $window_ref_ratio;
						if ($window_ref_piS > 0 && $window_ref_piS ne '*') {
							$window_ref_ratio = ($window_ref_piN / $window_ref_piS);
						} else {
							$window_ref_ratio = '*';
						}
						
						open(OUTFILE,">>$new_file_name");
						if ($first_print == 0) { 
							$first_print += 1; # Increment so header is only printed once
							print OUTFILE "file\tproduct\tfirst_site\tfirst_codon\tlast_site\tlast_codon\t".
								"sum_nonsyn_diffs\tsum_syn_diffs\t".
								#"Sum Nonsyn Diffs (Site-Based)\tSum Syn Diffs (Site-Based)\t".
								"sum_nonsyn_sites\tsum_syn_sites\t".
								#"Sum Nonsyn Sites (Site-Based)\tSum Syn Sites (Site-Based)\t".
								"piN\tpiS\t".
								#"πN (Site-Based)\tπS (Site-Based)\t".
								"piN/piS\t".
								#"πN/πS (Site-Based)\n";
								"piN_vs_ref\tpiS_vs_ref\tpiN/piS_vs_ref\n";
						}
						
						print OUTFILE "$file\t$product\t$site\t$first_codon\t$last_site\t$last_codon\t".
							"$sum_N_diffs_codon\t$sum_S_diffs_codon\t".
							#"$sum_N_diffs_site\t$sum_S_diffs_site\t".
							"$sum_N_sites_codon\t$sum_S_sites_codon\t".
							#"$sum_N_sites_site\t$sum_S_sites_site\t".
							"$window_codon_piN\t$window_codon_piS\t".
							#"$window_site_piN\t$window_site_piS\t".
							"$window_codon_based_ratio\t".
							#"$window_site_based_ratio\n";
							"$window_ref_piN\t$window_ref_piS\t$window_ref_ratio\n";
						
						close OUTFILE;
						
				} else {
					print "\n\n## WARNING:\n## For file $file, product $product,\n## there are ".
					"not enough sites to perform a sliding window of $sliding_window_size.\n\n";
				}
			}
		}
	}
	chdir('..');
}


#########################################################################################
# End the program by notifying the screen at command line
sub end_the_program {
	my $time2 = time;
	my $local_time2 = localtime;
	
	my $time_diff = ($time2 - $time1);
	my $time_diff_rounded = sprintf("%.2f",$time_diff);
	my $mins_elapsed = ($time_diff / 60);
	my $whole_mins_elapsed = int($mins_elapsed);
	my $whole_mins_in_secs = ($whole_mins_elapsed * 60);
	my $secs_remaining = ($time_diff - $whole_mins_in_secs);
	my $secs_remaining_rounded = sprintf("%.2f",$secs_remaining);
	
	chdir('SNPGenie_Results');
	open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
	print ERROR_FILE "NA\tNA\tNA\t".
			"SNPGenie completed at local time $local_time2. The process took $time_diff_rounded secs, i.e., ".
			"$whole_mins_elapsed mins and $secs_remaining_rounded secs\n";
	close ERROR_FILE;
	chdir('..');
	
	if($vcfformat == 4) { # remove the temp SNP reports
		foreach(@temp_vcf4_file_names) {
			unlink $_;
		}
	}
	
	print "\n################################################################################".
		"\n##                      SNPGenie completed successfully.                      ##".
		"\n##             Please find results in the SNPGenie_Results folder.            ##\n".
		"################################################################################".
		"\n\n\n"; 
}