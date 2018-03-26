#! /usr/bin/perl

#AMNH-HUXLEY-PBS is #! /usr/local/software/PERL/perl-5.26.0/bin/perl
#NYU-PRINCE-SLURM is #! /share/apps/perl/5.24.0/intel/bin/perl 
#PC is #! /usr/bin/perl

# POST-PROCESSING for PhyloGeneious.

# PROGRAM: Recalculate whole-product and sliding window dN/dS values FOR BETWEEN-GROUP 
# 		ANALYSES using a given window and step size for all families/products, and criteria 
# 		for number of taxa present in the alignment overall (total) or in each group (group).
# REQUIRES: the existence of "between_group_codon_results.txt" files, produced from an
#		initial run of SNPGenie, in all family directories (i.e., 
#		$OID_USER_DIR/data/*/between_group_codon_results.txt).

#########################################################################################
# CALL FORMAT:
#########################################################################################
# snpgenie_between_group_processor.pl <input_file> <codon_min_taxa_total> <codon_min_taxa_group> <sw_size> <sw_step_size> <num_bootstraps> <procs_per_node>
#########################################################################################

#########################################################################################
# EXAMPLE CALLS:
#########################################################################################
# snpgenie_between_group_processor.pl --between_group_codon_file=between_group_codon_results.txt --sliding_window_size=8 --sliding_window_step=1 --num_bootstraps=1000 --procs_per_node=4
# snpgenie_between_group_processor.pl --between_group_codon_file=between_group_codon_results.txt --sliding_window_size=1 --sliding_window_step=1 --num_bootstraps=1000 --procs_per_node=1 --codon_min_taxa_total=5 --codon_min_taxa_group=2
#########################################################################################

# OUTPUT:
# (1) between_group_sw_<#>codons_results.txt
#		If the filtering criteria (codon_min_taxa_total and/or codon_min_taxa_group) are
#		NOT met, then the columns N_sites, S_sites, N_diffs, S_diffs dN>total_dS, 
#		SE(dN-dS), and Z_value will be populated with NA values. REPLACES file by the same 
#		name if it already exists.
# (2) between_group_results_filtered.txt
#		Identical to the original between_group_results.txt, but with values calculated 
#		excluding codons which do not meet filtering criteria. REPLACES file by the same
#		name if it already exists.

# Copyright (C) 2017 Chase W. Nelson

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

# DATE CREATED: November 9, 2017

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
use Data::Dumper;
use Getopt::Long;
use List::Util qw(max);
use Parallel::ForkManager;


#########################################################################################
# INITIALIZE VARIABLES (all optional and/or defaulted)
my $between_group_codon_file;
my $sliding_window_size;
my $sliding_window_step;
my $codon_min_taxa_total;
my $codon_min_taxa_group;
my $num_bootstraps;
my $procs_per_node;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "between_group_codon_file=s" => \$between_group_codon_file,
			"sliding_window_size=i" => \$sliding_window_size,
			"sliding_window_step=i" => \$sliding_window_step,
			"codon_min_taxa_total=i" => \$codon_min_taxa_total,
			"codon_min_taxa_group=i" => \$codon_min_taxa_group,
			"num_bootstraps=i" => \$num_bootstraps,
			"procs_per_node:i" => \$procs_per_node)
			
			or die "\n### WARNING: Error in command line arguments. SNPGenie terminated.\n\n";

# Get the time
my $time1 = time;
my $local_time1 = localtime;

print "\n#########################################################################################";
print "\n##### SNPGenie filtering, sliding windows, and bootstrapping initiated at local time $local_time1\n";
print "#########################################################################################";
print "\n";

# Test and/or reset OPTIONS given the user's INPUT and RECORD PARAMETERS	
if(! $between_group_codon_file) {
	die "\n### WARNING: The --between_group_codon_file option must be provided\n".
		"### SNPGenie terminated.\n\n";
}

# Test and/or reset OPTIONS given the user's INPUT and RECORD PARAMETERS	
if(! $sliding_window_size) {
	if($sliding_window_size != 0) { # Called as a flag, but given no value
		$sliding_window_size = 10; # default behavior: 10-mer peptide
	}
} elsif($sliding_window_size < 1) {
	die "\n### WARNING: The --sliding_window_size option must be an integer ≥1\n".
		"### SNPGenie terminated.\n\n";
}

if(! $sliding_window_step) {
	if($sliding_window_step != 0) { # Called as a flag, but given no value
		$sliding_window_step = 1; # default behavior: one codon
	}
} elsif($sliding_window_step < 1) {
	die "\n### WARNING: The --sliding_window_step option must be an integer ≥1\n".
		"### SNPGenie terminated.\n\n";
}

if(! $codon_min_taxa_total) { # null or 0
	$codon_min_taxa_total = 0; # default behavior
} elsif(! $codon_min_taxa_total >= 1) {
	die "\n### WARNING: The --codon_min_taxa_total option must ≥1\n".
		"### Script terminated.\n\n";
}

if(! $codon_min_taxa_group) { # null or 0
	$codon_min_taxa_group = 0; # default behavior
} elsif(! $codon_min_taxa_group >= 1) {
	die "\n### WARNING: The --codon_min_taxa_group option must ≥1\n".
		"### Script terminated.\n\n";
}

if(! $num_bootstraps) {
	$num_bootstraps = 0; # DEFAULT is 0
	print "### BOOTSTRAP VALUE NOT GIVEN. DEFAULT IS NO BOOTSTRAPPING.\n" . 
		"### (You sure you don't want our statistics?)\n\n";
} elsif($num_bootstraps < 100) {
	die "\n### WARNING: The --num_bootstraps option must be an integer ≥100\n".
		"### SNPGenie terminated.\n\n";
}

if(! $procs_per_node) {
	$procs_per_node = 1; # DEFAULT is 1
	print "### PROCS PER NODE VALUE NOT GIVEN. DEFAULT USED: JUST 1.\n" . 
		"### (You sure you don't want our parallelization?)\n\n";
} elsif($procs_per_node < 1) {
	die "\n### WARNING: The --procs_per_node option must be an integer ≥1\n".
		"### SNPGenie terminated.\n\n";
}


#STDOUT->autoflush(1);	

#########################################################################################
# Input and output file names; same for every family directory
my $between_group_codon_file = 'between_group_codon_results.txt';
my $between_group_sw_output_file = 'between_group_sw_' . $sliding_window_size . 'codons_results.txt';
my $between_group_product_results_file = 'between_group_product_results_filtered.txt'; # later, looks like between\_group\_results\_filtered\.txt


#########################################################################################
# BEGIN PROCESSING


print "\n#########################################################################################";
print "\n##### Working on file $between_group_codon_file\n";
print "#########################################################################################";
print "\n";


if(-f "$between_group_codon_file") { # IF THE INPUT FILE EXISTS
	
	##################################################################################
	# READY THE OUTPUT FILE FOR SLIDING WINDOWS
	if(-f "$between_group_sw_output_file") { # file already exists
		warn "\n### WARNING: $between_group_sw_output_file already exists in this directory.\n".
			"### REMOVING AND REPLACING using current criteria.\n";
			
		unlink $between_group_sw_output_file;
	}
	
	open(OUTFILE_BETWEEN_SW,">$between_group_sw_output_file");
#		print OUTFILE_BETWEEN_SW "analysis\tfamily\tproduct\tgroup_1\tgroup_2\twindow\tfirst_codon\tlast_codon\t". ###???
	print OUTFILE_BETWEEN_SW "analysis\tproduct\tgroup_1\tgroup_2\twindow\tfirst_codon\tlast_codon\t". ###???
		"num_defined_codons_g1\tnum_defined_codons_g2\tmean_num_defined_codons\tmin_num_defined_codons\t".
		"N_sites\tS_sites\t".
		"N_diffs\tS_diffs\t".
		"dN\tdS\t";

	if($num_bootstraps > 1) {
		print OUTFILE_BETWEEN_SW "SE_dN\tSE_dS\t";
	}
	
	print OUTFILE_BETWEEN_SW "dN_minus_dS\tdN_over_dS";
		
	if($num_bootstraps > 1) {
		print OUTFILE_BETWEEN_SW "\tSE_dN_minus_dS\tZ_value\tsignificance\n";
	} else {
		print OUTFILE_BETWEEN_SW "\n";
	}
		
	close OUTFILE_BETWEEN_SW;
	
	##################################################################################
	# READY THE OUTPUT FILE FOR PRODUCT
	if(-f "$between_group_product_results_file") { # file already exists
		warn "\n### WARNING: $between_group_product_results_file already exists in this directory.\n".
			"### REMOVING AND REPLACING using current criteria.\n";
			
		unlink "$between_group_product_results_file";
	}
	
	open(BETWEEN_OUTFILE,">>$between_group_product_results_file");
#		print BETWEEN_OUTFILE "analysis\tfamily\tproduct\tgroup_1\tgroup_2\tN_sites\tS_sites\tN_diffs\tS_diffs\t". ###???
	print BETWEEN_OUTFILE "analysis\tproduct\tgroup_1\tgroup_2\tN_sites\tS_sites\tN_diffs\tS_diffs\t". ###???
		"dN\tdS\tJC_dN\tJC_dS\t";
		
	if($num_bootstraps > 1) {
		print BETWEEN_OUTFILE "SE_dN\tSE_dS\tSE_JC_dN\tSE_JC_dS\t";
	}

	print BETWEEN_OUTFILE "dN_minus_dS\tdNdS\tJC_dN_minus_dS\tJC_dNdS";
	
	if($num_bootstraps > 1) {
		print BETWEEN_OUTFILE "\tSE_dN_minus_dS\tZ_value\tsignificance\t" . 
			"SE_JC_dN_minus_dS\tJC_Z_value\tJC_significance\n";
	} else {
		print BETWEEN_OUTFILE "\n";
	}
	
	close BETWEEN_OUTFILE;
	
	
	##################################################################################
	# READY THE VARIABLES FOR THE WHOLE PARTITION
#		my %site_data_hh;
#		# $site_data_hh{<site>}->{<N_diffs>};
#		# $site_data_hh{<site>}->{<S_diffs>};
#		# $site_data_hh{<site>}->{<N_sites>};
#		# $site_data_hh{<site>}->{<S_sites>};
#		# $site_data_hh{<site>}->{<group_1_#_taxa>};
#		# $site_data_hh{<site>}->{<group_2_#_taxa>};
	
	##################################################################################
	# STORE BETWEEN-GROUP CODON DATA
	
	# Get header names
	my @header_names_arr = &get_header_names($between_group_codon_file,$between_group_codon_file);
	#print "@header_names_arr";
	
	my %header_indices;
	
	# Determine the index of each column
	for (my $i = 0; $i < scalar(@header_names_arr); $i++) {
		my $curr_header = $header_names_arr[$i];
		$header_indices{$curr_header} = $i;
	}
	
	# Make sure we got all the expected headers
#		my @required_between_group_headers = qw/analysis family product group_1 group_2 codon variability num_defined_codons_g1 num_defined_codons_g2 comparisons N_sites S_sites N_diffs S_diffs/; ###???
	my @required_between_group_headers = qw/analysis product group_1 group_2 codon variability num_defined_codons_g1 num_defined_codons_g2 comparisons N_sites S_sites N_diffs S_diffs/; ###???
	
	foreach(@required_between_group_headers) {
		unless(exists $header_indices{$_}) {
			die "### DIE: the header name $_ is not present in the between-group codon results file.\n\n";
		}
	}
	
	##################################################################################
	# STORE ALL FAMILY BETWEEN-GROUP CODON DATA
	my %between_group_codon_data_hh;
	my %between_group_a_p_groups;
	
	##################################################################################
	# STORE overall product values
	my %between_group_a_p_sums;
	
###		my $between_product_N_sites_sum = 0;
##		my $between_product_S_sites_sum = 0;
###		my $between_product_N_diffs_sum = 0;
##		my $between_product_S_diffs_sum = 0;
	
	##################################################################################
	# LOOP CODONS TO STORE DATA IN %between_group_codon_data_hh
	open(IN_BETWEEN_GROUP_CODON_FILE, "$between_group_codon_file") or die "\nCould not open the file $between_group_codon_file - $!\n\n";
	my $line = 0;
	while(<IN_BETWEEN_GROUP_CODON_FILE>) { # record
		if($line == 0) {
			$line++;
		} else {
			chomp;
			my @line_arr = split(/\t/,$_,-1);
			
			# Store data here
			my $analysis = $line_arr[$header_indices{'analysis'}];
#			my $analysis = "only_one";
#				my $family = $line_arr[$header_indices{'family'}];
#				my $product = $line_arr[$header_indices{'product'}];
			my $product = $line_arr[$header_indices{'product'}];
			my $codon = $line_arr[$header_indices{'codon'}];
			$codon =~ s/codon_//; # in case it's the old format
			$codon =~ s/\s//g; # in case it's the old format
			
			# Store the names of the groups in this pair
			my @curr_groups = sort ($line_arr[$header_indices{'group_1'}], $line_arr[$header_indices{'group_2'}]);
			my $this_group_1 = $curr_groups[0];
			my $this_group_2 = $curr_groups[1];
			
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{variability} = $line_arr[$header_indices{'variability'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{num_defined_codons_g1} = $line_arr[$header_indices{'num_defined_codons_g1'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{num_defined_codons_g2} = $line_arr[$header_indices{'num_defined_codons_g2'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{comparisons} = $line_arr[$header_indices{'comparisons'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{N_sites} = $line_arr[$header_indices{'N_sites'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{S_sites} = $line_arr[$header_indices{'S_sites'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{N_diffs} = $line_arr[$header_indices{'N_diffs'}];
			$between_group_codon_data_hh{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{$codon}->{S_diffs} = $line_arr[$header_indices{'S_diffs'}];
			
#				print "\n\nVariability: " . $between_group_codon_data_hh{$analysis}->{$family}->{$product}->{$codon}->{variability} . "\n\n";
#				print "\n\nN_diffs: " . $between_group_codon_data_hh{$analysis}->{$family}->{$product}->{$codon}->{N_diffs} . "\n\n";
			
#				print "\nSTORING:\nanalysis: $analysis\nfamily: $family\nproduct: $product\ncodon: $codon\nNsites: ".
#					$line_arr[$header_indices{'N_sites'}] . "\n\n";
			
			
			
			# THERE ARE MORE THAN 2.....
			# Store group names here
#			$between_group_a_p_groups{$analysis}->{$product}->{group_1} = $line_arr[$header_indices{'group_1'}];
#			$between_group_a_p_groups{$analysis}->{$product}->{group_2} = $line_arr[$header_indices{'group_2'}];
			# No longer needed because used as keys in the hh
			
			# Sum product values here, NONFILTERED
			$between_group_a_p_sums{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{N_sites} += $line_arr[$header_indices{'N_sites'}];
			$between_group_a_p_sums{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{S_sites} += $line_arr[$header_indices{'S_sites'}];
			$between_group_a_p_sums{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{N_diffs} += $line_arr[$header_indices{'N_diffs'}];
			$between_group_a_p_sums{$analysis}->{$this_group_1}->{$this_group_2}->{$product}->{S_diffs} += $line_arr[$header_indices{'S_diffs'}];
		}
	}
	close IN_BETWEEN_GROUP_CODON_FILE;
	
#### THIS DIFFERS FOR EACH BETWEEN-GROUP ANALYSIS, SO MUST BE RE-CALCULATED FOR EACH

#		my $between_product_dS;
#		if($between_product_S_sites_sum_filtered > 0) {
#			$between_product_dS = $between_product_S_diffs_sum_filtered / $between_product_S_sites_sum_filtered;
#		} else {
#			$between_product_dS = '*';
#		}
	
	##################################################################################
	# LOOP THROUGH EACH ANALYSIS, FAMILY, AND PARTITION FOR SLIDING WINDOW ANALYSIS
	foreach my $this_analysis (sort keys %between_group_codon_data_hh) {
		foreach my $group_name_i (sort keys %{$between_group_codon_data_hh{$this_analysis}}) {
			foreach my $group_name_j (sort keys %{$between_group_codon_data_hh{$this_analysis}->{$group_name_i}}) {
			
				foreach my $this_product (sort {$a <=> $b} keys %{$between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}}) {
					
					# Get group names
#					my $group_name_i = $between_group_a_p_groups{$this_analysis}->{$this_product}->{group_1};
#					my $group_name_j = $between_group_a_p_groups{$this_analysis}->{$this_product}->{group_2};
					
					
		#					print "\nTEST RETRIEVAL:\nanalysis: $this_analysis\nfamily: $this_family\nproduct: $this_product\n\n";
					
					#############################################################
					# SLIDING WINDOW WITH BOOTSTRAPS 
					# Generate and store sliding window and bootstrap analyses
					my %sliding_window_results_hh; # Will store first_codon->N_diffs_window/S_diffs_window/...
					
				#	my $sw_file_name = "product$product_id\_sw_" . $sliding_window_size . "codons\_results.txt";
				#	
				#	if (-e "$sw_file_name") {
				#		die "\n\n### WARNING:\n## The file $sw_file_name already exists.\n## Please ".
				#			"rename or move this file so a new one ".
				#			"can be created.\n\n";
				#	}
					
					my @codons = keys (%{$between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}});
					my @sorted_codons = sort {$a <=> $b} (@codons);
					my $num_sorted_codons = scalar @sorted_codons;
					
		#					print "\n\ncodons: @sorted_codons\n\n";
					
					if($sliding_window_size > $num_sorted_codons) {
						die "\n\n### WARNING:\n## For product $this_product\,\n## there are ".
							"not enough codons to perform a sliding window of $sliding_window_size.\n\n";
					}
					
					######################################################################
					# STEP THROUGH CODONS AND FILTER VALUES
					my $between_product_N_sites_sum_filtered = 0;
					my $between_product_S_sites_sum_filtered = 0;
					my $between_product_N_diffs_sum_filtered = 0;
					my $between_product_S_diffs_sum_filtered = 0;
					
					for (my $codon_index = 0; $codon_index <= $num_sorted_codons; $codon_index++) {
						my $codon_num = $codon_index + 1;
						
						my $num_defined_codons_g1 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{'num_defined_codons_g1'};
						my $num_defined_codons_g2 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{'num_defined_codons_g2'};
						
		#						print "\nnum_defined_codons_g1 is $num_defined_codons_g1\n";
						
						if($codon_min_taxa_total > 0 && $codon_min_taxa_group > 0) { # BOTH criteria have been specified
							my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
							my $taxa_group_min = $num_defined_codons_g1;
							
							if($num_defined_codons_g2 < $taxa_group_min) {
								$taxa_group_min = $num_defined_codons_g2;
							}
							
							if($taxa_total >= $codon_min_taxa_total && $taxa_group_min >= $codon_min_taxa_group) { # see if criteria are met by this codon
								$between_product_N_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_sites};
								$between_product_S_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_sites};
								$between_product_N_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_diffs};
								$between_product_S_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_diffs};
							}
						} elsif($codon_min_taxa_total > 0) { # only TOTAL min has been specified
							my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
							
							if($taxa_total >= $codon_min_taxa_total) { # see if criterion is met by this codon
								$between_product_N_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_sites};
								$between_product_S_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_sites};
								$between_product_N_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_diffs};
								$between_product_S_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_diffs};
							}
						} elsif($codon_min_taxa_group > 0) { # only GROUP min has been specified
							my $taxa_group_min = $num_defined_codons_g1;
							
							if($num_defined_codons_g2 < $taxa_group_min) {
								$taxa_group_min = $num_defined_codons_g2;
							}
							
							if($taxa_group_min >= $codon_min_taxa_group) { # see if criterion is met by this codon
								$between_product_N_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_sites};
								$between_product_S_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_sites};
								$between_product_N_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_diffs};
								$between_product_S_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_diffs};
							}
						} else { # BOTH CRITERIA 0 OR UNSPECIFIED
							$between_product_N_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_sites};
							$between_product_S_sites_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_sites};
							$between_product_N_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_diffs};
							$between_product_S_diffs_sum_filtered += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{S_diffs};
						}
					}
					
					my $between_product_dS_filtered;
					if($between_product_S_sites_sum_filtered > 0) {
						$between_product_dS_filtered = $between_product_S_diffs_sum_filtered / $between_product_S_sites_sum_filtered;
					} else {
						$between_product_dS_filtered = '*';
					}
					
					my $between_product_dN_filtered;
					if($between_product_N_sites_sum_filtered > 0) {
						$between_product_dN_filtered = $between_product_N_diffs_sum_filtered / $between_product_N_sites_sum_filtered;
					} else {
						$between_product_dN_filtered = '*';
					}
					
					# Jukes-Cantor Correction
					my $between_JC_dN = -3/4 * log(1 - ((4/3) * $between_product_dN_filtered));
					my $between_JC_dS = -3/4 * log(1 - ((4/3) * $between_product_dS_filtered));
					
					
					
		#					print "\nANALYSIS $this_analysis FAMILY $this_family PARTITION $this_product\n".
		#						"filtered dN = $between_product_dN_filtered filtered dS = $between_product_dS_filtered\n";
					
					######################################################################
					# STEP THROUGH WINDOWS, IGNORE CODONS THAT DON'T FIT CRITERIA
					# FOR EACH WINDOW
					my $window_num = 0;
					my $first_print = 0;
					
		#					print "\nProcessing sliding windows...\n";
					
					# CODON_INDEX must begin one before the smallest actual codon number... don't assume it's 1
					for (my $codon_index = 0; $codon_index <= ($num_sorted_codons - $sliding_window_size); $codon_index += $sliding_window_step) {
						$window_num++;
						
						my $first_codon = $codon_index + $sorted_codons[0];
						my $last_codon = $first_codon + $sliding_window_size - 1;
						
		#						print "#### PROCESSING WINDOW $window_num\, codons $first_codon\-$last_codon\... \n";
						
						if ($last_codon <= ($num_sorted_codons + $sorted_codons[0])) { # as long as we're in range 
							
							my $window_N_diffs;
							my $window_S_diffs;
							my $window_N_sites;
							my $window_S_sites;
											
							for (my $i = $first_codon; $i <= $last_codon; $i++) {
		#								print "\n\nTo codon $i, we add: " . $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_diffs} . "\n\n";
		#								print "\nTEST RETRIEVAL:\nanalysis: $this_analysis\nfamily: $this_family\nproduct: $this_product\n\n";
								
								my $num_defined_codons_g1 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{'num_defined_codons_g1'};
								my $num_defined_codons_g2 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{'num_defined_codons_g2'};
								
								if($codon_min_taxa_total > 0 && $codon_min_taxa_group > 0) { # BOTH criteria have been specified
									my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
									my $taxa_group_min = $num_defined_codons_g1;
									
									if($num_defined_codons_g2 < $taxa_group_min) {
										$taxa_group_min = $num_defined_codons_g2;
									}
									
									if($taxa_total >= $codon_min_taxa_total && $taxa_group_min >= $codon_min_taxa_group) { # see if criteria are met by this codon
										$window_N_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_diffs};
										$window_S_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_diffs};
										$window_N_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_sites};
										$window_S_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_sites};
									}
								} elsif($codon_min_taxa_total > 0) { # only TOTAL min has been specified
									my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
									
									if($taxa_total >= $codon_min_taxa_total) { # see if criterion is met by this codon
										$window_N_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_diffs};
										$window_S_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_diffs};
										$window_N_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_sites};
										$window_S_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_sites};
									}
								} elsif($codon_min_taxa_group > 0) { # only GROUP min has been specified
									my $taxa_group_min = $num_defined_codons_g1;
									
									if($num_defined_codons_g2 < $taxa_group_min) {
										$taxa_group_min = $num_defined_codons_g2;
									}
									
									if($taxa_group_min >= $codon_min_taxa_group) { # see if criterion is met by this codon
										$window_N_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_diffs};
										$window_S_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_diffs};
										$window_N_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_sites};
										$window_S_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_sites};
									}
								} else { # BOTH CRITERIA 0 OR UNSPECIFIED
									$window_N_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_diffs};
									$window_S_diffs += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_diffs};
									$window_N_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{N_sites};
									$window_S_sites += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$i}->{S_sites};
								}
								
							}
										
							my $window_dN;
							if ($window_N_sites > 0) {
								$window_dN = ($window_N_diffs / $window_N_sites);
							} else {
								$window_dN = '*';
							}
											
							my $window_dS;
							if ($window_S_sites > 0) {
								$window_dS = ($window_S_diffs / $window_S_sites);
							} else {
								$window_dS = '*';
							}
					
							my $window_w;
							if ($window_dS > 0 && $window_dS ne '*') {
								$window_w = ($window_dN / $window_dS);
							} else {
								$window_w = '*';
							}
							
							my $window_dN_minus_dS;
							if($window_dN >= 0 && $window_dS >= 0) {
								$window_dN_minus_dS = $window_dN - $window_dS;
							}
							
					#		print "\n\n";
							
							#################################################################################
							# BOOTSTRAP the WINDOW
							my $window_SE_dN_minus_dS = 0;
							my $window_boot_Z = 'NA';
	
							my $window_SE_dN;
							my $window_SE_dS;
							
							if($num_bootstraps > 1) {
								# MAKE BOOTSTRAP FILE
								
								#print "num_codons_in_product: $num_codons_in_product\n";
								
								my @sim_N_sites_arr;
								my @sim_S_sites_arr;
								my @sim_N_diffs_arr;
								my @sim_S_diffs_arr;
								
								for(my $bootstrap_num = 1; $bootstrap_num <= $num_bootstraps; $bootstrap_num++) {
									srand();
									
									# New bootstrap run, of $num_bootstraps
									#if($this_family == 460) {
									#	print "\nPartition $this_product \| BOOTSTRAP REPLICATE $bootstrap_num\n";
									#}
									
									my $sim_N_sites_sum = 0;
									my $sim_S_sites_sum = 0;
									my $sim_N_diffs_sum = 0;
									my $sim_S_diffs_sum = 0;
									
									# METHOD
									# (1) loop through $num_codons_in_product
									# (2) for each, select a random codon position's results: Nd, Sd, Ns, Ss. 
									# where are the aforementioned stored? As before:
									# $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$codon_num}->{N_diffs}, etc.
									
									# SAMPLE codon sites, up to the actual number of codons in window
									for(my $window_codon_index = 0; $window_codon_index < $sliding_window_size; $window_codon_index++) {
										# Choose a random codon between this window $first_codon and $last_codon
										my $random_codon_num = int(rand($last_codon - $first_codon + 1)) + $first_codon;
										
										my $num_defined_codons_g1 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{'num_defined_codons_g1'};
										my $num_defined_codons_g2 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{'num_defined_codons_g2'};
										
										if($codon_min_taxa_total > 0 && $codon_min_taxa_group > 0) { # BOTH criteria have been specified
											my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
											my $taxa_group_min = $num_defined_codons_g1;
											
											if($num_defined_codons_g2 < $taxa_group_min) {
												$taxa_group_min = $num_defined_codons_g2;
											}
											
											if($taxa_total >= $codon_min_taxa_total && $taxa_group_min >= $codon_min_taxa_group) { # see if criteria are met by this codon
												$sim_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
												$sim_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
												$sim_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
												$sim_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
											}
										} elsif($codon_min_taxa_total > 0) { # only TOTAL min has been specified
											my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
											
											if($taxa_total >= $codon_min_taxa_total) { # see if criterion is met by this codon
												$sim_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
												$sim_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
												$sim_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
												$sim_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
											}
										} elsif($codon_min_taxa_group > 0) { # only GROUP min has been specified
											my $taxa_group_min = $num_defined_codons_g1;
											
											if($num_defined_codons_g2 < $taxa_group_min) {
												$taxa_group_min = $num_defined_codons_g2;
											}
											
											if($taxa_group_min >= $codon_min_taxa_group) { # see if criterion is met by this codon
												$sim_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
												$sim_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
												$sim_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
												$sim_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
											}
										} else { # BOTH CRITERIA 0 OR UNSPECIFIED
											$sim_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
											$sim_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
											$sim_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
											$sim_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
										}
										
										
									} # finished compiling all sampled codons (cols in alignment)
									
									push(@sim_N_diffs_arr,$sim_N_diffs_sum);
									push(@sim_S_diffs_arr,$sim_S_diffs_sum);
									push(@sim_N_sites_arr,$sim_N_sites_sum);
									push(@sim_S_sites_arr,$sim_S_sites_sum);
									
									#print $random_sum . "\n";
									
								} # end last bootstrap
											
								# Recall we have $window_dN_minus_dS above
					
								# CALCULATE BOOTSTRAP STANDARD ERROR HERE; NEI & KUMAR (2000) METHOD
								my @sim_dN_minus_dS;
								my @sim_dN;
								my @sim_dS;
								
								for(my $sim_num = 0; $sim_num < scalar(@sim_N_sites_arr); $sim_num++) {
									my $this_round_dN = '*';
									if($sim_N_sites_arr[$sim_num] > 0) {
										$this_round_dN = $sim_N_diffs_arr[$sim_num] / $sim_N_sites_arr[$sim_num];
									}
									
									my $this_round_dS = '*';
									if($sim_S_sites_arr[$sim_num] > 0) {
										$this_round_dS = $sim_S_diffs_arr[$sim_num] / $sim_S_sites_arr[$sim_num];
									}
									
									my $this_round_dN_minus_dS = '*';
									if($this_round_dN >= 0 && $this_round_dS >= 0) {
										$this_round_dN_minus_dS = $this_round_dN - $this_round_dS;
									}
									
									push(@sim_dN_minus_dS, $this_round_dN_minus_dS);
									push(@sim_dN, $this_round_dN);
									push(@sim_dS, $this_round_dS);
								}
								
								$window_SE_dN_minus_dS = &standard_deviation(@sim_dN_minus_dS);
								$window_SE_dN = &standard_deviation(@sim_dN);
								$window_SE_dS = &standard_deviation(@sim_dS);
		
								if($window_SE_dN_minus_dS > 0) {
									$window_boot_Z = $window_dN_minus_dS / $window_SE_dN_minus_dS;
								}
								
								if($window_SE_dN_minus_dS == 0) {
									$window_SE_dN_minus_dS = 'NA';
								}
							
							} # END BOOTSTRAPS for WINDOW
							
							my $product_between_group_dS = '*';
							if($between_group_a_p_sums{$this_analysis}->{$this_product}->{S_sites} > 0) {
								$product_between_group_dS = $between_group_a_p_sums{$this_analysis}->{$this_product}->{S_diffs} / 
									$between_group_a_p_sums{$this_analysis}->{$this_product}->{S_sites};
							}
							
							my $window_dN_exceeds_product_dS = '';
							if($window_dN > $between_product_dS_filtered) {
								$window_dN_exceeds_product_dS = '*';
							}
							
		#							print "\nbetween_product_dS_filtered is $between_product_dS_filtered\n\n";
							
							my $min_defined_codons_in_window = 1000; # or something big
							
							my $window_num_defined_codons_g1 = 0;
							for(my $j = $first_codon; $j <= $last_codon; $j++) {
								$window_num_defined_codons_g1 += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$j}->{'num_defined_codons_g1'};
								
								if($between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$j}->{'num_defined_codons_g1'} < $min_defined_codons_in_window) {
									$min_defined_codons_in_window = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$j}->{'num_defined_codons_g1'};
								}
							}
							
							my $window_num_defined_codons_g2 = 0;
							for(my $j = $first_codon; $j <= $last_codon; $j++) {
								$window_num_defined_codons_g2 += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$j}->{'num_defined_codons_g2'};
								
								if($between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$j}->{'num_defined_codons_g2'} < $min_defined_codons_in_window) {
									$min_defined_codons_in_window = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$j}->{'num_defined_codons_g2'};
								}
							}
							
							if($min_defined_codons_in_window==0 || $min_defined_codons_in_window eq '' || ! defined($min_defined_codons_in_window)) {
								$min_defined_codons_in_window = 0;
							}
							
							my $window_mean_num_defined_codons = ($window_num_defined_codons_g1 + $window_num_defined_codons_g2)/2;
							
							# In case filters not passed by ANY codons therein, all results blank
							if(! $window_N_sites && ! $window_S_sites && ! $window_N_diffs && ! $window_S_diffs) {
								$window_N_sites = 'NA';
								$window_S_sites = 'NA';
								$window_N_diffs = 'NA';
								$window_S_diffs = 'NA';
								
								$window_dN_exceeds_product_dS = 'NA';
							}
							
							my $out_line_between_sw = "$this_analysis\t$this_product\t$group_name_i\t$group_name_j\t".
								"$window_num\t$first_codon\t$last_codon\t" .
								"$window_num_defined_codons_g1\t$window_num_defined_codons_g2\t$window_mean_num_defined_codons\t$min_defined_codons_in_window\t".
								"$window_N_sites\t$window_S_sites\t$window_N_diffs\t$window_S_diffs\t".
								"$window_dN\t$window_dS\t"; # $window_dN_minus_dS\t$window_w\t$window_dN_exceeds_product_dS";
							
							
						
							if($num_bootstraps > 1) {
								$out_line_between_sw .= "$window_SE_dN\t$window_SE_dS\t"
							}
							
							$out_line_between_sw .= "$window_dN_minus_dS\t$window_w";
							
							if($num_bootstraps > 1) {
								my $significance = '';
								
								if($window_SE_dN_minus_dS > 0) {
									
									if(abs($window_boot_Z) > 2.81) {
										$significance = '***';
									} elsif(abs($window_boot_Z) > 1.96) {
										$significance = '**';
									} elsif(abs($window_boot_Z) > 1.64) {
										$significance = '*';
									}
								
								} else {
									$significance = 'NA';
								}
								
								# Z-tests worthless if only one codon in window
								if($sliding_window_size == 1) {
									$window_SE_dN_minus_dS = 'NA';
									$window_boot_Z = 'NA';
									$significance = 'NA';
								}
								
								$out_line_between_sw .= "\t$window_SE_dN_minus_dS\t$window_boot_Z\t$significance\n";
								
							} else { # no bootstraps
								$out_line_between_sw .= "\n";
							}
							
							
							open(OUTFILE_BETWEEN_SW,">>between\_group\_sw_" . $sliding_window_size . "codons\_results.txt");
							print OUTFILE_BETWEEN_SW "$out_line_between_sw";
							close OUTFILE_BETWEEN_SW;
											
						} # we're in range of the product
					} # finish sliding window step-through
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					#############################
					# BOOTSTRAP THE WHOLE PRODUCT
					my $SE_dN_minus_dS = 0;
					my $product_boot_Z = 'NA';
					
					my $SE_dN;
					my $SE_dS;
					
					my $SE_JC_dN_minus_dS = 0;
					my $JC_product_boot_Z = 'NA';
					
					my $SE_JC_dN;
					my $SE_JC_dS;
					
					if($num_bootstraps > 1) {
						# MAKE BOOTSTRAP FILE
						open(BOOT_FILE_PRODUCT,">>$this_product\_bootstrap\_results\.txt");
						print BOOT_FILE_PRODUCT "analysis\tgroup1\tgroup2\tproduct_name\tbootstrap_num\tsim_product_N_sites_sum\t".
								"sim_product_S_sites_sum\tsim_product_N_diffs_sum\tsim_product_S_diffs_sum\n";
						
						#################
						# BOOTSTRAPPING TO CALCULATE STANDARD ERROR HERE?
						print "\n#########################################################################################\n" .
							"Bootstrapping $this_product between-group $group_name_i vs. $group_name_j for dN/dS standard error...\n";
						
						#print "num_codons_in_product: $num_codons_in_product\n";
						
						my @sim_N_sites_arr;
						my @sim_S_sites_arr;
						my @sim_N_diffs_arr;
						my @sim_S_diffs_arr;
						
				##		# PARALLELIZE IT
				##		mkdir("$product_name\_bootstrap_temp_files");
				##		chdir("$product_name\_bootstrap_temp_files");
				##		#my $procs = 60; # number to do at once, cores to assign simultaneously, machine-limited. CUVIER has 80
				##		my $pm = Parallel::ForkManager->new($procs_per_node);
							
						for(my $bootstrap_num = 1; $bootstrap_num <= $num_bootstraps; $bootstrap_num++) {
				##			$pm->start and next; # this is IT
							srand();
							
							# New bootstrap run, of $num_bootstraps
							#print "bootstrap $bootstrap_num\n";
							
							my $sim_product_N_sites_sum = 0;
							my $sim_product_S_sites_sum = 0;
							my $sim_product_N_diffs_sum = 0;
							my $sim_product_S_diffs_sum = 0;
							
							#my $random_sum = 0;
							
							my %sim_poly_site_data_hh; # $sim_poly_site_data_hh{codon_index}->{unique_variants}
							
							# SO, HERE'S the problem: we shouldn't sample WITHIN SITES, but rather SAMPLE
							# THE SITES THEMSELVES!!
							
							# (1) indeed, go through $num_codons_in_product
							# (2) for each, selection a random codon position's results: Nd, Sd, Ns, Ss. 
							# where are the aforementioned stored? We've just added the following to do so: 
							
				#				$codon_to_results_hh{$codon_num}->{N_diffs}
				#				$codon_to_results_hh{$codon_num}->{S_diffs}
				#				$codon_to_results_hh{$codon_num}->{N_sites}
				#				$codon_to_results_hh{$codon_num}->{S_sites}
							
							# SAMPLE codon sites, up to the actual number of codons
							for(my $codon_index = 0; $codon_index < $num_sorted_codons; $codon_index++) {
								# Which codon site do we choose?
								my $random_codon_num = int(rand($num_sorted_codons + 1));
								
								my $num_defined_codons_g1 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{'num_defined_codons_g1'};
								my $num_defined_codons_g2 = $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{'num_defined_codons_g2'};
								
								if($codon_min_taxa_total > 0 && $codon_min_taxa_group > 0) { # BOTH criteria have been specified
									my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
									my $taxa_group_min = $num_defined_codons_g1;
									
									if($num_defined_codons_g2 < $taxa_group_min) {
										$taxa_group_min = $num_defined_codons_g2;
									}
									
									if($taxa_total >= $codon_min_taxa_total && $taxa_group_min >= $codon_min_taxa_group) { # see if criteria are met by this codon
										$sim_product_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
										$sim_product_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
										$sim_product_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
										$sim_product_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
									}
								} elsif($codon_min_taxa_total > 0) { # only TOTAL min has been specified
									my $taxa_total = $num_defined_codons_g1 + $num_defined_codons_g2;
									
									if($taxa_total >= $codon_min_taxa_total) { # see if criterion is met by this codon
										$sim_product_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
										$sim_product_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
										$sim_product_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
										$sim_product_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
									}
								} elsif($codon_min_taxa_group > 0) { # only GROUP min has been specified
									my $taxa_group_min = $num_defined_codons_g1;
									
									if($num_defined_codons_g2 < $taxa_group_min) {
										$taxa_group_min = $num_defined_codons_g2;
									}
									
									if($taxa_group_min >= $codon_min_taxa_group) { # see if criterion is met by this codon
										$sim_product_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
										$sim_product_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
										$sim_product_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
										$sim_product_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
									}
								} else { # BOTH CRITERIA 0 OR UNSPECIFIED
									$sim_product_N_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_diffs};
									$sim_product_S_diffs_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_diffs};
									$sim_product_N_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{N_sites};
									$sim_product_S_sites_sum += $between_group_codon_data_hh{$this_analysis}->{$group_name_i}->{$group_name_j}->{$this_product}->{$random_codon_num}->{S_sites};
								}
								
							} # finished compiling all sampled codons (cols in alignment)
							
							# Print product totals for BOOTSTRAP
							my $out_line_boot = "ANALYSIS\tgroup_name_i\t$group_name_j\t$this_product\t$bootstrap_num\t$sim_product_N_sites_sum\t".
								"$sim_product_S_sites_sum\t$sim_product_N_diffs_sum\t$sim_product_S_diffs_sum\n";
							
							print BOOT_FILE_PRODUCT "$out_line_boot";
							
							#print "$out_line_boot";
							
							push(@sim_N_diffs_arr,$sim_product_N_diffs_sum);
							push(@sim_S_diffs_arr,$sim_product_S_diffs_sum);
							push(@sim_N_sites_arr,$sim_product_N_sites_sum);
							push(@sim_S_sites_arr,$sim_product_S_sites_sum);
							
							#print $random_sum . "\n";
							
						} # end last bootstrap
						
						close BOOT_FILE_PRODUCT;
						
						# Calculate ACTUAL overall values for the product
#						my $between_product_dN_filtered;
#						if($between_product_N_sites_sum_filtered > 0) {
#							$between_product_dN_filtered = $between_product_N_diffs_sum_filtered / $between_product_N_sites_sum_filtered;
#						} else {
#							$between_product_dN_filtered = '*';
#						}
#						
#						my $between_product_dS_filtered;
#						if($between_product_S_sites_sum_filtered > 0) {
#							$between_product_dS_filtered = $between_product_S_diffs_sum_filtered / $between_product_S_sites_sum_filtered;
#						} else {
#							$between_product_dS_filtered = '*';
#						}
						
						my $between_product_w_filtered;
						if ($between_product_dS_filtered > 0 && $between_product_dS_filtered ne '*') {
							$between_product_w_filtered = $between_product_dN_filtered / $between_product_dS_filtered;
						} else {
							$between_product_w_filtered = '*';
						}
						
						my $between_product_dN_minus_dS_filtered;
						if($between_product_dN_filtered >= 0 && $between_product_dS_filtered >= 0) {
							$between_product_dN_minus_dS_filtered = $between_product_dN_filtered - $between_product_dS_filtered;
						}
						
						my $JC_w;
						if ($between_JC_dS > 0 && $between_JC_dS ne '*') {
							$JC_w = $between_JC_dN / $between_JC_dS;
						} else {
							$JC_w = '*';
						}
						
						my $JC_dN_minus_dS;
						if($between_JC_dN >= 0 && $between_JC_dS >= 0) {
							$JC_dN_minus_dS = $between_JC_dN - $between_JC_dS;
						}
						
						# CALCULATE BOOTSTRAP STANDARD ERROR HERE! NEI & KUMAR (2000) EQUATIONS
						my @sim_dN_minus_dS;
						my @sim_dN;
						my @sim_dS;
						
						my @sim_JC_dN_minus_dS;
						my @sim_JC_dN;
						my @sim_JC_dS;
						
						for(my $sim_num = 0; $sim_num < scalar(@sim_N_sites_arr); $sim_num++) {
							my $this_round_dN = '*';
							if($sim_N_sites_arr[$sim_num] > 0) {
								$this_round_dN = $sim_N_diffs_arr[$sim_num] / $sim_N_sites_arr[$sim_num];
							}
							
							my $this_round_JC_dN = -3/4 * log(1 - ((4/3) * $this_round_dN));
							
							my $this_round_dS = '*';
							if($sim_S_sites_arr[$sim_num] > 0) {
								$this_round_dS = $sim_S_diffs_arr[$sim_num] / $sim_S_sites_arr[$sim_num];
							}
							
							my $this_round_JC_dS = -3/4 * log(1 - ((4/3) * $this_round_dS));
							
							my $this_round_dN_minus_dS = '*';
							if($this_round_dN >= 0 && $this_round_dS >= 0) {
								$this_round_dN_minus_dS = $this_round_dN - $this_round_dS;
							}
							
							my $this_round_JC_dN_minus_dS = '*';
							if($this_round_JC_dN >= 0 && $this_round_JC_dS >= 0) {
								$this_round_JC_dN_minus_dS = $this_round_JC_dN - $this_round_JC_dS;
							}
							
							push(@sim_dN_minus_dS, $this_round_dN_minus_dS);
							push(@sim_dN, $this_round_dN);
							push(@sim_dS, $this_round_dS);
							
							push(@sim_JC_dN_minus_dS, $this_round_JC_dN_minus_dS);
							push(@sim_JC_dN, $this_round_JC_dN);
							push(@sim_JC_dS, $this_round_JC_dS);
						}
						
						$SE_dN_minus_dS = &standard_deviation(@sim_dN_minus_dS);
						$SE_dN = &standard_deviation(@sim_dN);
						$SE_dS = &standard_deviation(@sim_dS);
						
						$SE_JC_dN_minus_dS = &standard_deviation(@sim_JC_dN_minus_dS);
						$SE_JC_dN = &standard_deviation(@sim_JC_dN);
						$SE_JC_dS = &standard_deviation(@sim_JC_dS);
						
						if($SE_dN_minus_dS > 0) {
							$product_boot_Z = $between_product_dN_minus_dS_filtered / $SE_dN_minus_dS;
						}
						
						if($SE_dN_minus_dS == 0) {
							$SE_dN_minus_dS = 'NA';
						}
						
						if($SE_JC_dN_minus_dS > 0) {
							$JC_product_boot_Z = $JC_dN_minus_dS / $SE_JC_dN_minus_dS; #?#
						}
						
						if($SE_JC_dN_minus_dS == 0) {
							$SE_JC_dN_minus_dS = 'NA';
						}
						
				#		print "\nproduct SE(dN-dS) = $SE_dN_minus_dS\n".
				#			"product Z-value = $product_boot_Z\n\n";
						
					} # END BOOTSTRAPS
					
					
					# Print product totals
					my $product_dN;
					my $product_dS;
					my $product_dN_over_dS;
					
					if($between_product_N_sites_sum_filtered > 0) {
						$product_dN = $between_product_N_diffs_sum_filtered / $between_product_N_sites_sum_filtered;
					} else {
						$product_dN = '*';
					}
					
					if($between_product_S_sites_sum_filtered > 0) {
						$product_dS = $between_product_S_diffs_sum_filtered / $between_product_S_sites_sum_filtered;
					} else {
						$product_dS = '*';
					}
					
					my $product_dN_minus_dS = $product_dN - $product_dS;
					
					if($product_dS > 0) {
						$product_dN_over_dS = $product_dN / $product_dS;
					} else {
						$product_dN_over_dS = '*';
					}
					
					
					# Jukes-Cantor Correction
					my $product_JC_dN_minus_dS = $between_JC_dN - $between_JC_dS;
					
					my $product_JC_w;
					if($between_JC_dS > 0) {
						$product_JC_w = $between_JC_dN / $between_JC_dS;
					} else {
						$product_JC_w = '*';
					}
					
					
					my $out_line = "$this_analysis\t$this_product\t$group_name_i\t$group_name_j\t$between_product_N_sites_sum_filtered\t".
						"$between_product_S_sites_sum_filtered\t$between_product_N_diffs_sum_filtered\t$between_product_S_diffs_sum_filtered\t".
						"$between_product_dN_filtered\t$between_product_dS_filtered\t" .
						"$between_JC_dN\t$between_JC_dS\t";
						
					if($num_bootstraps > 1) {
						$out_line .= "$SE_dN\t$SE_dS\t$SE_JC_dN\t$SE_JC_dS\t";
					}
					
					$out_line .= "$product_dN_minus_dS\t$product_dN_over_dS" .
						"\t$product_JC_dN_minus_dS\t$product_JC_w";
					
					my $rounded_dN = sprintf("%.5f",$product_dN);
					my $rounded_dS = sprintf("%.5f",$product_dS);
					my $rounded_dNdS = sprintf("%.5f",$product_dN_over_dS);
					
					if($product_dN_over_dS eq '*') {
						$rounded_dNdS = '*';
					}
					
					print "\ndN=$rounded_dN\ndS=$rounded_dS\ndN\/dS=$rounded_dNdS\n";
					
					if($num_bootstraps > 1) {
				#		my $z_value;
						my $rounded_SE = sprintf("%.5f",$SE_dN_minus_dS);
						my $rounded_Z = sprintf("%.5f",$product_boot_Z);
						
						print "\nSE(dN-dS) = $rounded_SE\n".
							"Z-value = $rounded_Z\n\n";
						
						my $significance = '';
						
						if($SE_dN_minus_dS > 0) {
						
				#			$z_value = $product_dN_minus_dS / $SE_dN_minus_dS; 
							
							if(abs($product_boot_Z) > 2.81) {
								$significance = '***';
							} elsif(abs($product_boot_Z) > 1.96) {
								$significance = '**';
							} elsif(abs($product_boot_Z) > 1.64) {
								$significance = '*';
							}
						
						}
						
						my $JC_significance = '';
						
						if($SE_JC_dN_minus_dS > 0) {
						
				#			$z_value = $product_dN_minus_dS / $SE_dN_minus_dS; 
							
							if(abs($JC_product_boot_Z) > 2.81) {
								$JC_significance = '***';
							} elsif(abs($JC_product_boot_Z) > 1.96) {
								$JC_significance = '**';
							} elsif(abs($JC_product_boot_Z) > 1.64) {
								$JC_significance = '*';
							}
						
						}
						
						$out_line .= "\t$SE_dN_minus_dS\t$product_boot_Z\t$significance" . 
							"\t$SE_JC_dN_minus_dS\t$JC_product_boot_Z\t$JC_significance\n";
						
						
					} else {
						$out_line .= "\n";
					}
					
					# PRINT OUT RE-CALCULATED PARTITION VALUES
					open(BETWEEN_OUTFILE,">>$between_group_product_results_file");
		
					print BETWEEN_OUTFILE "$out_line";
					#print "$out_line";
		
					close BETWEEN_OUTFILE;
					
				} # done looping products in the family
			} # done looping group 2
		} # done looping group 1
		
	} # done looping analyses in the file
	
} # finish case in which between_group_codon_results.txt exists
	

print "\n\n";


#########################################################################################
# END PROGRAM and print a completion message to screen
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
sub mean {
    my @values = @_;
    my $length = @values;
    my $sum;
    
    foreach (@values) {
    	$sum += $_;
    }
    
    return($sum / $length);
}


#########################################################################################
sub standard_deviation {
    my @values = @_;
    my $length = @values;
	my $mean_of_values = &mean(@values);
	my $sum_squared_deviations;
	
	foreach (@values) {
		$sum_squared_deviations += ($_ - $mean_of_values)**2;
	}
	
	my $variance = ($sum_squared_deviations) / ($length - 1);
	
	return(sqrt($variance));
	
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
				@line_arr = split(/\t/,$_);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,\w+,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_);
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
				@line_arr = split(/\t/,$_);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_);
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
	
	print "SNPGenie completed at local time $local_time2. The process took $time_diff_rounded secs, i.e., ".
			"$whole_mins_elapsed mins and $secs_remaining_rounded secs\n";

	print "\n################################################################################".
		"\n##          SNPGenie between-group processor completed successfully.          ##".
		"\n##        Please find results in the \/data\/ dir and subdirs (families).       ##\n".
		"################################################################################".
		"\n\n\n"; 
}
