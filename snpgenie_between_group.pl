#!/usr/bin/env perl

# PROGRAM: SNPGenie, between-group analysis, for two or more aligned multi-FASTA files.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# snpgenie_between_group.pl --gtf_file=<CDS_annotations>.gtf --num_bootstraps=10000 --procs_per_node=16
#########################################################################################

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

# DATE CREATED: April 10, 2015
# DATE MODIFIED: December 16, 2016
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com
# AFFILIATION1: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA
# AFFILIATION2: BigPlant Consortium, Center for Genomics and Systems Biology, New York 
#     University, New York, NY 10003, USA

# ADD:
# —SORT codons before entry into the NNN-NNN comps hh, so we don't have, e.g., both AAG-AAA and AAA-AAG
# —ADD MAJORITY CODON to codon file
# —PROBLEM with min-codon-count in variant results

use strict;
#use warnings;
use Data::Dumper;
use List::Util qw(max);
use Parallel::ForkManager;
use Getopt::Long;

my $time1 = time;
my $local_time1 = localtime;

#########################################################################################
# INITIALIZE (OPTIONAL) INPUT VARIABLES

# AUTOMATICALLY DETECT FASTA FILES IN WORKING DIRECTORY
my @fasta_files = &get_fasta_file_names;
foreach(@fasta_files) {
	unless($_ =~ /.fa/) { die "\n\n# FASTA files must contain .fa, .fas, or .fasta extension. TERMINATED\n\n"; }
}

# Initialize others
my $gtf_file;
my $procs_per_node;
my $num_bootstraps;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "gtf_file=s" => \$gtf_file,
			"procs_per_node=i" => \$procs_per_node,
			"num_bootstraps=i" => \$num_bootstraps)
			
			or die "\n### WARNING: Error in command line arguments. Script terminated.\n\n";
			# If an argument is called as a flag, its value is 0; if not called, it's null


unless($gtf_file =~ /.gtf/) {
	die "\n### WARNING: The --gtf_file option must be a file with a .gtf or .gtf extension\n".
		"### SNPGenie terminated.\n\n";
}

if(! $procs_per_node) {
	if($procs_per_node != 0) {
		$procs_per_node = 1; # DEFAULT is 1
	}
} elsif($procs_per_node < 1) {
	die "\n### WARNING: The --procs_per_node option must be an integer ≥1\n".
		"### SNPGenie terminated.\n\n";
} 

if(! $num_bootstraps) {
	if($num_bootstraps != 0) { # Called as a flag, but given no value
		$num_bootstraps = 1000; # default behavior: 1000 bootstraps
	}
} else {
	print "\n### No bootstrapping will occur.\n\n";
}

print "\n# Number of parallel processes to be invoked is $procs_per_node\.\n";

print "\n################################################################################".
	"\n##                                                                            ##".
	"\n##                      Between-Group SNPGenie Initiated!                     ##".
	"\n##                                                                            ##".
	"\n################################################################################\n";


print "\nSNPGenie initiated at local time $local_time1\n";

#########################################################################################
# Store product information
my @product_names_arr = &get_product_names_from_gtf($gtf_file);
@product_names_arr = sort {$a <=> $b} @product_names_arr;
#print "\nProduct names: @product_names_arr\n";

# Determine all product coordinates
my %product_coordinates_ha;
foreach my $product_name (@product_names_arr) {
	#print "product name is: $product_name\n";
	my @product_coord_arr = &get_product_coordinates($product_name);
	#print "\n\n$product_name product_coord arr: @product_coord_arr\n\n";
	
	# Save to a hash of arrays
	push(@{$product_coordinates_ha{$product_name}->{product_coord_arr}},@product_coord_arr);
	#print "\n\n$product_name product_coord arr in harr: @{$product_coordinates_ha{$product_name}->{product_coord_arr}} \n\n";
}

my @seqs_aa;

#print "\nfasta file 0 (now b) is: $fasta_files[0]\n";

my @fasta_file_names;

# This proves too prohibitive for big sequence data
# Instead, store just for the products and groups individually, below
foreach(@fasta_files) { # GROUP NAMES ARE THE FASTA FILES
	# Read in the sequence from the file
	my $seq = '';
	my @seqs_arr;
	#my $header = '';
	#my @headers_arr;
	my $seq_num = 0;
	my $last_seq_length;
	
	open(IN, "$_") or die "Could not open file $_\n";
	
	print "\nRecording coding sequence data for $_...\n";
	
	push(@fasta_file_names,$_);
	
	while(<IN>) {
		chomp;
		if(/>/) {
			if($seq_num == 0) {
				#$header = $_;
				$seq_num ++;
			} else {
				$seq = uc($seq);
				$seq =~ tr/U/T/;
				
				push(@seqs_arr,$seq);
				#push(@headers_arr,$header);
				#$header = $_;
				$seq_num ++;
				
				my $this_seq_length = length($seq);
				
				#print "\nseq $seq_num is of length $this_seq_length\n";
				
				if($last_seq_length && ($last_seq_length != $this_seq_length)) {
					die "\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
				} else {
					$last_seq_length = $this_seq_length;
					#print "\nseq: $seq\n";
					$seq = '';
				}
			}
		} else {
			$seq .= $_;
		}
	}
	
	close IN;
	
	$seq = uc($seq);
	$seq =~ tr/U/T/;
	push(@seqs_arr,$seq);
	#push(@headers_arr,$header);
	push(@seqs_aa,\@seqs_arr); # ??? work this way. Yes, seems to
}


print "\n";

#print "\nfasta file 0 (now c) is: $fasta_files[0]\n";

# Go through each product and DO DIS THANG
# HERE'S HOW to get the first array:
#print "@{$seqs_aa[0]}\n";

#my $counter=1;
#foreach(@{$seqs_aa[0]}) {
#	print "$counter: $_\n\n";
#	$counter++;
#}

my %products_groups_seqs_hh;
# FORMAT:
#%products_groups_seqs_hh = {
#	'ORF1a' => {
#		'group_1' => {
#			'seq_1' => 'ACGT',
#			'seq_2' => 'ACGT',
#			... ,
#			'seq_n' => 'ACGT',
#		},
#		'group_2' => {
#			'seq_1' => 'ACGT',
#			'seq_2' => 'ACGT',
#			... ,
#			'seq_n' => 'ACGT',
#		}
#	},
#	
#	'ORF1b' => {
#		'group_1' => {
#			'seq_1' => 'ACGT',
#			'seq_2' => 'ACGT',
#			... ,
#			'seq_n' => 'ACGT',
#		},
#		'group_2' => {
#			'seq_1' => 'ACGT',
#			'seq_2' => 'ACGT',
#			... ,
#			'seq_n' => 'ACGT',
#		}
#	}
#};

foreach(@product_names_arr) {
	my $product_name = $_;
	my @product_coord_arr = @{$product_coordinates_ha{$product_name}->{product_coord_arr}};
	
	# New segments approach
	my %product_starts;
	my %product_stops;
	
	my $num_segments = (@product_coord_arr / 2);
	for(my $i=1; $i<=$num_segments; $i++) { # $i<=scalar(@product_coord_arr)
		$product_starts{$i} = $product_coord_arr[2*$i-2];
		$product_stops{$i} = $product_coord_arr[2*$i-1];
	}
	
	# Build an array of all the sequences for each group
	my $group_num = 0;
	foreach (@seqs_aa) { # for each group
		my @group = @{$_};
		
		
		
		my $group_id = $fasta_file_names[$group_num];
		
		
		$group_num ++;
#		my $group_id = 'group_' . $group_num;
		
		
		
		
		
		
		my $seq_num = 0;
		foreach my $sequence (@group) { # for each seq in alignment (max coverage)
			$seq_num ++;
			my $seq_id = 'seq_' . $seq_num;
			
			#build the product seq and add it to full group list
			my $this_product_seq = '';
			
			for(my $seg=1; $seg<=$num_segments; $seg++) { # for each segment
				#my $this_index = $seg - 1;
				
				#my $this_start = $product_starts[$this_index];
				my $this_start = $product_starts{$seg};
				my $this_start_index = $this_start - 1;
				
				#my $this_stop = $product_stop[$this_index];
				my $this_stop = $product_stops{$seg};
				my $this_stop_index = $this_stop - 1;
				
				my $this_seg_length = ($this_stop - $this_start + 1);
				
				$this_product_seq .= substr($sequence,$this_start_index,$this_seg_length);
			} # done building sequence with all segments	
			
			# Add the sequence to an array of all sequences for this product in the group
			$products_groups_seqs_hh{$product_name}->{$group_id}->{$seq_id} = $this_product_seq;
		}
	}
}

#my @test_keys = keys %{$products_groups_seqs_hh{'ORF5'}->{'group_1'}};
#@test_keys = keys %products_groups_seqs_hh;
#print "\n\n @test_keys \n\n";
#my $test_seq = $products_groups_seqs_hh{'ORF5'}->{'group_1'}->{'seq_1'};
#print "\n\n $test_seq \n\n";

STDOUT->autoflush(1);


# READY OUTPUT FILES
open(CODON_FILE,">>between\_group\_codon\_results\.txt");
print CODON_FILE "analysis\tproduct\tgroup_1\tgroup_2\tcodon\tvariability\t".
	"num_defined_codons_g1\tnum_defined_codons_g2\tcomparisons\tN_sites\tS_sites\tN_diffs\tS_diffs\n";
close CODON_FILE;

open(OUTFILE,">>between\_group\_product\_results\.txt");
print OUTFILE "analysis\tproduct\tgroup_1\tgroup_2\tN_sites\tS_sites\tN_diffs\tS_diffs\t".
	"dN\tdS\t".
	"dN\-dS\t".
	"dN\/dS\n";
#print "group_1\tgroup_2\tproduct\tN_sites\tS_sites\tN_diffs\tS_diffs\n";

#if($num_bootstraps > 1) {
#	print OUTFILE "\tSE_piN_minus_piS\tz_value\n";
#} else {
#	print OUTFILE "\n";
#}

# YES, it all works
#my %product_data_hh;
my $between_total_N_sites = 0;
my $between_total_S_sites = 0;
my $between_total_N_diffs = 0;
my $between_total_S_diffs = 0;

my @group_names_arr;

foreach my $product_name (keys %products_groups_seqs_hh) { # for each product
#foreach my $product_name (my @thisarray = ('ORF5')) { # for each product
	print "\n#########################################################\n".
		"################ Analyzing $product_name\n".
		"#########################################################\n";

	# I want, for this product:
	#groups_codons_aa[group_index]->[codon_index]->array of all codons in group
	my @groups_codons_aa;
	
	my $last_num_codons;
	my $group_index = 0;
	
	# Store every sequence for each group; this goes really fast
	foreach my $group_name (keys %{$products_groups_seqs_hh{$product_name}}) { # for each group (fasta file)
		print "\ngroup index $group_index is group $group_name\n";
		
		push(@group_names_arr, $group_name);
		
		foreach my $sequence_name (keys %{$products_groups_seqs_hh{$product_name}->{$group_name}}) { # for each seq in group
			my $sequence = $products_groups_seqs_hh{$product_name}->{$group_name}->{$sequence_name};
			
			# Make sure it's a complete set of codons
			my $seq_length = length($sequence);
			if(($seq_length % 3) != 0) {
				die "\n\nDIE: A sequence in $product_name is not a multiple of 3 (complete codon set). TERMINATED.\n\n";
			}
			
			# Build all codons, put in @aa with @codon_num->@codons
			my $num_codons = $seq_length / 3;
			
			if((! $last_num_codons) && ($last_num_codons ne '') && ($num_codons != $last_num_codons)) {
				die "\n\nDIE: In $product_name, there are sequences of different length ($last_num_codons and $num_codons). TERMINATED.\n\n";
			} else {
				$last_num_codons = $num_codons;
				#print "changed\n";
			}
			
			my $codon_index = 0; # so we're going to populate $groups_codons_aa[0]->[0]->@
			for(my $i=1; $i<=$num_codons; $i++) { # for each codon
				my $array_index = $i - 1;
				#my $codon = substr($sequence,$codon_index,3);
				# More efficient to call substr at beginning and gobble
				my $codon = substr($sequence,0,3,"");
				
				push(@{$groups_codons_aa[$group_index]->[$array_index]},$codon);
				
				#$codon_index+=3;
			}
		}
		$group_index++;
	} # Done storing all groups
	
	# Compare all codons in one group to all codons in another
	my $num_groups = scalar @groups_codons_aa;
#	my $num_codons_1 = scalar @{$groups_codons_aa[0]};
#	my $num_seqs_1 = scalar @{$groups_codons_aa[0]->[0]};
#	my $num_codons_2 = scalar @{$groups_codons_aa[1]};
#	my $num_seqs_2 = scalar @{$groups_codons_aa[1]->[0]};
#	
#	print "\n\n num groups is $num_groups \n";
#	print " num codons for $product_name in group 1 $num_codons_1 \n";
#	print " num seq for $product_name in group 1 $num_seqs_1 \n";
#	print " num codons for $product_name in group 2 $num_codons_2 \n";
#	print " num seq for $product_name in group 2 $num_seqs_2 \n\n";
	
	# COMEBACK: make a number of sites subroutine for the whole codon, not just the sites
	# OR, instead of calling the subroutine every time, just COUNT the comparisons!
	# $hh_comps{CGA}->{CGG}->256 YES!!
	
	# COMPARE ALL THIS POSITION'S CODONS IN ONE GROUP TO ALL CODONS IN ANOTHER
	
	# GROUP i
	for(my $i=0; $i<$num_groups; $i++) {
		
		
		
		# STORE THE SEQUENCES FOR THIS PRODUCT FROM GROUP I ONLY
		# the fasta file corresponding to group i is:
		# push(@{$groups_codons_aa[$group_index]->[$array_index]},$codon);
		# foreach my $group_name (keys %{$products_groups_seqs_hh{$product_name}}) { # for each group (fasta file) # <- order in which stored
		
		
		
		
		# GROUP j
		for(my $j=$i+1; $j<$num_groups; $j++) { # comparing groups i and j
			
			# Get species/sequence names for each group.
			#?my @group_1_spp;
			#?my @group_2_spp;
			#?push(@group_1_spp,$_);
			
			my $group_name_i = $group_names_arr[$i];
			my $group_name_j = $group_names_arr[$j];
			
			#print "\nmy group name i is: " . $fasta_files[0] . "\n";
			#print "\nmy group name i is: $group_name_i\nmy group name j is: $group_name_j\n";
			
			
			# STORE THE SEQUENCES FOR THIS PRODUCT FROM GROUP J ONLY
			# the fasta file corresponding to group j is:
			
			
			my $num_codons_in_product = scalar @{$groups_codons_aa[$i]};
			#my $num_seqs_per_codon = scalar @{$groups_codons_aa[$i]->[0]}; # but this is different for the two groups COMEBACK#####
			# THIS NEEDS TO BE DIFFERENT FOR TWO GROUPS
			
			# Updated here:
			my $num_seqs_per_codon_groupi = scalar @{$groups_codons_aa[$i]->[0]};
			my $num_seqs_per_codon_groupj = scalar @{$groups_codons_aa[$j]->[0]};
			
			print "\nComparing $num_codons_in_product codons in $product_name between groups $i and $j\...\n";
			
			#groups_codons_aa[group_index]->[codon_index]->array of all codons in group
			
			
			
			
			
			
			
			
			# Quickly see if it's polymorphic between-group to save time. It's faster just to check,
			# EVEN IF not polymorphic
			my @polymorphic_codons_arr;
			my %codonIndex_conserved;
			
			print "\nDetermining which codons are polymorphic...\n";
			
			
			
			
			
			
			# PARALLELIZE IT
##
			mkdir("$product_name\_polymorphic_codons_arr");
			chdir("$product_name\_polymorphic_codons_arr");
#			my $procs_poly = 10; # number to do at once, cores to assign simultaneously, machine-limited. CUVIER has 80
			my $pm_poly = Parallel::ForkManager->new($procs_per_node);


			# DETERMINE WHICH CODONS ARE POLYMORPHIC
			for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) { # for each codon in product
##
				$pm_poly->start and next; # this is IT
				
				#print "codon index $codon_index\n";
				#my $codon_to_print = $groups_codons_aa[$i]->[$codon_index]->[1];
				#$| = 1;
				#print "$codon_to_print";
				my $between_g_comparisons = 0;
				my $this_codon_poly = 0;
				
				OUTER: for (my $gi_seq_index = 0; $gi_seq_index < $num_seqs_per_codon_groupi; $gi_seq_index++) { # for each gi sequence at this codon
					
					# Codon to compare against
					
					INNER: for (my $gj_seq_index = 0; $gj_seq_index < $num_seqs_per_codon_groupj; $gj_seq_index++) { # for each gj sequence at this codon
						#?my $group_1_species = $group_1_spp[$gi_spp_index];
						#?my $group_2_species = $group_2_spp[$gj_spp_index];
						
						my $codon_gi = $groups_codons_aa[$i]->[$codon_index]->[$gi_seq_index];
						my $codon_gj = $groups_codons_aa[$j]->[$codon_index]->[$gj_seq_index];
						#print "Comparing codons $codon_gi and $codon_gj\n";
						
						#?my $codon_gi = $data_codonnum_spp_codon_hh{$codon_num}->{$group_1_species};
						#?my $codon_gj = $data_codonnum_spp_codon_hh{$codon_num}->{$group_2_species};
						
						if(($codon_gi ne $codon_gj) && ! ($codon_gi =~ 'N' || $codon_gi =~ '-' || $codon_gi eq '' || $codon_gj =~ 'N' || $codon_gj =~ '-' || $codon_gj eq '')) {
							$this_codon_poly = 1;
							#print "\ncodon $codon_index is polymorphic\n";
							last OUTER; # gotta break of out TWO loops here and go to next codon in product
						}
					} # INNER
				} # OUTER
##				push(@polymorphic_codons_arr,$this_codon_poly);	
##
				open(THIS_POLY_TEMP_FILE,">>$codon_index");
				print THIS_POLY_TEMP_FILE $this_codon_poly;
				close THIS_POLY_TEMP_FILE;
				
				$pm_poly->finish; # special name
			}

##				
			$pm_poly->wait_all_children; # special name, methods within module
			
			my @polymorphic_codons_FILES_arr = glob "*";
			@polymorphic_codons_FILES_arr = sort {$a <=> $b} @polymorphic_codons_FILES_arr;
			
			#print "sorted polymorphic_codons_FILE_arr: @polymorphic_codons_FILES_arr\n";
			
			foreach(@polymorphic_codons_FILES_arr) { # file names, actually
				my $poly_file = $_;
				
				open(CURR_POLY_FILE, $poly_file) or die "\n## Cannot open $poly_file. TERMINATED.\n\n";
				while(<CURR_POLY_FILE>) {
					chomp;
					push(@polymorphic_codons_arr,$_);	
				}
				close CURR_POLY_FILE;
				unlink $poly_file;
			}
			
			chdir("..");
			rmdir("$product_name\_polymorphic_codons_arr");
				
				
			
			open(POLY_FILE,">>between\_group\_polymorphism.txt");
			print POLY_FILE "\nDone recording polymorphism in $product_name between groups $i and $j\nIt is: @polymorphic_codons_arr\n\n";
			close POLY_FILE;
			
			print "\nDone recording polymorphism in $product_name between groups $i and $j\n";
			#?print "\nDone recording polymorphism in partition $partition_id between the groups:\nGROUP 1: $group_name_i\nGROUP 2: $group_name_j\n";
			#print "It is: @polymorphic_codons_arr\n\n";
			
			print "\nAnalyzing polymorphic codons in $product_name\...\n";
			
			# calculate number of total codons so that we can track % progress?
			# INTRO: skip if not poly, just do values for the first codon observed
			
#			my $between_product_N_sites_sum = 0;
#			my $between_product_S_sites_sum = 0;
#			my $between_product_N_diffs_sum = 0;
#			my $between_product_S_diffs_sum = 0;
			
			
			print "\nAnalyzing each codon...\n";
			
			
			
			
			# PARALLEL
			mkdir("$product_name\_codon\_analysis");
			chdir("$product_name\_codon\_analysis");
			my $pm_codons = Parallel::ForkManager->new($procs_per_node);
			# To do this, must write all output to files during the loop
			# after all codons finished, then we must read through and scoop results
			# while storing in the appropriate places
			
			
			### NEW POLYMORPHIC SITE ANALYSIS BOOTSTRAP ?
			#my %poly_site_data_hh; # $poly_site_data_hh{codon_index}->{unique_variants}
			
			my %codonum_group_numDefinedCodons_hh;
			
# BOOTSTRAP STORAGE	APPROACH	
##			my %gi_gj_codon_comp_stat_hh; # To store all comparisons for later bootstraps
			
			# STORE FOR SLIDING WINDOWS AND BOOTSTRAPPING
			my %codon_data_hh; # Will store codon->N_diffs_codon/S_diffs_codon/...
			
			PAIRWISE_DELETION: for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) { # for each codon in product
				$pm_codons->start and next; # this is IT
				
				my $codon_num = $codon_index+1;
				my $this_codon_output = '';
				
				#print "codon $codon_num at time " . time . "\n";
				
				#print CODON_FILE "$product_name\t$i\t$j\tcodon_" . $codon_num . " \t";
#				print CODON_FILE "$product_name\t$group_name_i\t$group_name_j\tcodon_" . $codon_num . " \t";
				$this_codon_output .= "ANALYSIS\t$product_name\t$group_name_i\t$group_name_j\tcodon_" . $codon_num . " \t";
				
#				my $curr_comp_num = 0;
				
				
				
				
				
				
				# Count number of defined codons in group i
				for (my $gi_spp_index = 0; $gi_spp_index < $num_seqs_per_codon_groupi; $gi_spp_index++) { # for each gj sequence at this codon
#					my $group_1_species = $group_1_spp[$gi_spp_index];
					my $this_species_codon = $groups_codons_aa[$i]->[$codon_index]->[$gi_spp_index];
#					my $this_species_codon = $data_codonnum_spp_codon_hh{$codon_num}->{$group_1_species};
					
					if(! ($this_species_codon =~ /N/ || $this_species_codon =~ /-/ || $this_species_codon eq '')) {
						$codonum_group_numDefinedCodons_hh{$codon_num}->{group1}++;
					}
				}
				
				# Count number of defined codons in group j
				for (my $gj_spp_index = 0; $gj_spp_index < $num_seqs_per_codon_groupj; $gj_spp_index++) { # for each gj sequence at this codon
#					my $group_2_species = $group_2_spp[$gj_spp_index];
					my $this_species_codon = $groups_codons_aa[$j]->[$codon_index]->[$gj_spp_index];
#					my $this_species_codon = $data_codonnum_spp_codon_hh{$codon_num}->{$group_2_species};
					
					if(! ($this_species_codon =~ /N/ || $this_species_codon =~ /-/ || $this_species_codon eq '')) {
						$codonum_group_numDefinedCodons_hh{$codon_num}->{group2}++;
					}
				}
				
				
				
				
				
				
				
				# ONLY IF POLYMORPHIC between-group (WILL THERE BE DIFFERENCES)
				if($polymorphic_codons_arr[$codon_index]) { # polymorphic codon
#					print CODON_FILE "polymorphic\t";
#					$this_codon_output .= "polymorphic\t";
					$this_codon_output .= "polymorphic\t" . $codonum_group_numDefinedCodons_hh{$codon_num}->{group1} . 
						"\t" . $codonum_group_numDefinedCodons_hh{$codon_num}->{group2} . "\t";
					
					my %comps_hh;
					
					# Sequence (codon) from group 1
					for (my $gi_seq_index = 0; $gi_seq_index < $num_seqs_per_codon_groupi; $gi_seq_index++) { # for each gi sequence at this codon
						
							
							# Sequence (codon) to compare against, group 2
							for (my $gj_seq_index = 0; $gj_seq_index < $num_seqs_per_codon_groupj; $gj_seq_index++) { # for each gj sequence at this codon
								#?my $group_1_species = $group_1_spp[$gi_spp_index];
								#?my $group_2_species = $group_2_spp[$gj_spp_index];
							
							
								my $codon_gi = $groups_codons_aa[$i]->[$codon_index]->[$gi_seq_index];
								my $codon_gj = $groups_codons_aa[$j]->[$codon_index]->[$gj_seq_index];
								
								if(! ($codon_gj =~ 'N' || $codon_gj =~ '-' || $codon_gj eq '' || $codon_gi =~ 'N' || $codon_gi =~ '-' || $codon_gi eq '')) { # they don't have to be the same; syn still stored here
#								if(! ($codon_gj =~ 'N') && ! ($codon_gj =~ '-')) { # they don't have to be the same; syn still stored here
									$comps_hh{$codon_gi}->{$codon_gj}+=1;
								}
							}
					
					}
					
					# ^ do this_actual_count here above?
					
					
					
					
					
					
					
					
					
					
					# SUM UP STUFF HERE
					my $between_g_comparisons = 0;
					my $sum_N_sites = 0;
					my $sum_S_sites = 0;
					my $sum_N_diffs = 0;
					my $sum_S_diffs = 0;
					
					#print"\nBetween groups $i and $j we compare";
					
					foreach my $codon_gi (keys %comps_hh) {
						foreach my $codon_gj (keys %{$comps_hh{$codon_gi}}) {
							
							if(($codon_gi =~ /-/ || $codon_gi =~ /N/) && ($codon_gj =~ /-/ || $codon_gj =~ /N/)) { # these will skip all additions
#							if($codon_gi =~ /-/ || $codon_gi =~ /N/ || $codon_gj =~ /-/ || $codon_gj =~ /N/) { # these will skip all additions
#								print CODON_FILE "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n";
#								$this_codon_output .= "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n"; # this was pointless because it's a private string
								next PAIRWISE_DELETION;
							} else {
								my $weight = $comps_hh{$codon_gi}->{$codon_gj};
								
								### NEW POLYMORPHIC SITE ANALYSIS ADDED
#								$poly_site_data_hh{$codon_index}->{$codon_si} += $weight;
#								$poly_site_data_hh{$codon_index}->{$codon_sj} += $weight;
								###
								
#								print CODON_FILE "($codon_gi\-$codon_gj)";
								$this_codon_output .= "($codon_gi\-$codon_gj)";
								
								# Sites codon gi
								my @codon_gi_sites_1_arr = &get_number_of_sites($codon_gi,1);
								my @codon_gi_sites_2_arr = &get_number_of_sites($codon_gi,2);
								my @codon_gi_sites_3_arr = &get_number_of_sites($codon_gi,3);
							
								my $codon_gi_N_sites = ($codon_gi_sites_1_arr[0] + $codon_gi_sites_2_arr[0] + $codon_gi_sites_3_arr[0]);
								my $codon_gi_S_sites = ($codon_gi_sites_1_arr[1] + $codon_gi_sites_2_arr[1] + $codon_gi_sites_3_arr[1]);							
								
								# Site codon gj
								my @codon_gj_sites_1_arr = &get_number_of_sites($codon_gj,1);
								my @codon_gj_sites_2_arr = &get_number_of_sites($codon_gj,2);
								my @codon_gj_sites_3_arr = &get_number_of_sites($codon_gj,3);
							
								my $codon_gj_N_sites = ($codon_gj_sites_1_arr[0] + $codon_gj_sites_2_arr[0] + $codon_gj_sites_3_arr[0]);
								my $codon_gj_S_sites = ($codon_gj_sites_1_arr[1] + $codon_gj_sites_2_arr[1] + $codon_gj_sites_3_arr[1]);
								
								#print "\nComparing codons $codon_gi and $codon_gj:\n";
								
								my $mean_comp_N_sites = ($codon_gi_N_sites + $codon_gj_N_sites) / 2;
								my $mean_comp_S_sites = ($codon_gi_S_sites + $codon_gj_S_sites) / 2;
								
								$sum_N_sites += $weight * $mean_comp_N_sites;
								$sum_S_sites += $weight * $mean_comp_S_sites;
								
								# Differences
								my $N_diffs = 0;
								my $S_diffs = 0;
								
								if($codon_gi ne $codon_gj) {
									my @diffs_arr = &return_avg_diffs($codon_gi,$codon_gj);
									$N_diffs = $diffs_arr[0];
									$S_diffs = $diffs_arr[1];
								}
								
								$sum_N_diffs += $weight * $N_diffs;
								$sum_S_diffs += $weight * $S_diffs;							
	
								$between_g_comparisons += $weight;
								
#								$this_codon_data .= "$weight\t$mean_comp_N_sites\t$mean_comp_S_sites\t" .
#									"$N_diffs\t$S_diffs";
									
									
									
								
# BOOTSTRAP STORAGE	APPROACH								
##								#my @gi_gj_codon_mm_aa; #gi_gj_codon_mm_aa[group1_index]->[group2_index]->[codon_index]->[N_diffs/S_diffs/N_sites/S_sites]->array of all values for that stat
##								if($num_bootstraps > 1) {
##									print "\ngroup $i / group $j / codon $codon_index $codon_gi\-$codon_gj / weight $weight starting with comp $curr_comp_num and ending with comp ";
##									
##									# STORE HERE FOR BOOTSTRAPS
##									for(my $k = 1; $k <= $weight; $k++) {
##										
##										
##										
##										$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{N_diffs} = $N_diffs; # curr_comp_num starts at 0
##										$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{S_diffs} = $S_diffs;
##										$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{N_sites} = $mean_comp_N_sites;
##										$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{S_sites} = $mean_comp_S_sites;
##										
##										$curr_comp_num ++;
##										
##										#push(@{$gi_gj_codon_mm_aa[$i]->[$j]->[$codon_index]->[N_diffs]},$N_diffs);
##										#push(@{$gi_gj_codon_mm_aa[$i]->[$j]->[$codon_index]->[S_diffs]},$S_diffs);
##										#push(@{$gi_gj_codon_mm_aa[$i]->[$j]->[$codon_index]->[N_sites]},$mean_comp_N_sites);
##										#push(@{$gi_gj_codon_mm_aa[$i]->[$j]->[$codon_index]->[S_sites]},$mean_comp_S_sites);
##									}
##									
##									print "$curr_comp_num\n";
##								}
								#my $random_codon_index = $gi_gj_codon_mm_aa[$i]->[$j]->[int(rand(@{$gi_gj_codon_mm_aa[$i]->[$j]}))];
								#my $random_N_diffs = $gi_gj_codon_mm_aa[$i]->[$j]->[$random_codon_index];
								#my $random_S_diffs;
								#my $random_N_sites;
								#my $random_S_sites;
								
								#my $RAND_COMP_NUM = int(rand($THIS_CODON_NUM_COMPS));
								#my $random_N_diffs = $gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$RAND_COMP_NUM}->{N_diffs}
								
									##groups_codons_aa[group_index]->[codon_index]->array of all codons in group
									#my @groups_codons_aa;
									#push(@{$groups_codons_aa[$group_index]->[$array_index]},$codon);
									#my $codon_gi = $groups_codons_aa[$i]->[$codon_index]->[$gi_seq_index];
									#my $random_codon = $codonum_codon_aa[$codon_index]->[int(rand(@{$codonum_codon_aa[$codon_index]}))];
								
														
							}
						}
					}
#					print CODON_FILE "\t";
					$this_codon_output .= "\t";
					
					#print "\nfor a total of $between_g_comparisons comparisons\n";
					my $mean_N_sites = 0;
					my $mean_S_sites = 0;
					my $mean_N_diffs = 0;
					my $mean_S_diffs = 0;	
					
					if($between_g_comparisons > 0) {
						$mean_N_sites = ($sum_N_sites / $between_g_comparisons);
						$mean_S_sites = ($sum_S_sites / $between_g_comparisons);
						$mean_N_diffs = ($sum_N_diffs / $between_g_comparisons);
						$mean_S_diffs = ($sum_S_diffs / $between_g_comparisons);					
					}
					
#					print CODON_FILE "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";
					$this_codon_output .= "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs";
					
					# SLIDING WINDOW & BOOTSTRAP
					$codon_data_hh{$codon_num}->{N_sites} = $mean_N_sites;
					$codon_data_hh{$codon_num}->{S_sites} = $mean_S_sites;
					$codon_data_hh{$codon_num}->{N_diffs} = $mean_N_diffs;
					$codon_data_hh{$codon_num}->{S_diffs} = $mean_S_diffs;
					
					
					###
#					$between_product_N_sites_sum += $mean_N_sites;
#					$between_product_S_sites_sum += $mean_S_sites;
#					$between_product_N_diffs_sum += $mean_N_diffs;
#					$between_product_S_diffs_sum += $mean_S_diffs;
					###
					
					
					
					#print ".. and it's all done.\n";
					
				} else { # not polymorphic; just use first DEFINED codon from the first group because it's conserved
					#my $conserved_codon = $groups_codons_aa[$i]->[$codon_index]->[0]; # grab first codon10
					
					my $conserved_codon = '';
##					FIND_CONSERVED_CODON: for(my $k = 0; $k < $num_seqs_per_codon_groupi; $k++) { # for each codon in product
##						
##						my $curr_codon_rep = $groups_codons_aa[$i]->[$codon_index]->[$k]; # grab first codon THAT DOESN'T CONTAIN N or GAP
##						if($curr_codon_rep =~ 'N' || $curr_codon_rep =~ '-') {
##							next FIND_CONSERVED_CODON;
##						} else {
##							$conserved_codon = $curr_codon_rep;
##							last FIND_CONSERVED_CODON;
##						}
##					}
					
					
					
					
					
					FIND_CONSERVED_CODON: for(my $gi_spp_index = 0; $gi_spp_index < $num_seqs_per_codon_groupi; $gi_spp_index++) {
						
						FIND_CONSERVED_G2: for (my $gj_spp_index = 0; $gj_spp_index < $num_seqs_per_codon_groupj; $gj_spp_index++) { # for each gj sequence at this codon
							
##							my $group_1_species = $group_1_spp[$gi_spp_index];
##							my $group_2_species = $group_2_spp[$gj_spp_index];
##							
##							my $codon_gi_rep = $data_codonnum_spp_codon_hh{$codon_num}->{$group_1_species};
##							my $codon_gj_rep = $data_codonnum_spp_codon_hh{$codon_num}->{$group_2_species};
							
							my $codon_gi_rep = $groups_codons_aa[$i]->[$codon_index]->[$gi_spp_index];
							my $codon_gj_rep = $groups_codons_aa[$j]->[$codon_index]->[$gj_spp_index];
							
#											if($codon_num == 2 && $partition_id == 1) {
#												print "\ncodon 2\ngroup_1_species=$group_1_species codon=$codon_gi_rep\ngroup_2_species=$group_2_species codon=$codon_gj_rep\n";
#											}
							
							# grab first codon THAT DOESN'T CONTAIN N or GAP
							if($codon_gi_rep =~ 'N' || $codon_gi_rep =~ '-' || $codon_gi_rep eq '') {
								if($codon_gj_rep =~ 'N' || $codon_gj_rep =~ '-' || $codon_gj_rep eq '') {
									#print "\ncodon $codon_num codons are codon_gi_rep $codon_gi_rep and codon_gj_rep $codon_gj_rep and neither work\n";
									next FIND_CONSERVED_G2;
								} else {
									$conserved_codon = $codon_gj_rep;
#													if($codon_num == 2 && $partition_id == 1) {
#														print "\ncodon $codon_num codons are codon_gi_rep $codon_gi_rep and codon_gj_rep $codon_gj_rep and the last works\n";
#													}
									last FIND_CONSERVED_CODON;
								}
							} else {
								$conserved_codon = $codon_gi_rep;
#												if($codon_num == 2 && $partition_id == 1) {
#													print "\ncodon $codon_num codons are codon_gi_rep $codon_gi_rep and codon_gj_rep $codon_gj_rep and the first works\n";
#												}
								last FIND_CONSERVED_CODON;
							}
						}
					}
					
					# Print a warning if no conserved codon was found
					if($conserved_codon eq '') {
						warn "\n### WARNING: no codons found in between-group analysis at codon position $codon_num\.\n".
							"### Most likely, the selected groups consist of all gaps at this codon, which is not necessarily a problem.\n".
							"### Inserting gaps and proceeding.\n\n";
						
						$conserved_codon = '---';
					}
					
					
					
					#STORE THIS CONSERVED CODON FOR LATER USE IN BOOTSTRAPPING
#					$codonIndex_conserved{$codon_index} = $conserved_codon;
					#HERE
					
					my @codon_sites_1_arr = &get_number_of_sites($conserved_codon,1);
					my @codon_sites_2_arr = &get_number_of_sites($conserved_codon,2);
					my @codon_sites_3_arr = &get_number_of_sites($conserved_codon,3);
					
					my $codon_N_sites = ($codon_sites_1_arr[0] + $codon_sites_2_arr[0] + $codon_sites_3_arr[0]);
					my $codon_S_sites = ($codon_sites_1_arr[1] + $codon_sites_2_arr[1] + $codon_sites_3_arr[1]);
					
#					print CODON_FILE "conserved\t$conserved_codon\t".
#						"$codon_N_sites\t$codon_S_sites\t0\t0\n";
##					$this_codon_output .= "conserved\t$conserved_codon\t".
##						"$codon_N_sites\t$codon_S_sites\t0\t0";
					
#					$between_product_N_sites_sum += $codon_N_sites;
#					$between_product_S_sites_sum += $codon_S_sites;
					
					
					$this_codon_output .= "conserved\t" . $codonum_group_numDefinedCodons_hh{$codon_num}->{group1} . 
						"\t" . $codonum_group_numDefinedCodons_hh{$codon_num}->{group2} . "\t$conserved_codon\t".
						"$codon_N_sites\t$codon_S_sites\t0\t0";
					
					# SLIDING WINDOW & BOOTSTRAP
					$codon_data_hh{$codon_num}->{N_sites} = $codon_N_sites;
					$codon_data_hh{$codon_num}->{S_sites} = $codon_S_sites;
					
					
					#if($product_name =~ /2K/) {
					#	print "\nProduct $product_name, conserved codon $conserved_codon N sites $codon_N_sites\n";
					#}
					
# BOOTSTRAP STORAGE APPROACH
##					# STORE HERE FOR BOOTSTRAPS: just one value	
##					if($num_bootstraps > 1) {
##						$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{N_diffs} = 0; # curr_comp_num starts at 0
##						$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{S_diffs} = 0;
##						$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{N_sites} = $codon_N_sites;
##						$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{S_sites} = $codon_S_sites;
##					}
									
					
					
				} # end not polymorphic
				
				open(THIS_CODON_TEMP_FILE,">>$codon_index");
				print THIS_CODON_TEMP_FILE "$this_codon_output";
				close THIS_CODON_TEMP_FILE;
				
				$pm_codons->finish; # special name
				
			} # end all codons in product

			# PARALLEL
			$pm_codons->wait_all_children; # special name, methods within module
	
			
			# Calculate overall values for the product
			my $between_product_N_sites_sum = 0;
			my $between_product_S_sites_sum = 0;
			my $between_product_N_diffs_sum = 0;
			my $between_product_S_diffs_sum = 0;
							
			
			# Vacuum up output and delete temp files	
			my @codon_lines_arr;	
			my @codon_analysis_FILES_arr = glob "*";
			@codon_analysis_FILES_arr = sort {$a <=> $b} @codon_analysis_FILES_arr;
			
			#print "sorted polymorphic_codons_FILE_arr: @polymorphic_codons_FILES_arr\n";
			
			foreach(@codon_analysis_FILES_arr) { # file names, actually
				my $codon_file = $_; # this is the INDEX
				
				open(CURR_CODON_FILE, $codon_file) or die "\n## Cannot open $codon_file. TERMINATED.\n\n";
				while(<CURR_CODON_FILE>) {
					chomp;
					
					push(@codon_lines_arr,$_);
					
					my @this_line_arr = split(/\t/,$_,-1);
					#analysis / product / group_1 / group_2 / codon / variability / comparisons / N_sites / S_sites / N_diffs / S_diffs
					my $poly = $this_line_arr[5]; #?
					
					# IF POLY
					if($poly eq 'polymorphic') {
						$between_product_N_sites_sum += $this_line_arr[9];
						$between_product_S_sites_sum += $this_line_arr[10];
						$between_product_N_diffs_sum += $this_line_arr[11];
						$between_product_S_diffs_sum += $this_line_arr[12];
					
					# IF CONSERVED
					} elsif($poly eq 'conserved') {
						$between_product_N_sites_sum += $this_line_arr[9];
						$between_product_S_sites_sum += $this_line_arr[10];
						
						#STORE THIS CONSERVED CODON FOR LATER USE IN BOOTSTRAPPING
##						$codonIndex_conserved{$codon_file} = $this_line_arr[5]; # COMPARISONS column
					# PROBLEM
					} else {
						die "\nPROBLEM: site neither conserved nor polymorphic!? DIE\n\n";
					}
					
					
					
		
		
				}
				close CURR_CODON_FILE;
				unlink $codon_file;
			}
			
			chdir("..");
			rmdir("$product_name\_codon\_analysis");
			# done sweeping output

			# Concatenate and print parallelized codon results to a results file
			open(CODON_FILE,">>between\_group\_codon\_results\.txt");
			foreach(@codon_lines_arr) {
				print CODON_FILE "$_\n";
			}
			close CODON_FILE;
			
			
			

			
			
			
# BOOTSTRAP PSEUDOSAMPLE APPROACH
##########################################################################################			
			#################
			# BOOTSTRAPPING: SIMPLY CREATE these groups' pseudosamples; calculate values;
			# store values and discard pseudosample; will compare later
			
			# WE WILL BOOTSTRAP FOR THIS (group i vs. group j) comparison
			# Each bootstrap will consist of two pseudosamples: one from group i, one from group j
			
			
#			my $var_piN_minus_piS;
#			my $SE_piN_minus_piS;
			
			
			# I DON'T THINK THIS IS RIGHT. USE WITHOUT BOOTSTRAPS, THEN RUN
			# snpgenie.sw.pl
			
##BOOTSTRAP##				if($num_bootstraps > 1) {
##BOOTSTRAP##				
##BOOTSTRAP###				# MAKE BOOTSTRAP FILE
##BOOTSTRAP###				open(BOOT_FILE_PRODUCT,">>$product_name\_bootstrap\_results\.txt");
##BOOTSTRAP###				print BOOT_FILE_PRODUCT "product_name\tbootstrap_num\tsim_product_N_sites_sum\t".
##BOOTSTRAP###						"sim_product_S_sites_sum\tsim_product_N_diffs_sum\tsim_product_S_diffs_sum\n";
##BOOTSTRAP##				
##BOOTSTRAP##				#################
##BOOTSTRAP##				# BOOTSTRAPPING TO CALCULATE STANDARD ERROR HERE?
##BOOTSTRAP##				print "\nBootstrapping for standard errors, n = $num_bootstraps\...\n\n";
##BOOTSTRAP##				
##BOOTSTRAP##				#my @simulated_sequences_aa; # first index is codon position; second index is random sequence
##BOOTSTRAP##				# $simulated_sequences_aa[CODON_INDEX]->[RANDSEQ_INDEX]
##BOOTSTRAP##				# @{$simulated_sequences_aa[CODON_INDEX]}
##BOOTSTRAP##				#print "num_codons_in_product: $num_codons_in_product\n";
##BOOTSTRAP##				
##BOOTSTRAP##				my @sim_N_sites_arr;
##BOOTSTRAP##				my @sim_S_sites_arr;
##BOOTSTRAP##				my @sim_N_diffs_arr;
##BOOTSTRAP##				my @sim_S_diffs_arr;
##BOOTSTRAP##				
##BOOTSTRAP##				
##BOOTSTRAP##				# PARALLELIZE IT
##BOOTSTRAP##				mkdir("$product_name\_bootstrap_temp_files");
##BOOTSTRAP##				chdir("$product_name\_bootstrap_temp_files");
##BOOTSTRAP##				#my $procs = 60; # number to do at once, cores to assign simultaneously, machine-limited. CUVIER has 80
##BOOTSTRAP##				my $pm = Parallel::ForkManager->new($procs_per_node);
##BOOTSTRAP##				
##BOOTSTRAP##				for(my $bootstrap_num = 1; $bootstrap_num <= $num_bootstraps; $bootstrap_num++) {
##BOOTSTRAP##					$pm->start and next; # this is IT
##BOOTSTRAP##					
##BOOTSTRAP##					srand();
##BOOTSTRAP##					
##BOOTSTRAP##					# New bootstrap run, of $num_bootstraps
##BOOTSTRAP##					#print "bootstrap $bootstrap_num\n";
##BOOTSTRAP##					
##BOOTSTRAP##					#my @simulated_sequences_aa; # first index is codon position; second index is random sequence
##BOOTSTRAP##					# $simulated_sequences_aa[CODON_INDEX]->[RANDSEQ_INDEX]
##BOOTSTRAP##					# @{$simulated_sequences_aa[CODON_INDEX]}
##BOOTSTRAP##					#print "num_codons_in_product: $num_codons_in_product\n";
##BOOTSTRAP##					
##BOOTSTRAP##					my $sim_product_N_sites_sum = 0;
##BOOTSTRAP##					my $sim_product_S_sites_sum = 0;
##BOOTSTRAP##					my $sim_product_N_diffs_sum = 0;
##BOOTSTRAP##					my $sim_product_S_diffs_sum = 0;
##BOOTSTRAP##					
##BOOTSTRAP##					#my $random_sum = 0;
##BOOTSTRAP##					
##BOOTSTRAP##					SIM_PAIRWISE_DELETION: for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) { # for each codon in product
##BOOTSTRAP##						my @simulated_codons_groupi;
##BOOTSTRAP##						my @simulated_codons_groupj;
##BOOTSTRAP##						#my $codon_num = $codon_index+1;
##BOOTSTRAP##						#print "\tcodon_index: $codon_index\n";
##BOOTSTRAP##						
##BOOTSTRAP##						if($polymorphic_codons_arr[$codon_index]) { # if codon is polymorphic in ACTUAL DATA, value is 1
##BOOTSTRAP##						
##BOOTSTRAP##							##### HERE, WHERE ARE ALREADY WITHIN A PARTICULAR GROUP PAIR (i & j) COMPARISON
##BOOTSTRAP##							# index for group i is i, j is j. Easy	
##BOOTSTRAP##						
##BOOTSTRAP##							# GENERATE PSEUDOSAMPLE FOR GROUP I
##BOOTSTRAP##							for(my $seq_index=0; $seq_index<$num_seqs_per_codon_groupi; $seq_index++) {
##BOOTSTRAP##							#for(my $i=0; $i<$num_bootstraps; $i++) { ### <-- this actually should be num SEQUENCES; this whole thing is JUST ONE bootstrap!
##BOOTSTRAP##								#print "i: $i\n";
##BOOTSTRAP##								# Choose random codon from all sequences. I think it should contain N, since that's how the real data are
##BOOTSTRAP##								# ACTUAL codons are stored in @{$codonum_codon_aa{codon_index}}->array of all codons at that site
##BOOTSTRAP##								
##BOOTSTRAP##								#srand();
##BOOTSTRAP##								my $rand_seq_index = int(rand($num_seqs_per_codon_groupi));
##BOOTSTRAP##								#my $random_codon = $groups_codons_aa[$i]->[$codon_index]->[int(rand(@{$groups_codons_aa[$i]->[$codon_index]}))];
##BOOTSTRAP##								my $random_codon = $groups_codons_aa[$i]->[$codon_index]->[$rand_seq_index];
##BOOTSTRAP##				
##BOOTSTRAP##								#if($product_name eq '2K peptide') {
##BOOTSTRAP##								#	print "\nFor bootstrapping group i=$i, we randomly chose index $rand_seq_index with codon $random_codon\n";
##BOOTSTRAP##								#}
##BOOTSTRAP##								
##BOOTSTRAP##								#print "\t\tFor bootstrapping group i=$i with N=$num_seqs_per_codon_groupi\,\n\t\t we randomly chose index $rand_seq_index with codon $random_codon\n";
##BOOTSTRAP##								
##BOOTSTRAP##								push(@simulated_codons_groupi,$random_codon);
##BOOTSTRAP##							}
##BOOTSTRAP##							
##BOOTSTRAP##							# GENERATE PSEUDOSAMPLE FOR GROUP J
##BOOTSTRAP##							for(my $seq_index=0; $seq_index<$num_seqs_per_codon_groupj; $seq_index++) {
##BOOTSTRAP##								#print "i: $i\n";
##BOOTSTRAP##								# Choose random codon from all sequences. I think it should contain N, since that's how the real data are
##BOOTSTRAP##								
##BOOTSTRAP##								#srand();
##BOOTSTRAP##								my $rand_seq_index = int(rand($num_seqs_per_codon_groupj));
##BOOTSTRAP##								my $random_codon = $groups_codons_aa[$j]->[$codon_index]->[$rand_seq_index];			
##BOOTSTRAP##								
##BOOTSTRAP##								#print "\t\tFor bootstrapping group j=$j with N=$num_seqs_per_codon_groupj\,\n\t\t we randomly chose index $rand_seq_index with codon $random_codon\n";
##BOOTSTRAP##								
##BOOTSTRAP##								push(@simulated_codons_groupj,$random_codon);
##BOOTSTRAP##							}
##BOOTSTRAP##						
##BOOTSTRAP##						}
##BOOTSTRAP##						
##BOOTSTRAP##						
##BOOTSTRAP##						
##BOOTSTRAP##						
##BOOTSTRAP##						
##BOOTSTRAP##						
##BOOTSTRAP##						
##BOOTSTRAP##						# Now perform the actual calculation
##BOOTSTRAP##						# We need mean dN and dS between all pairs of simulated codons
##BOOTSTRAP##						
##BOOTSTRAP##						if($polymorphic_codons_arr[$codon_index]) { # if codon is polymorphic in ACTUAL DATA, value is 1
##BOOTSTRAP##						# There is a POSSIBILITY it's polymorphic between the pseudosamples
##BOOTSTRAP##						#	print CODON_FILE "polymorphic\t";
##BOOTSTRAP##				
##BOOTSTRAP##							my %sim_comps_hh;
##BOOTSTRAP##							for (my $si_seq_index = 0; $si_seq_index < $num_seqs_per_codon_groupi; $si_seq_index++) { # for each si sequence at this codon
##BOOTSTRAP##								my $codon_si = $simulated_codons_groupi[$si_seq_index];
##BOOTSTRAP##								
##BOOTSTRAP##								#print "\t\tThere are num_seqs_per_codon_groupi=$num_seqs_per_codon_groupi and scalar(simulated_codons_groupi)=".scalar(@simulated_codons_groupi)."\n";
##BOOTSTRAP##								
##BOOTSTRAP##								# Sequence (codon) to compare against
##BOOTSTRAP##								if(! ($codon_si =~ 'N') && ! ($codon_si =~ '-')) {
##BOOTSTRAP##									for (my $sj_seq_index = 0; $sj_seq_index < $num_seqs_per_codon_groupj; $sj_seq_index++) { # for each sj sequence at this codon
##BOOTSTRAP##										my $codon_sj = $simulated_codons_groupj[$sj_seq_index];
##BOOTSTRAP##										
##BOOTSTRAP##										#print "\t\tThere are num_seqs_per_codon_groupj=$num_seqs_per_codon_groupj and scalar(simulated_codons_groupj)=".scalar(@simulated_codons_groupj)."\n";
##BOOTSTRAP##										
##BOOTSTRAP##										if(! ($codon_sj =~ 'N') && ! ($codon_sj =~ '-')) { # they don't have to be the same; syn still stored here
##BOOTSTRAP##											$sim_comps_hh{$codon_si}->{$codon_sj}+=1;
##BOOTSTRAP##										}
##BOOTSTRAP##									}
##BOOTSTRAP##								}
##BOOTSTRAP##							}
##BOOTSTRAP##													
##BOOTSTRAP##							# SUM UP STUFF HERE
##BOOTSTRAP##							my $between_s_comparisons = 0;
##BOOTSTRAP##							my $sim_sum_N_sites = 0;
##BOOTSTRAP##							my $sim_sum_S_sites = 0;
##BOOTSTRAP##							my $sim_sum_N_diffs = 0;
##BOOTSTRAP##							my $sim_sum_S_diffs = 0;
##BOOTSTRAP##							
##BOOTSTRAP##							foreach my $codon_si (keys %sim_comps_hh) {
##BOOTSTRAP##								foreach my $codon_sj (keys %{$sim_comps_hh{$codon_si}}) {
##BOOTSTRAP##									
##BOOTSTRAP##									if($codon_si =~ /-/ || $codon_si =~ /N/ || $codon_sj =~ /-/ || $codon_sj =~ /N/) { # these will skip all additions
##BOOTSTRAP##				#						print CODON_FILE "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n";
##BOOTSTRAP##										next SIM_PAIRWISE_DELETION;
##BOOTSTRAP##									} else {						
##BOOTSTRAP##										my $weight = $sim_comps_hh{$codon_si}->{$codon_sj};
##BOOTSTRAP##										
##BOOTSTRAP##										### NEW POLYMORPHIC SITE ANALYSIS ADDED
##BOOTSTRAP###										$poly_site_data_hh{$codon_index}->{$codon_si} += $weight;
##BOOTSTRAP###										$poly_site_data_hh{$codon_index}->{$codon_sj} += $weight;
##BOOTSTRAP##										###
##BOOTSTRAP##										
##BOOTSTRAP##					#					print CODON_FILE "($codon_si\-$codon_sj)";
##BOOTSTRAP##										
##BOOTSTRAP##										# Sites codon gi
##BOOTSTRAP##										my @codon_si_sites_1_arr = &get_number_of_sites($codon_si,1);
##BOOTSTRAP##										my @codon_si_sites_2_arr = &get_number_of_sites($codon_si,2);
##BOOTSTRAP##										my @codon_si_sites_3_arr = &get_number_of_sites($codon_si,3);
##BOOTSTRAP##									
##BOOTSTRAP##										my $codon_si_N_sites = ($codon_si_sites_1_arr[0] + $codon_si_sites_2_arr[0] + $codon_si_sites_3_arr[0]);
##BOOTSTRAP##										my $codon_si_S_sites = ($codon_si_sites_1_arr[1] + $codon_si_sites_2_arr[1] + $codon_si_sites_3_arr[1]);							
##BOOTSTRAP##										
##BOOTSTRAP##										# Site codon gj
##BOOTSTRAP##										my @codon_sj_sites_1_arr = &get_number_of_sites($codon_sj,1);
##BOOTSTRAP##										my @codon_sj_sites_2_arr = &get_number_of_sites($codon_sj,2);
##BOOTSTRAP##										my @codon_sj_sites_3_arr = &get_number_of_sites($codon_sj,3);
##BOOTSTRAP##									
##BOOTSTRAP##										my $codon_sj_N_sites = ($codon_sj_sites_1_arr[0] + $codon_sj_sites_2_arr[0] + $codon_sj_sites_3_arr[0]);
##BOOTSTRAP##										my $codon_sj_S_sites = ($codon_sj_sites_1_arr[1] + $codon_sj_sites_2_arr[1] + $codon_sj_sites_3_arr[1]);
##BOOTSTRAP##										
##BOOTSTRAP##										#print "\nComparing codons $codon_si and $codon_sj:\n";
##BOOTSTRAP##										
##BOOTSTRAP##										my $mean_comp_N_sites = ($codon_si_N_sites + $codon_sj_N_sites) / 2;
##BOOTSTRAP##										my $mean_comp_S_sites = ($codon_si_S_sites + $codon_sj_S_sites) / 2;
##BOOTSTRAP##										
##BOOTSTRAP##										#print "\t\tN_sites here=$mean_comp_N_sites weight here=$weight\n";
##BOOTSTRAP##										
##BOOTSTRAP##										$sim_sum_N_sites += $weight * $mean_comp_N_sites;
##BOOTSTRAP##										$sim_sum_S_sites += $weight * $mean_comp_S_sites;
##BOOTSTRAP##										
##BOOTSTRAP##										# Differences
##BOOTSTRAP##										my $N_diffs = 0;
##BOOTSTRAP##										my $S_diffs = 0;
##BOOTSTRAP##										
##BOOTSTRAP##										if($codon_si ne $codon_sj) {
##BOOTSTRAP##											my @diffs_arr = &return_avg_diffs($codon_si,$codon_sj);
##BOOTSTRAP##											$N_diffs = $diffs_arr[0];
##BOOTSTRAP##											$S_diffs = $diffs_arr[1];
##BOOTSTRAP##										}
##BOOTSTRAP##										
##BOOTSTRAP##										$sim_sum_N_diffs += $weight * $N_diffs;
##BOOTSTRAP##										$sim_sum_S_diffs += $weight * $S_diffs;							
##BOOTSTRAP##				#	
##BOOTSTRAP##										$between_s_comparisons += $weight;
##BOOTSTRAP##									}
##BOOTSTRAP##								}
##BOOTSTRAP##							}
##BOOTSTRAP##				#			print CODON_FILE "\t";
##BOOTSTRAP##							
##BOOTSTRAP##							#print "\nfor a total of $between_s_comparisons comparisons\n";
##BOOTSTRAP##							
##BOOTSTRAP##							my $mean_N_sites = ($sim_sum_N_sites / $between_s_comparisons);
##BOOTSTRAP##							my $mean_S_sites = ($sim_sum_S_sites / $between_s_comparisons);
##BOOTSTRAP##							my $mean_N_diffs = ($sim_sum_N_diffs / $between_s_comparisons);
##BOOTSTRAP##							my $mean_S_diffs = ($sim_sum_S_diffs / $between_s_comparisons);
##BOOTSTRAP##							
##BOOTSTRAP##							#print "\t\tMean_N_sites here=$mean_N_sites between_s_comparisons here=$between_s_comparisons\n";
##BOOTSTRAP##							
##BOOTSTRAP##				#			print CODON_FILE "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";
##BOOTSTRAP##							
##BOOTSTRAP##							$sim_product_N_sites_sum += $mean_N_sites;
##BOOTSTRAP##							$sim_product_S_sites_sum += $mean_S_sites;
##BOOTSTRAP##							$sim_product_N_diffs_sum += $mean_N_diffs;
##BOOTSTRAP##							$sim_product_S_diffs_sum += $mean_S_diffs;
##BOOTSTRAP##							
##BOOTSTRAP##							#print ".. and it's all done.\n";
##BOOTSTRAP##							
##BOOTSTRAP##						} else { # not polymorphic; just use first codon from the first group because it's conserved
##BOOTSTRAP##							
##BOOTSTRAP##							#my $conserved_codon = '';
##BOOTSTRAP##							
##BOOTSTRAP##							my $conserved_codon = $codonIndex_conserved{$codon_index};
##BOOTSTRAP##							
##BOOTSTRAP##							#print "\t\tThe conserved codon is: $conserved_codon\n";
##BOOTSTRAP##							
##BOOTSTRAP###							FIND_CONSERVED_CODON_SIM: for(my $k = 0; $k < $num_seqs_per_codon_groupi; $k++) { # for each codon in product
##BOOTSTRAP###								my $curr_codon_rep = $simulated_codons_groupi[$k]; # grab first codon THAT DOESN'T CONTAIN N/-
##BOOTSTRAP###								if($curr_codon_rep =~ "N" || $curr_codon_rep =~ "-") {
##BOOTSTRAP###									next FIND_CONSERVED_CODON_SIM;
##BOOTSTRAP###								} else {
##BOOTSTRAP###									$conserved_codon = $curr_codon_rep;
##BOOTSTRAP###									last FIND_CONSERVED_CODON_SIM;
##BOOTSTRAP###								}
##BOOTSTRAP###							}
##BOOTSTRAP##							
##BOOTSTRAP##							my @codon_sites_1_arr = &get_number_of_sites($conserved_codon,1);
##BOOTSTRAP##							my @codon_sites_2_arr = &get_number_of_sites($conserved_codon,2);
##BOOTSTRAP##							my @codon_sites_3_arr = &get_number_of_sites($conserved_codon,3);
##BOOTSTRAP##							
##BOOTSTRAP##							my $codon_N_sites = ($codon_sites_1_arr[0] + $codon_sites_2_arr[0] + $codon_sites_3_arr[0]);
##BOOTSTRAP##							my $codon_S_sites = ($codon_sites_1_arr[1] + $codon_sites_2_arr[1] + $codon_sites_3_arr[1]);
##BOOTSTRAP##							
##BOOTSTRAP##			#				print CODON_FILE "conserved\t$conserved_codon\t".
##BOOTSTRAP##			#					"$codon_N_sites\t$codon_S_sites\t0\t0\n";
##BOOTSTRAP##							
##BOOTSTRAP##							$sim_product_N_sites_sum += $codon_N_sites;
##BOOTSTRAP##							$sim_product_S_sites_sum += $codon_S_sites;
##BOOTSTRAP##							
##BOOTSTRAP##							#if($product_name =~ /2K/) {
##BOOTSTRAP##							#	print "\nProduct $product_name, conserved codon $conserved_codon N sites $codon_N_sites\n";
##BOOTSTRAP##							#}
##BOOTSTRAP##							
##BOOTSTRAP##						}
##BOOTSTRAP##					}
##BOOTSTRAP##					
##BOOTSTRAP##					# Print product totals for BOOTSTRAP
##BOOTSTRAP##					my $out_line_boot = "$product_name\t$bootstrap_num\t$sim_product_N_sites_sum\t".
##BOOTSTRAP##						"$sim_product_S_sites_sum\t$sim_product_N_diffs_sum\t$sim_product_S_diffs_sum\n";
##BOOTSTRAP##					
##BOOTSTRAP###					print BOOT_FILE_PRODUCT "$out_line_boot";
##BOOTSTRAP##					
##BOOTSTRAP##					#print "$out_line_boot";
##BOOTSTRAP##					
##BOOTSTRAP##					push(@sim_N_sites_arr,$sim_product_N_sites_sum);
##BOOTSTRAP##					push(@sim_S_sites_arr,$sim_product_S_sites_sum);
##BOOTSTRAP##					push(@sim_N_diffs_arr,$sim_product_N_diffs_sum);
##BOOTSTRAP##					push(@sim_S_diffs_arr,$sim_product_S_diffs_sum);
##BOOTSTRAP##					
##BOOTSTRAP##					#print $random_sum . "\n";
##BOOTSTRAP##					
##BOOTSTRAP##					open(THIS_BOOTSTRAP_TEMP_FILE,">>temp\_between\_group\_bootstrap\_$bootstrap_num\.txt");
##BOOTSTRAP##					#PRINT
##BOOTSTRAP##					#print THIS_BOOTSTRAP_TEMP_FILE "$bootstrap_num\t$i\t$j\t".
##BOOTSTRAP##					print THIS_BOOTSTRAP_TEMP_FILE "$bootstrap_num\t$group_name_i\t$group_name_j\t".
##BOOTSTRAP##						"$sim_product_N_diffs_sum\t$sim_product_S_diffs_sum\t".
##BOOTSTRAP##						"$sim_product_N_sites_sum\t$sim_product_S_sites_sum\n";
##BOOTSTRAP##					close THIS_BOOTSTRAP_TEMP_FILE;
##BOOTSTRAP##					
##BOOTSTRAP##					$pm->finish; # special name
##BOOTSTRAP##
##BOOTSTRAP##				} # end last bootstrap
##BOOTSTRAP##				
##BOOTSTRAP##				$pm->wait_all_children; # special name, methods within module
##BOOTSTRAP##				chdir("..");
##BOOTSTRAP##				
##BOOTSTRAP##				#close BOOT_FILE_PRODUCT;
##BOOTSTRAP##				
##BOOTSTRAP##				my $actual_piN = $between_product_N_diffs_sum / $between_product_N_sites_sum;
##BOOTSTRAP##				my $actual_piS = $between_product_S_diffs_sum / $between_product_S_sites_sum;
##BOOTSTRAP##				my $actual_piN_minus_piS = $actual_piN - $actual_piS;
##BOOTSTRAP##				
##BOOTSTRAP##				my $sum_sq_diffs;
##BOOTSTRAP##				# CALCULATE BOOTSTRAP VALUES HERE! NEI & KUMAR EQUATIONS
##BOOTSTRAP##				for(my $sim_num = 0; $sim_num < scalar(@sim_N_sites_arr); $sim_num++) {
##BOOTSTRAP##					
##BOOTSTRAP##					#if($product_name eq '2K peptide') {
##BOOTSTRAP##					#	print "\nThis $product_name number of N sites is: " . $sim_N_sites_arr[$sim_num] . "\n";
##BOOTSTRAP##					#}
##BOOTSTRAP##					
##BOOTSTRAP##					my $this_round_piN = $sim_N_diffs_arr[$sim_num] / $sim_N_sites_arr[$sim_num];
##BOOTSTRAP##					my $this_round_piS = $sim_S_diffs_arr[$sim_num] / $sim_S_sites_arr[$sim_num];
##BOOTSTRAP##					my $this_round_piN_minus_piS = $this_round_piN - $this_round_piS;
##BOOTSTRAP##					
##BOOTSTRAP##					$sum_sq_diffs += ($this_round_piN_minus_piS - $actual_piN_minus_piS)**2;
##BOOTSTRAP##				}
##BOOTSTRAP##			
##BOOTSTRAP###				$var_piN_minus_piS = $sum_sq_diffs / (scalar(@sim_N_sites_arr));
##BOOTSTRAP###				$SE_piN_minus_piS = $var_piN_minus_piS ** (1/2);
##BOOTSTRAP##				
##BOOTSTRAP##				#print "\nSE for this product is: $SE_piN_minus_piS\n\n";
##BOOTSTRAP##				
##BOOTSTRAP##			} # end case in which BOOTSTRAPS are calculated using pseudosamples
##########################################################################################








# BOOTSTRAP STORAGE APPROACH
##			if($num_bootstraps > 1) {
##				#MKDIR, CHDIR?
##				mkdir("$product_name\_bootstrap_temp_files");
##				chdir("$product_name\_bootstrap_temp_files");
##	
##				# TEST BOOTSTRAPPING WORKED
##				
##				#my $RAND_COMP_NUM = int(rand($THIS_CODON_NUM_COMPS));
##				#my $num_bootstrap_replicates = 100;
##				
##				# PARALLELIZE IT
##				my $procs = 60; # number to do at once, cores to assign simultaneously, machine-limited. CUVIER has 80
##				my $pm = Parallel::ForkManager->new($procs); # num_parallel_procs
##				
##				for(my $bootstrap_num = 1; $bootstrap_num <= $num_bootstraps; $bootstrap_num++) {
##					$pm->start and next; # this is IT
##					
##					my $this_replicate_N_diffs_sum = 0;
##					my $this_replicate_S_diffs_sum = 0;
##					my $this_replicate_N_sites_sum = 0;
##					my $this_replicate_S_sites_sum = 0;
##					
##					
##					#print "\ncodon_index\tgroupi\tgroupj\tnum_pw_comps\trand_comparison_index\tsome_N_diffs\tsome_S_diffs\tsome_N_sites\tsome_S_sites\n";
##					
##					for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) { # for each codon in product
##						my $codon_num = $codon_index+1;
##		
##		#				for (my $gi_seq_index = 0; $gi_seq_index < $num_seqs_per_codon_groupi; $gi_seq_index++) { # for each gi sequence at this codon
##		
##		#					for (my $gj_seq_index = 0; $gj_seq_index < $num_seqs_per_codon_groupj; $gj_seq_index++) { # for each gj sequence at this codon
##								
##							my $num_pw_comps = max(keys %{$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}}); # started at 0, up to nC2
##								
##								
##								for(my $comparison_num = 0; $comparison_num <= $num_pw_comps; $comparison_num++) {	
##									srand();
##									my $rand_comparison_index = int(rand($num_pw_comps));
##									
##									$this_replicate_N_diffs_sum += ($gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{N_diffs} / ($num_pw_comps + 1));
##									$this_replicate_S_diffs_sum += ($gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{S_diffs} / ($num_pw_comps + 1));
##									$this_replicate_N_sites_sum += ($gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{N_sites} / ($num_pw_comps + 1));
##									$this_replicate_S_sites_sum += ($gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{S_sites} / ($num_pw_comps + 1));
##									
##									#print THIS_BOOTSTRAP_TEMP_FILE "$codon_index\t$i\t$j\t$num_pw_comps\t$rand_comparison_index\t". 
##									#	$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{N_diffs} ."\t".
##									#	$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{S_diffs} ."\t".
##									#	$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{N_sites} ."\t".
##									#	$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$rand_comparison_index}->{S_sites} ."\n";
##									
##									
##									
##									#$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{N_diffs} = $N_diffs; # curr_comp_num starts at 0
##									#$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{S_diffs} = $S_diffs;
##									#$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{N_sites} = $mean_comp_N_sites;
##									#$gi_gj_codon_comp_stat_hh{$i}->{$j}->{$codon_index}->{$curr_comp_num}->{S_sites} = $mean_comp_S_sites;
##								
##								}
##			
##		#					}
##		#				}
##					}
##					
##					open(THIS_BOOTSTRAP_TEMP_FILE,">>temp\_between\_group\_bootstrap\_$bootstrap_num\.txt");
##					#PRINT
##					print THIS_BOOTSTRAP_TEMP_FILE "$bootstrap_num\t$i\t$j\t".
##						"$this_replicate_N_diffs_sum\t$this_replicate_S_diffs_sum\t".
##						"$this_replicate_N_sites_sum\t$this_replicate_S_sites_sum\n";
##					close THIS_BOOTSTRAP_TEMP_FILE;
##					
##					$pm->finish; # special name
##				}
##				
##				$pm->wait_all_children; # special name, methods within module
##				
##				chdir("..");
##			} # end case in which BOOTSTRAPS are calculated by storing all pw comps






			
			
			# Calculate overall values for the product
			my $between_product_dN;
			if($between_product_N_sites_sum > 0) {
				$between_product_dN = $between_product_N_diffs_sum / $between_product_N_sites_sum;
			} else {
				$between_product_dN = '*';
			}
			
			my $between_product_dS;
			if($between_product_S_sites_sum > 0) {
				$between_product_dS = $between_product_S_diffs_sum / $between_product_S_sites_sum;
			} else {
				$between_product_dS = '*';
			}
			
			my $between_product_w;
			if ($between_product_dS > 0 && $between_product_dS ne '*') {
				$between_product_w = $between_product_dN / $between_product_dS;
			} else {
				$between_product_w = '*';
			}
			
			my $between_product_dN_minus_dS;
			if($between_product_dN >= 0 && $between_product_dS >= 0) {
				$between_product_dN_minus_dS = $between_product_dN - $between_product_dS;
			}
			
			$between_total_N_sites += $between_product_N_sites_sum;
			$between_total_S_sites += $between_product_S_sites_sum;
			$between_total_N_diffs += $between_product_N_diffs_sum;
			$between_total_S_diffs += $between_product_S_diffs_sum;	
			
##			# Print product totals
##			my $product_piN;
##			my $product_piS;
##			my $product_piN_over_piS;
##			
##			if($between_product_N_sites_sum > 0) {
##				$product_piN = $between_product_N_diffs_sum / $between_product_N_sites_sum;
##			} else {
##				$product_piN = '*';
##			}
##			
##			if($between_product_S_sites_sum > 0) {
##				$product_piS = $between_product_S_diffs_sum / $between_product_S_sites_sum;
##			} else {
##				$product_piS = '*';
##			}
##			
##			my $product_piN_minus_piS = $product_piN - $product_piS;
##			
##			if($product_piS > 0) {
##				$product_piN_over_piS = $product_piN / $product_piS;
##			} else {
##				$product_piN_over_piS = '*';
##			}
			
#			my $out_line = "$product_name\t$between_product_N_sites_sum\t".
#				"$between_product_S_sites_sum\t$between_product_N_diffs_sum\t$between_product_S_diffs_sum\t".
#				"$product_piN\t$product_piS\t$product_piN_minus_piS\t$product_piN_over_piS";
			
#			if($num_bootstraps > 1) {
#				my $z_value;
#			
#				if($SE_piN_minus_piS > 0) {
#					$z_value = $product_piN_minus_piS / $SE_piN_minus_piS; 
#				} else {
#					$z_value = '*';
#				} 
#				
#				$out_line .= "\t$SE_piN_minus_piS\t$z_value\n"; #####
#				$out_line .= "\n";
#			} else {
#				$out_line .= "\n";
#			}
			
#			print OUTFILE "$out_line";
#			print "ACTUAL:\n$out_line\n\n";
			
##			$total_N_sites += $between_product_N_sites_sum;
##			$total_S_sites += $between_product_S_sites_sum;
##			$total_N_diffs += $between_product_N_diffs_sum;
##			$total_S_diffs += $between_product_S_diffs_sum;
			
			# Print product totals
			#my $out_line = "$i\t$j\t$product_name\t$between_product_N_sites_sum\t".
			my $out_line = "ANALYSIS\t$product_name\t" .
				"$group_name_i\t$group_name_j\t$between_product_N_sites_sum\t" .
				"$between_product_S_sites_sum\t$between_product_N_diffs_sum\t$between_product_S_diffs_sum\t" .
				"$between_product_dN\t$between_product_dS\t$between_product_dN_minus_dS\t$between_product_w";
			
			print OUTFILE "$out_line\n";
			#print "$out_line";
			
		} # end group i vs. j
	} # end group i; move on to i+1
} # end this product

close OUTFILE;

print "\nTotals:\nN_sites: $between_total_N_sites\nS_sites: $between_total_S_sites\n".
	"N_diffs: $between_total_N_diffs\nS_diffs: $between_total_S_diffs\n\n";






# PROCESS BOOTSTRAPS: same for storage and/or pseudosample approaches
#if($num_bootstraps > 1) {
#	open(BOOTSTRAP_RESULTS_FILE,">>between\_group\_bootstrap\_results\.txt");
#	print BOOTSTRAP_RESULTS_FILE "product\treplicate\tgroup1\tgroup2\tN_diffs\tS_diffs\tN_sites\tS_sites\tdN\tdS\n";
#	close BOOTSTRAP_RESULTS_FILE;
#	my %product_bootstrap_results;
#	my @bootstrap_directories = glob "*_bootstrap_temp_files";
#	#print "\nMy bootstrap directories are: @bootstrap_directories\n\n";
#	foreach(@bootstrap_directories) { # one directory for each product
#		my $directory_name = $_;
#		my $product_name = $directory_name;
#		$product_name =~ s/_bootstrap_temp_files//;
#		chdir("$directory_name");
#		
#		my @bootstrap_files = glob "*";
#		
#		my @bootstrap_num_array;
#		my @group_1_id_array;
#		my @group_2_id_array;
#		my @N_diffs_array;
#		my @S_diffs_array;
#		my @N_sites_array;
#		my @S_sites_array;
#		
#		foreach(@bootstrap_files) {
#			my $bootstrap_file = $_;
#			
#			open (CURR_BOOT_FILE, $bootstrap_file) or die "\n## Cannot open $bootstrap_file. TERMINATED.\n\n";
#			while (<CURR_BOOT_FILE>) {
#				chomp;
#				my $this_line = $_;
#				my @line_arr = split(/\t/,$_,-1);
#				push(@bootstrap_num_array,$line_arr[0]);
#				push(@group_1_id_array,$line_arr[1]);
#				push(@group_2_id_array,$line_arr[2]);
#				push(@N_diffs_array,$line_arr[3]);
#				push(@S_diffs_array,$line_arr[4]);
#				push(@N_sites_array,$line_arr[5]);
#				push(@S_sites_array,$line_arr[6]);
#				
#			}
#			close CURR_BOOT_FILE;
#			unlink $bootstrap_file;
#		}
#		
##		print "\n\nFor $product_name, my\nN_diffs_array: @N_diffs_array\nS_diffs_array: @S_diffs_array\n".
##			"N_sites_array: @N_sites_array\nS_sites_array: @S_sites_array\n\n";
#		
#		chdir('..');
#		rmdir("$directory_name");
#		
#		open(BOOTSTRAP_RESULTS_FILE,">>between\_group\_bootstrap\_results\.txt");
#		for(my $index = 0; $index < scalar(@N_diffs_array); $index++) {
#			my $bootstrap_number = $index+1;
#			
#			my $this_replicate_dN;
#			my $this_replicate_dS;
#			
#			if($N_sites_array[$index] > 0) {
#				$this_replicate_dN = $N_diffs_array[$index] / $N_sites_array[$index]; 
#			} else {
#				$this_replicate_dN = '*';
#			}
#			
#			if($S_sites_array[$index] > 0) {
#				$this_replicate_dS = $S_diffs_array[$index] / $S_sites_array[$index]; 
#			} else {
#				$this_replicate_dS = '*';
#			} 
#			
#			print BOOTSTRAP_RESULTS_FILE "$product_name\t".
#				"$bootstrap_num_array[$index]\t".
#				"$group_1_id_array[$index]\t$group_2_id_array[$index]\t".
#				"$N_diffs_array[$index]\t$S_diffs_array[$index]\t".
#				"$N_sites_array[$index]\t$S_sites_array[$index]\t".
#				"$this_replicate_dN\t$this_replicate_dS\n";
#		}
#		close BOOTSTRAP_RESULTS_FILE;
#		
#	}
#} # end processing bootstraps

&end_the_program;


#########################################################################################
###############################                         #################################
###############################       SUBROUTINES       #################################
###############################                         #################################
#########################################################################################


#########################################################################################
# Obtains all file names in current directory ending in .fa and/or .fasta
sub get_fasta_file_names { 
	my @fasta_file_names = glob "*.fa";
	my @other_fasta_file_names;

	@other_fasta_file_names = glob "*.fasta";
	push (@fasta_file_names,@other_fasta_file_names);

	if (scalar(@fasta_file_names) == 0) {
		die "\n\n## WARNING: There are no .fa or .fasta files. SNPGenie terminated.\n\n";
	}
	return 	@fasta_file_names;
}

#########################################################################################
sub get_product_names_from_gtf {
	my ($cds_file) = @_;
	#print "\n\n$cds_file\n\n";
	my %products_hash;
	open (CURRINFILE, $cds_file) or die "\n## Cannot open $cds_file. TERMINATED.\n\n";
	while (<CURRINFILE>) {
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			if($_ =~/gene_id \"gene\:([\w\s\.\-\:']+)\"/) { # transcript_id not a problem
				$products_hash{$1} = 1;
			} elsif($_ =~ /gene_id \"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
			} elsif($_ =~/gene_id \"([\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
			} else {
				die "\n\n## WARNING: CDS annotation(s) in $gtf_file does not have a ".
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
	
	open (CURRINFILE, $gtf_file);
	while (<CURRINFILE>) { # go through the GTF file
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			my $product_present = 0;
			my $this_line = $_;
			my $this_start;
			my $this_stop;
			
			# incomplete segments
			if ($_ =~/gene_id "$product";/) {
				if ($_ =~/CDS\t(\d+)\t(\d+)/) {
					$product_present = 1;
					$this_start = $1;
					$this_stop = $2;
				}
			} elsif($_ =~/gene_id "gene\:$product";/) {
				if ($_ =~/CDS\t(\d+)\t(\d+)/) {
					$product_present = 1;
					$this_start = $1;
					$this_stop = $2;
				}
			}
			
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
# Awesome subroutine that returns the number of differences between two codons,
# averaged over all paths
# We have altered this subroutine to return 0 differences for STOP codons,
# but synonymous differences if it's one STOP compared to another.
sub return_avg_diffs {
	my ($codon1,$codon2) = @_;
	
	my %all_diffs_hh = (
		'AAA' => {
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 0,
				'S' => 1
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 0
			},
			'ACC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 1,
				'S' => 0
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 1,
				'S' => 1
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 1,
				'S' => 1
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 0
			},
			'CCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 2,
				'S' => 0
			},
			'CTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 1,
				'S' => 1
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 3,
				'S' => 0
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 3,
				'S' => 0
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTT' => {
				'N' => 2.75,
				'S' => 0.25
			},
		},
		'AAC' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 0,
				'S' => 1
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 1,
				'S' => 0
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 1
			},
			'CCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 1
			},
			'GCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 3,
				'S' => 0
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 3,
				'S' => 0
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'AAG' => {
			'AAA' => {
				'N' => 0,
				'S' => 1
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 1,
				'S' => 0
			},
			'ACT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 1,
				'S' => 1
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 1,
				'S' => 1
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCG' => {
				'N' => 2,
				'S' => 0
			},
			'CCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CTG' => {
				'N' => 2,
				'S' => 0
			},
			'CTT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GAA' => {
				'N' => 1,
				'S' => 1
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 3,
				'S' => 0
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 3,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTC' => {
				'N' => 3,
				'S' => 0
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 3,
				'S' => 0
			},
		},
		'AAT' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAC' => {
				'N' => 0,
				'S' => 1
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 0
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 1,
				'S' => 1
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 1,
				'S' => 1
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 3,
				'S' => 0
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 3,
				'S' => 0
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'ACA' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 1,
				'S' => 1
			},
			'AAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 0,
				'S' => 1
			},
			'ACG' => {
				'N' => 0,
				'S' => 1
			},
			'ACT' => {
				'N' => 0,
				'S' => 1
			},
			'AGA' => {
				'N' => 1,
				'S' => 0
			},
			'AGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 1,
				'S' => 1
			},
			'AGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 1
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 0
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 0
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 0
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1,
				'S' => 0
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
		},
		'ACC' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 1,
				'S' => 1
			},
			'ACA' => {
				'N' => 0,
				'S' => 1
			},
			'ACG' => {
				'N' => 0,
				'S' => 1
			},
			'ACT' => {
				'N' => 0,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 1,
				'S' => 1
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 0
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 0
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 0
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'ACG' => {
			'AAA' => {
				'N' => 1,
				'S' => 1
			},
			'AAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 0,
				'S' => 1
			},
			'ACC' => {
				'N' => 0,
				'S' => 1
			},
			'ACT' => {
				'N' => 0,
				'S' => 1
			},
			'AGA' => {
				'N' => 1,
				'S' => 1
			},
			'AGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 2,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 0
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTG' => {
				'N' => 2,
				'S' => 0
			},
			'CTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 0
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 0
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
		},
		'ACT' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 1,
				'S' => 1
			},
			'AAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 0,
				'S' => 1
			},
			'ACC' => {
				'N' => 0,
				'S' => 1
			},
			'ACG' => {
				'N' => 0,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 1,
				'S' => 1
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 0
			},
			'CGA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 0
			},
			'GGA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'AGA' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 1
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 0
			},
			'ACC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 0,
				'S' => 1
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 0,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 0,
				'S' => 2
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.66666666666667,
				'S' => 1.33333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGA' => {
				'N' => 1,
				'S' => 0
			},
			'GGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 3,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 3,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 1,
				'S' => 1
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTT' => {
				'N' => 2.75,
				'S' => 0.25
			},
		},
		'AGC' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 1
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 1,
				'S' => 0
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 0,
				'S' => 1
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 0
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGC' => {
				'N' => 1,
				'S' => 0
			},
			'GGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 3,
				'S' => 0
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'AGG' => {
			'AAA' => {
				'N' => 1,
				'S' => 1
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 1,
				'S' => 0
			},
			'ACT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 0,
				'S' => 1
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 0,
				'S' => 2
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 0,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 1.66666666666667,
				'S' => 1.33333333333333
			},
			'CTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 1,
				'S' => 0
			},
			'GGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 3,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 3,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTC' => {
				'N' => 3,
				'S' => 0
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 3,
				'S' => 0
			},
		},
		'AGT' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 1
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 0
			},
			'AGA' => {
				'N' => 1,
				'S' => 0
			},
			'AGC' => {
				'N' => 0,
				'S' => 1
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 0
			},
			'CTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 1,
				'S' => 0
			},
			'GTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 3,
				'S' => 0
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'ATA' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 1,
				'S' => 0
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 1,
				'S' => 0
			},
			'AGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 0,
				'S' => 1
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 0,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 2,
				'S' => 0
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 0
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTA' => {
				'N' => 1,
				'S' => 0
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TCT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTA' => {
				'N' => 1,
				'S' => 0
			},
			'TTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
		},
		'ATC' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 1
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 0
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 0,
				'S' => 1
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 0,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 0
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 0
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 3,
				'S' => 0
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'ATG' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 1,
				'S' => 0
			},
			'ACT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CAC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'CCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCG' => {
				'N' => 2,
				'S' => 0
			},
			'CCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTG' => {
				'N' => 1,
				'S' => 0
			},
			'CTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GAC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 1,
				'S' => 0
			},
			'GTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 3,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 3,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TGC' => {
				'N' => 3,
				'S' => 0
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 3,
				'S' => 0
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 1,
				'S' => 0
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'ATT' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 1,
				'S' => 1
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1,
				'S' => 0
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 0,
				'S' => 1
			},
			'ATC' => {
				'N' => 0,
				'S' => 1
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 1,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 3,
				'S' => 0
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'CAA' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 1
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 0
			},
			'ACC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 2,
				'S' => 1
			},
			'ACT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'AGT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'ATA' => {
				'N' => 2,
				'S' => 0
			},
			'ATC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 0,
				'S' => 1
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 0
			},
			'CCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGA' => {
				'N' => 1,
				'S' => 0
			},
			'CGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 1,
				'S' => 0
			},
			'CTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 1,
				'S' => 1
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1,
				'S' => 2
			},
			'TTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
		},
		'CAC' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 1
			},
			'ACA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 0,
				'S' => 1
			},
			'CCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 1,
				'S' => 0
			},
			'CCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 1,
				'S' => 0
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 1,
				'S' => 0
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 1
			},
			'GCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'CAG' => {
			'AAA' => {
				'N' => 1,
				'S' => 1
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'AGC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'ATA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'CAA' => {
				'N' => 0,
				'S' => 1
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCG' => {
				'N' => 1,
				'S' => 0
			},
			'CCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 1,
				'S' => 0
			},
			'CGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTG' => {
				'N' => 1,
				'S' => 0
			},
			'CTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAA' => {
				'N' => 1,
				'S' => 1
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTA' => {
				'N' => 1,
				'S' => 2
			},
			'TTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1,
				'S' => 1
			},
			'TTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
		},
		'CAT' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 1
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAC' => {
				'N' => 0,
				'S' => 1
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 1,
				'S' => 0
			},
			'CGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 1,
				'S' => 0
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 1,
				'S' => 0
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 1,
				'S' => 1
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'CCA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 2,
				'S' => 1
			},
			'AAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 1,
				'S' => 0
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'AGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATA' => {
				'N' => 2,
				'S' => 0
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1,
				'S' => 1
			},
			'CAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 0,
				'S' => 1
			},
			'CCG' => {
				'N' => 0,
				'S' => 1
			},
			'CCT' => {
				'N' => 0,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 0
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 0
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 0
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1,
				'S' => 0
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'TTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
		},
		'CCC' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 0
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2,
				'S' => 1
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 2,
				'S' => 1
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 1,
				'S' => 1
			},
			'CCA' => {
				'N' => 0,
				'S' => 1
			},
			'CCG' => {
				'N' => 0,
				'S' => 1
			},
			'CCT' => {
				'N' => 0,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 0
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 0
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 0
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 0
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2,
				'S' => 1
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'CCG' => {
			'AAA' => {
				'N' => 2,
				'S' => 1
			},
			'AAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 0
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'AGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CAA' => {
				'N' => 1,
				'S' => 1
			},
			'CAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 0,
				'S' => 1
			},
			'CCC' => {
				'N' => 0,
				'S' => 1
			},
			'CCT' => {
				'N' => 0,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 0
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 0
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 0
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 0
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'TTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
		},
		'CCT' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 0
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2,
				'S' => 1
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2,
				'S' => 1
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 1
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 0,
				'S' => 1
			},
			'CCC' => {
				'N' => 0,
				'S' => 1
			},
			'CCG' => {
				'N' => 0,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 0
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 0
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2,
				'S' => 1
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'CGA' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AAG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'AAT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'ACT' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'AGA' => {
				'N' => 0,
				'S' => 1
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 0,
				'S' => 2
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'ATG' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'ATT' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1,
				'S' => 1
			},
			'CAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 0
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 0,
				'S' => 1
			},
			'CGG' => {
				'N' => 0,
				'S' => 1
			},
			'CGT' => {
				'N' => 0,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 0
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 0
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 1,
				'S' => 1
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 1.25,
				'S' => 1.75
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'CGC' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2,
				'S' => 1
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 1,
				'S' => 1
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 1,
				'S' => 1
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 2,
				'S' => 1
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 1,
				'S' => 1
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 0
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 0,
				'S' => 1
			},
			'CGG' => {
				'N' => 0,
				'S' => 1
			},
			'CGT' => {
				'N' => 0,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 0
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 0
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'CGG' => {
			'AAA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'AAC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ACA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'ACC' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'ACG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'AGA' => {
				'N' => 0,
				'S' => 2
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 0,
				'S' => 1
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'ATC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CAA' => {
				'N' => 1,
				'S' => 1
			},
			'CAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 0
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 0,
				'S' => 1
			},
			'CGC' => {
				'N' => 0,
				'S' => 1
			},
			'CGT' => {
				'N' => 0,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 0
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 0
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTA' => {
				'N' => 1.25,
				'S' => 1.75
			},
			'TTC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
		},
		'CGT' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2,
				'S' => 1
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 1,
				'S' => 1
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 1,
				'S' => 1
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 2,
				'S' => 1
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 1
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 0
			},
			'CGA' => {
				'N' => 0,
				'S' => 1
			},
			'CGC' => {
				'N' => 0,
				'S' => 1
			},
			'CGG' => {
				'N' => 0,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 0
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'CTA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 2,
				'S' => 0
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGG' => {
				'N' => 1.66666666666667,
				'S' => 1.33333333333333
			},
			'AGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1,
				'S' => 1
			},
			'CAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 0
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 0
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 0,
				'S' => 1
			},
			'CTG' => {
				'N' => 0,
				'S' => 1
			},
			'CTT' => {
				'N' => 0,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 0
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 1.5
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCC' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'TCG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'TCT' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 0,
				'S' => 1
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 0,
				'S' => 2
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'CTC' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 1,
				'S' => 1
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 1,
				'S' => 1
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 0
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 0
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 0,
				'S' => 1
			},
			'CTG' => {
				'N' => 0,
				'S' => 1
			},
			'CTT' => {
				'N' => 0,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 0
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 1,
				'S' => 1
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'CTG' => {
			'AAA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AAC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'ACA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGA' => {
				'N' => 1.66666666666667,
				'S' => 1.33333333333333
			},
			'AGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 1,
				'S' => 1
			},
			'CAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 0
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 0
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 0,
				'S' => 1
			},
			'CTC' => {
				'N' => 0,
				'S' => 1
			},
			'CTT' => {
				'N' => 0,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 0
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 1.5
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'TCC' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'TCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.83333333333333
				'S' => 0 # 1.16666666666667
			},
			'TGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TTA' => {
				'N' => 0,
				'S' => 2
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 0,
				'S' => 1
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'CTT' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 1,
				'S' => 1
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 1,
				'S' => 1
			},
			'CAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 0
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 0
			},
			'CTA' => {
				'N' => 0,
				'S' => 1
			},
			'CTC' => {
				'N' => 0,
				'S' => 1
			},
			'CTG' => {
				'N' => 0,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 1,
				'S' => 1
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'GAA' => {
			'AAA' => {
				'N' => 1,
				'S' => 0
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 1
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 0
			},
			'ACC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 2,
				'S' => 1
			},
			'ACT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'AGG' => {
				'N' => 2,
				'S' => 1
			},
			'AGT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATA' => {
				'N' => 2,
				'S' => 0
			},
			'ATC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 1,
				'S' => 0
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 1,
				'S' => 1
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 0
			},
			'CCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGA' => {
				'N' => 2,
				'S' => 0
			},
			'CGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 2,
				'S' => 0
			},
			'CTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 0,
				'S' => 1
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 0
			},
			'GCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGA' => {
				'N' => 1,
				'S' => 0
			},
			'GGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 1,
				'S' => 0
			},
			'GTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTG' => {
				'N' => 2,
				'S' => 1
			},
			'TTT' => {
				'N' => 2.75,
				'S' => 0.25
			},
		},
		'GAC' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 1
			},
			'ACA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 1
			},
			'CCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 0,
				'S' => 1
			},
			'GCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 1,
				'S' => 0
			},
			'GCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGC' => {
				'N' => 1,
				'S' => 0
			},
			'GGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 1,
				'S' => 0
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'GAG' => {
			'AAA' => {
				'N' => 1,
				'S' => 1
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 1,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATC' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'CAA' => {
				'N' => 1,
				'S' => 1
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 1,
				'S' => 0
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCG' => {
				'N' => 2,
				'S' => 0
			},
			'CCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 2,
				'S' => 0
			},
			'CGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTG' => {
				'N' => 2,
				'S' => 0
			},
			'CTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAA' => {
				'N' => 0,
				'S' => 1
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCG' => {
				'N' => 1,
				'S' => 0
			},
			'GCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 1,
				'S' => 0
			},
			'GGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 1,
				'S' => 0
			},
			'GTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 2.75,
				'S' => 0.25
			},
		},
		'GAT' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 1
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.83333333333333,
				'S' => 0.166666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 1,
				'S' => 1
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAC' => {
				'N' => 0,
				'S' => 1
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCT' => {
				'N' => 1,
				'S' => 0
			},
			'GGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 1,
				'S' => 0
			},
			'GTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'GCA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 2,
				'S' => 1
			},
			'AAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 1,
				'S' => 0
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGG' => {
				'N' => 2,
				'S' => 1
			},
			'AGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATA' => {
				'N' => 2,
				'S' => 0
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 1
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 0
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 0
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 0
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 1,
				'S' => 1
			},
			'GAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCC' => {
				'N' => 0,
				'S' => 1
			},
			'GCG' => {
				'N' => 0,
				'S' => 1
			},
			'GCT' => {
				'N' => 0,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 0
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 0
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1,
				'S' => 0
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTG' => {
				'N' => 2,
				'S' => 1
			},
			'TTT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
		},
		'GCC' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 0
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 2,
				'S' => 1
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 0
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 1,
				'S' => 1
			},
			'GCA' => {
				'N' => 0,
				'S' => 1
			},
			'GCG' => {
				'N' => 0,
				'S' => 1
			},
			'GCT' => {
				'N' => 0,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 0
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 0
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 0
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'GCG' => {
			'AAA' => {
				'N' => 2,
				'S' => 1
			},
			'AAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 0
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CAA' => {
				'N' => 2,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 0
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 0
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2,
				'S' => 0
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 1
			},
			'GAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 0,
				'S' => 1
			},
			'GCC' => {
				'N' => 0,
				'S' => 1
			},
			'GCT' => {
				'N' => 0,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 0
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 0
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 0
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
		},
		'GCT' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 0
			},
			'AGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2,
				'S' => 1
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 0
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 1
			},
			'GAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 0,
				'S' => 1
			},
			'GCC' => {
				'N' => 0,
				'S' => 1
			},
			'GCG' => {
				'N' => 0,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 0
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'GGA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAG' => {
				'N' => 2,
				'S' => 1
			},
			'AAT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'ACA' => {
				'N' => 2,
				'S' => 0
			},
			'ACC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACG' => {
				'N' => 2,
				'S' => 1
			},
			'ACT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGA' => {
				'N' => 1,
				'S' => 0
			},
			'AGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 1,
				'S' => 1
			},
			'AGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATA' => {
				'N' => 2,
				'S' => 0
			},
			'ATC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 1
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 2,
				'S' => 0
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 0
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 0
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 1,
				'S' => 1
			},
			'GAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 0
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 0,
				'S' => 1
			},
			'GGG' => {
				'N' => 0,
				'S' => 1
			},
			'GGT' => {
				'N' => 0,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 0
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 1,
				'S' => 1
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTG' => {
				'N' => 2,
				'S' => 1
			},
			'TTT' => {
				'N' => 2.25,
				'S' => 0.75
			},
		},
		'GGC' => {
			'AAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 0
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 1,
				'S' => 1
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 0
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 0,
				'S' => 1
			},
			'GGG' => {
				'N' => 0,
				'S' => 1
			},
			'GGT' => {
				'N' => 0,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 0
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 2,
				'S' => 1
			},
		},
		'GGG' => {
			'AAA' => {
				'N' => 2,
				'S' => 1
			},
			'AAC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGA' => {
				'N' => 1,
				'S' => 1
			},
			'AGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 2,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2,
				'S' => 0
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 0
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2,
				'S' => 0
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 1
			},
			'GAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 0
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 0,
				'S' => 1
			},
			'GGC' => {
				'N' => 0,
				'S' => 1
			},
			'GGT' => {
				'N' => 0,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 0
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTA' => {
				'N' => 2,
				'S' => 1
			},
			'TTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
		},
		'GGT' => {
			'AAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1,
				'S' => 1
			},
			'CGT' => {
				'N' => 1,
				'S' => 0
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2,
				'S' => 1
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 1
			},
			'GAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 0
			},
			'GGA' => {
				'N' => 0,
				'S' => 1
			},
			'GGC' => {
				'N' => 0,
				'S' => 1
			},
			'GGG' => {
				'N' => 0,
				'S' => 1
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 2,
				'S' => 1
			},
			'TTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'GTA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACA' => {
				'N' => 2,
				'S' => 0
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 1
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 2,
				'S' => 0
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 0
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 0
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 0
			},
			'GAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 1,
				'S' => 1
			},
			'GAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 0
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 0
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 0,
				'S' => 1
			},
			'GTG' => {
				'N' => 0,
				'S' => 1
			},
			'GTT' => {
				'N' => 0,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCA' => {
				'N' => 2,
				'S' => 0
			},
			'TCC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCG' => {
				'N' => 2,
				'S' => 1
			},
			'TCT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TGG' => {
				'N' => 2,
				'S' => 1
			},
			'TGT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'TTA' => {
				'N' => 1,
				'S' => 0
			},
			'TTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1,
				'S' => 1
			},
			'TTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
		},
		'GTC' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 1,
				'S' => 1
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 0
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 1,
				'S' => 1
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 0
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 0
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 0,
				'S' => 1
			},
			'GTG' => {
				'N' => 0,
				'S' => 1
			},
			'GTT' => {
				'N' => 0,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAT' => {
				'N' => 2,
				'S' => 1
			},
			'TCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCC' => {
				'N' => 2,
				'S' => 0
			},
			'TCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCT' => {
				'N' => 2,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 2,
				'S' => 1
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'GTG' => {
			'AAA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AAC' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'ACA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'AGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 2,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2,
				'S' => 0
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 0
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 0
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 1,
				'S' => 1
			},
			'GAC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAG' => {
				'N' => 1,
				'S' => 0
			},
			'GAT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 0
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 0
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 0,
				'S' => 1
			},
			'GTC' => {
				'N' => 0,
				'S' => 1
			},
			'GTT' => {
				'N' => 0,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TCA' => {
				'N' => 2,
				'S' => 1
			},
			'TCC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCG' => {
				'N' => 2,
				'S' => 0
			},
			'TCT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'TGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1,
				'S' => 0
			},
			'TTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
		},
		'GTT' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 1,
				'S' => 1
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2,
				'S' => 1
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 0
			},
			'GAA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 1,
				'S' => 1
			},
			'GAG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 0
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1,
				'S' => 1
			},
			'GGT' => {
				'N' => 1,
				'S' => 0
			},
			'GTA' => {
				'N' => 0,
				'S' => 1
			},
			'GTC' => {
				'N' => 0,
				'S' => 1
			},
			'GTG' => {
				'N' => 0,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAC' => {
				'N' => 2,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCC' => {
				'N' => 2,
				'S' => 1
			},
			'TCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TCT' => {
				'N' => 2,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 1
			},
			'TGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'TAA' => { # altered because STOP
			'AAA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'AAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'AAG' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'AAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ACA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ACC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'ACG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'ACT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'AGA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'AGC' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'AGG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'AGT' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'ATA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ATC' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'ATG' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'ATT' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'CAA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'CAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CAG' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'CAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CCA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CCC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'CCG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'CCT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'CGA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CGC' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'CGG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'CGT' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'CTA' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'CTC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'CTG' => {
				'N' => 0, # 1.5
				'S' => 0 # 1.5
			},
			'CTT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GAA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'GAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GAG' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'GAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GCA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GCC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GCG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'GCT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GGA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GGC' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'GGG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'GGT' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'GTA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GTC' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'GTG' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'GTT' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAC' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 0
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TCA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TCC' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TCG' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TCT' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 0
				'S' => 0 # 1
			},
			'TGC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGG' => {
				'N' => '*',
				'S' => '*'
			},
			'TGT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TTA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TTC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TTG' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TTT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
		},
		'TAC' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 1
			},
			'ACA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 3,
				'S' => 0
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 3,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 3,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 1,
				'S' => 0
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 1
			},
			'CCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 1,
				'S' => 0
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 1
			},
			'GCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 0,
				'S' => 1
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 0
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'TAG' => { # altered because STOP
			'AAA' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'AAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'AAG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'AAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ACA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'ACC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'ACG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ACT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'AGA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'AGC' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'AGG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'AGT' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'ATA' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'ATC' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'ATG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ATT' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'CAA' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'CAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CAG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'CAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CCA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'CCC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'CCG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CCT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'CGA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'CGC' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'CGG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CGT' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'CTA' => {
				'N' => 0, # 1.5
				'S' => 0 # 1.5
			},
			'CTC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'CTG' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'CTT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GAA' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'GAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GAG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'GAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GCA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'GCC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GCG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GCT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GGA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'GGC' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'GGG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GGT' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'GTA' => {
				'N' => 0, # 2
				'S' => 0 # 1
			},
			'GTC' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'GTG' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GTT' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 0
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TCA' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TCC' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TCG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TCT' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TTA' => {
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TTC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TTG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TTT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
		},
		'TAT' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 1,
				'S' => 1
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 1,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 3,
				'S' => 0
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 3,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 3,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 1,
				'S' => 1
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 1,
				'S' => 0
			},
			'CCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 1,
				'S' => 1
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 1,
				'S' => 0
			},
			'GCA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 0,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'TCA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'AAG' => {
				'N' => 2,
				'S' => 1
			},
			'AAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACA' => {
				'N' => 1,
				'S' => 0
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'AGG' => {
				'N' => 2,
				'S' => 1
			},
			'AGT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ATA' => {
				'N' => 2,
				'S' => 0
			},
			'ATC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CAA' => {
				'N' => 2,
				'S' => 0
			},
			'CAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CAG' => {
				'N' => 2,
				'S' => 1
			},
			'CAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCA' => {
				'N' => 1,
				'S' => 0
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 0
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2,
				'S' => 1
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GCA' => {
				'N' => 1,
				'S' => 0
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2,
				'S' => 0
			},
			'GTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTG' => {
				'N' => 2,
				'S' => 1
			},
			'GTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 0,
				'S' => 1
			},
			'TCG' => {
				'N' => 0,
				'S' => 1
			},
			'TCT' => {
				'N' => 0,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 1,
				'S' => 1
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 1,
				'S' => 0
			},
			'TTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1,
				'S' => 1
			},
			'TTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
		},
		'TCC' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 0
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 0
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 0
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 1,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 0,
				'S' => 1
			},
			'TCG' => {
				'N' => 0,
				'S' => 1
			},
			'TCT' => {
				'N' => 0,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'TCG' => {
			'AAA' => {
				'N' => 2,
				'S' => 1
			},
			'AAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 0
			},
			'ACT' => {
				'N' => 1,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 1
			},
			'AGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ATC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 2,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 0
			},
			'CCT' => {
				'N' => 1,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CGG' => {
				'N' => 2,
				'S' => 0
			},
			'CGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 0
			},
			'GCT' => {
				'N' => 1,
				'S' => 1
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 0,
				'S' => 1
			},
			'TCC' => {
				'N' => 0,
				'S' => 1
			},
			'TCT' => {
				'N' => 0,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTG' => {
				'N' => 1,
				'S' => 0
			},
			'TTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
		},
		'TCT' => {
			'AAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 1,
				'S' => 1
			},
			'ACC' => {
				'N' => 1,
				'S' => 1
			},
			'ACG' => {
				'N' => 1,
				'S' => 1
			},
			'ACT' => {
				'N' => 1,
				'S' => 0
			},
			'AGA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 1,
				'S' => 1
			},
			'CCC' => {
				'N' => 1,
				'S' => 1
			},
			'CCG' => {
				'N' => 1,
				'S' => 1
			},
			'CCT' => {
				'N' => 1,
				'S' => 0
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 1.83333333333333,
				'S' => 1.16666666666667
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 1,
				'S' => 1
			},
			'GCC' => {
				'N' => 1,
				'S' => 1
			},
			'GCG' => {
				'N' => 1,
				'S' => 1
			},
			'GCT' => {
				'N' => 1,
				'S' => 0
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TAT' => {
				'N' => 1,
				'S' => 0
			},
			'TCA' => {
				'N' => 0,
				'S' => 1
			},
			'TCC' => {
				'N' => 0,
				'S' => 1
			},
			'TCG' => {
				'N' => 0,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'TGA' => { # altered because STOP
			'AAA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'AAC' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'AAG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'AAT' => {
				'N' => 0, # 3
				'S' => 0 # 0
			},
			'ACA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ACC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'ACG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'ACT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'AGA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'AGC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'AGG' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'AGT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ATA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'ATC' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'ATG' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'ATT' => {
				'N' => 0, # 2.66666666666667
				'S' => 0 # 0.333333333333333
			},
			'CAA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CAC' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'CAG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'CAT' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'CCA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'CCC' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'CCG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'CCT' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'CGA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'CGC' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'CGG' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'CGT' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'CTA' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'CTC' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'CTG' => {
				'N' => 0, # 1.83333333333333
				'S' => 0 # 1.16666666666667
			},
			'CTT' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'GAA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GAC' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'GAG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'GAT' => {
				'N' => 0, # 2.75
				'S' => 0 # 0.25
			},
			'GCA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GCC' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'GCG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'GCT' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'GGA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'GGC' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'GGG' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'GGT' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'GTA' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'GTC' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'GTG' => {
				'N' => 0, # 2.33333333333333
				'S' => 0 # 0.666666666666667
			},
			'GTT' => {
				'N' => 0, # 2.5
				'S' => 0 # 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 0
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TCA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TCC' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TCG' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TCT' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGG' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGT' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TTA' => {
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TTC' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TTG' => {
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TTT' => {
				'N' => 0, # 2
				'S' => 0 # 0
			},
		},
		'TGC' => {
			'AAA' => {
				'N' => 3,
				'S' => 0
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 3,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 1,
				'S' => 0
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 1
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 3,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 0
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 1,
				'S' => 1
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 0
			},
			'CTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 0
			},
			'GGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 1,
				'S' => 1
			},
			'GTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GTC' => {
				'N' => 2,
				'S' => 0
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 0
			},
			'TCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 0,
				'S' => 1
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 1
			},
		},
		'TGG' => {
			'AAA' => {
				'N' => 2,
				'S' => 1
			},
			'AAC' => {
				'N' => 3,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 3,
				'S' => 0
			},
			'ACA' => {
				'N' => 2,
				'S' => 1
			},
			'ACC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 1,
				'S' => 1
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 1,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 3,
				'S' => 0
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 3,
				'S' => 0
			},
			'CAA' => {
				'N' => 2,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'CAG' => {
				'N' => 2,
				'S' => 0
			},
			'CAT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCG' => {
				'N' => 2,
				'S' => 0
			},
			'CCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGG' => {
				'N' => 1,
				'S' => 0
			},
			'CGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CTC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CTT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 1,
				'S' => 0
			},
			'GGT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 2,
				'S' => 1
			},
			'GTC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 2,
				'S' => 0
			},
			'GTT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => '*',
				'S' => '*'
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 1,
				'S' => 0
			},
			'TCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 1,
				'S' => 1
			},
			'TTC' => {
				'N' => 2,
				'S' => 0
			},
			'TTG' => {
				'N' => 1,
				'S' => 0
			},
			'TTT' => {
				'N' => 2,
				'S' => 0
			},
		},
		'TGT' => {
			'AAA' => {
				'N' => 3,
				'S' => 0
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 3,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 1,
				'S' => 1
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 1,
				'S' => 0
			},
			'ATA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 1
			},
			'ATG' => {
				'N' => 3,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2,
				'S' => 1
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 1,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 1,
				'S' => 0
			},
			'CTA' => {
				'N' => 2,
				'S' => 1
			},
			'CTC' => {
				'N' => 2,
				'S' => 1
			},
			'CTG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTT' => {
				'N' => 2,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.66666666666667,
				'S' => 0.333333333333333
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 1,
				'S' => 1
			},
			'GGC' => {
				'N' => 1,
				'S' => 1
			},
			'GGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 1,
				'S' => 0
			},
			'GTA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GTC' => {
				'N' => 2,
				'S' => 1
			},
			'GTG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 2,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 1,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 0,
				'S' => 1
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 1,
				'S' => 1
			},
			'TTG' => {
				'N' => 2,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'TTA' => {
			'AAA' => {
				'N' => 2,
				'S' => 0
			},
			'AAC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'AAG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AAT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'ACA' => {
				'N' => 2,
				'S' => 0
			},
			'ACC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ACG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AGA' => {
				'N' => 2,
				'S' => 0
			},
			'AGC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'AGG' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'AGT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'ATA' => {
				'N' => 1,
				'S' => 0
			},
			'ATC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CAA' => {
				'N' => 1,
				'S' => 1
			},
			'CAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CAG' => {
				'N' => 1,
				'S' => 2
			},
			'CAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 1.25,
				'S' => 1.75
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 0,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 0,
				'S' => 2
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 0
			},
			'GAC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAG' => {
				'N' => 2,
				'S' => 1
			},
			'GAT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GCA' => {
				'N' => 2,
				'S' => 0
			},
			'GCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCG' => {
				'N' => 2,
				'S' => 1
			},
			'GCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGA' => {
				'N' => 2,
				'S' => 0
			},
			'GGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 2,
				'S' => 1
			},
			'GGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 1,
				'S' => 0
			},
			'GTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 1,
				'S' => 1
			},
			'GTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 0
			},
			'TCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 1,
				'S' => 1
			},
			'TCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 1,
				'S' => 1
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 0,
				'S' => 1
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'TTC' => {
			'AAA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'AAC' => {
				'N' => 2,
				'S' => 0
			},
			'AAG' => {
				'N' => 3,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 1
			},
			'ACA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ACC' => {
				'N' => 2,
				'S' => 0
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 1
			},
			'AGA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'AGC' => {
				'N' => 2,
				'S' => 0
			},
			'AGG' => {
				'N' => 3,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 1
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1,
				'S' => 0
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 1
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 0
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 1
			},
			'CCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCC' => {
				'N' => 2,
				'S' => 0
			},
			'CCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 0
			},
			'CGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 1
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 0
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAC' => {
				'N' => 2,
				'S' => 0
			},
			'GAG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAT' => {
				'N' => 2,
				'S' => 1
			},
			'GCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCC' => {
				'N' => 2,
				'S' => 0
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 1
			},
			'GGA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GGC' => {
				'N' => 2,
				'S' => 0
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 1
			},
			'GTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 1,
				'S' => 0
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 1
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 1
			},
			'TCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCC' => {
				'N' => 1,
				'S' => 0
			},
			'TCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 1,
				'S' => 1
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 0
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 1
			},
			'TTA' => {
				'N' => 1,
				'S' => 0
			},
			'TTG' => {
				'N' => 1,
				'S' => 0
			},
			'TTT' => {
				'N' => 0,
				'S' => 1
			},
		},
		'TTG' => {
			'AAA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'AAC' => {
				'N' => 3,
				'S' => 0
			},
			'AAG' => {
				'N' => 2,
				'S' => 0
			},
			'AAT' => {
				'N' => 3,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'ACC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACG' => {
				'N' => 2,
				'S' => 0
			},
			'ACT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'AGA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'AGC' => {
				'N' => 3,
				'S' => 0
			},
			'AGG' => {
				'N' => 2,
				'S' => 0
			},
			'AGT' => {
				'N' => 3,
				'S' => 0
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 2,
				'S' => 0
			},
			'ATG' => {
				'N' => 1,
				'S' => 0
			},
			'ATT' => {
				'N' => 2,
				'S' => 0
			},
			'CAA' => {
				'N' => 1,
				'S' => 2
			},
			'CAC' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CAG' => {
				'N' => 1,
				'S' => 1
			},
			'CAT' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'CCA' => {
				'N' => 1.5,
				'S' => 1.5
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CCT' => {
				'N' => 2,
				'S' => 1
			},
			'CGA' => {
				'N' => 1.25,
				'S' => 1.75
			},
			'CGC' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CGG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'CGT' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CTA' => {
				'N' => 0,
				'S' => 2
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 0,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 1
			},
			'GAA' => {
				'N' => 2,
				'S' => 1
			},
			'GAC' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAG' => {
				'N' => 2,
				'S' => 0
			},
			'GAT' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GCA' => {
				'N' => 2,
				'S' => 1
			},
			'GCC' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCG' => {
				'N' => 2,
				'S' => 0
			},
			'GCT' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GGA' => {
				'N' => 2,
				'S' => 1
			},
			'GGC' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGG' => {
				'N' => 2,
				'S' => 0
			},
			'GGT' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GTA' => {
				'N' => 1,
				'S' => 1
			},
			'GTC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTG' => {
				'N' => 1,
				'S' => 0
			},
			'GTT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 1
			},
			'TAC' => {
				'N' => 2,
				'S' => 0
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 1
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 2,
				'S' => 0
			},
			'TCA' => {
				'N' => 1,
				'S' => 1
			},
			'TCC' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCG' => {
				'N' => 1,
				'S' => 0
			},
			'TCT' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 1.5
				'S' => 0 # 0.5
			},
			'TGC' => {
				'N' => 2,
				'S' => 0
			},
			'TGG' => {
				'N' => 1,
				'S' => 0
			},
			'TGT' => {
				'N' => 2,
				'S' => 0
			},
			'TTA' => {
				'N' => 0,
				'S' => 1
			},
			'TTC' => {
				'N' => 1,
				'S' => 0
			},
			'TTT' => {
				'N' => 1,
				'S' => 0
			},
		},
		'TTT' => {
			'AAA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'AAC' => {
				'N' => 2,
				'S' => 1
			},
			'AAG' => {
				'N' => 3,
				'S' => 0
			},
			'AAT' => {
				'N' => 2,
				'S' => 0
			},
			'ACA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'ACC' => {
				'N' => 2,
				'S' => 1
			},
			'ACG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'ACT' => {
				'N' => 2,
				'S' => 0
			},
			'AGA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'AGC' => {
				'N' => 2,
				'S' => 1
			},
			'AGG' => {
				'N' => 3,
				'S' => 0
			},
			'AGT' => {
				'N' => 2,
				'S' => 0
			},
			'ATA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'ATC' => {
				'N' => 1,
				'S' => 1
			},
			'ATG' => {
				'N' => 2,
				'S' => 0
			},
			'ATT' => {
				'N' => 1,
				'S' => 0
			},
			'CAA' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAC' => {
				'N' => 2,
				'S' => 1
			},
			'CAG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'CAT' => {
				'N' => 2,
				'S' => 0
			},
			'CCA' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCC' => {
				'N' => 2,
				'S' => 1
			},
			'CCG' => {
				'N' => 2.16666666666667,
				'S' => 0.833333333333333
			},
			'CCT' => {
				'N' => 2,
				'S' => 0
			},
			'CGA' => {
				'N' => 2,
				'S' => 1
			},
			'CGC' => {
				'N' => 2,
				'S' => 1
			},
			'CGG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'CGT' => {
				'N' => 2,
				'S' => 0
			},
			'CTA' => {
				'N' => 1,
				'S' => 1
			},
			'CTC' => {
				'N' => 1,
				'S' => 1
			},
			'CTG' => {
				'N' => 1,
				'S' => 1
			},
			'CTT' => {
				'N' => 1,
				'S' => 0
			},
			'GAA' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAC' => {
				'N' => 2,
				'S' => 1
			},
			'GAG' => {
				'N' => 2.75,
				'S' => 0.25
			},
			'GAT' => {
				'N' => 2,
				'S' => 0
			},
			'GCA' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCC' => {
				'N' => 2,
				'S' => 1
			},
			'GCG' => {
				'N' => 2.33333333333333,
				'S' => 0.666666666666667
			},
			'GCT' => {
				'N' => 2,
				'S' => 0
			},
			'GGA' => {
				'N' => 2.25,
				'S' => 0.75
			},
			'GGC' => {
				'N' => 2,
				'S' => 1
			},
			'GGG' => {
				'N' => 2.5,
				'S' => 0.5
			},
			'GGT' => {
				'N' => 2,
				'S' => 0
			},
			'GTA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTC' => {
				'N' => 1,
				'S' => 1
			},
			'GTG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'GTT' => {
				'N' => 1,
				'S' => 0
			},
			'TAA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAC' => {
				'N' => 1,
				'S' => 1
			},
			'TAG' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TAT' => {
				'N' => 1,
				'S' => 0
			},
			'TCA' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCC' => {
				'N' => 1,
				'S' => 1
			},
			'TCG' => {
				'N' => 1.5,
				'S' => 0.5
			},
			'TCT' => {
				'N' => 1,
				'S' => 0
			},
			'TGA' => { # altered because STOP
				'N' => 0, # 2
				'S' => 0 # 0
			},
			'TGC' => {
				'N' => 1,
				'S' => 1
			},
			'TGG' => {
				'N' => 2,
				'S' => 0
			},
			'TGT' => {
				'N' => 1,
				'S' => 0
			},
			'TTA' => {
				'N' => 1,
				'S' => 0
			},
			'TTC' => {
				'N' => 0,
				'S' => 1
			},
			'TTG' => {
				'N' => 1,
				'S' => 0
			}
		}
	);
	
	my $num_N_diffs = $all_diffs_hh{$codon1}->{$codon2}->{N};
	my $num_S_diffs = $all_diffs_hh{$codon1}->{$codon2}->{S};
	
	my @diffs_arr = ($num_N_diffs,$num_S_diffs);
	
	return @diffs_arr;
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
		"\n##                      SNPGenie completed successfully.                      ##".
		"\n##                Please find results in the working directory.               ##\n".
		"################################################################################".
		"\n\n\n"; 
}
