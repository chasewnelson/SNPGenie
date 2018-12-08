#! /usr/bin/perl

# PROGRAM: SNPGenie for within-group analysis is a Perl program that calculates dN/dS, 
# piN/piS, and gene diversity for a single aligned FASTA file, with accompanying CDS 
# annotations in a GTF file.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# snpgenie_within_group.pl --fasta_file_name=<aligned_seqs>.fa --gtf_file_name=<CDS_annotations>.gtf --num_bootstraps=10000 --procs_per_node=16
#########################################################################################

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

# AUTHOR: Chase W. Nelson
# Copyright (C) 2017 Chase W. Nelson

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

# SNPGenie for use with HPV!
# WAS snpgenie_between_group.pl
# We have co-opted by natural selection. That's what this is.

# ADD:
# —SORT codons before entry into the NNN-NNN comps hh, so we don't have, e.g., both AAG-AAA and AAA-AAG
# —ADD MAJORITY CODON to codon file
# —PROBLEM with min-codon-count in variant results

use strict;
#use warnings;
use Data::Dumper;
use List::Util qw(max);
#use Math::Random::OO::Bootstrap;
#use Statistics::Basic qw(:all);
use Parallel::ForkManager;
#use Statistics::Basic qw(:all); # includes median(), mean(), variance(), stddev(), covariance(), correlation()
use Getopt::Long;

STDOUT->autoflush(1);

# Get the time
my $time1 = time;
my $local_time1 = localtime;

#########################################################################################
# INITIALIZE (OPTIONAL) INPUT VARIABLES
my $fasta_file_name;
my $gtf_file_name;
my $procs_per_node;
#my $sliding_window_size;
#my $sliding_window_step;
my $num_bootstraps;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "fasta_file_name=s" => \$fasta_file_name,
			"gtf_file_name=s" => \$gtf_file_name,
			"procs_per_node=i" => \$procs_per_node,
#			"sliding_window_size:i" => \$sliding_window_size, # optional integer parameter
#			"sliding_window_step:i" => \$sliding_window_step, # optional integer parameter
			"num_bootstraps=i" => \$num_bootstraps)
			
			or die "\n### WARNING: Error in command line arguments. Script terminated.\n\n";
			# If an argument is called as a flag, its value is 0; if not called, it's null

unless($fasta_file_name =~ /.fa/) {
	die "\n### WARNING: The --fasta_file_name option must be a file with a .fa or .fasta extension\n".
		"### SNPGenie terminated.\n\n";
}

unless($gtf_file_name =~ /.gtf/) {
	die "\n### WARNING: The --gtf_file_name option must be a file with a .gtf or .gtf extension\n".
		"### SNPGenie terminated.\n\n";
}

if(! $procs_per_node) {
	if($procs_per_node != 0) {
		$procs_per_node = 1; # DEFAULT is 1
	}
#	$param_file_contents .= "MINIMUM ALLELE FREQUENCY: Default used; all SNPs included\n";
} elsif($procs_per_node < 1) {
	die "\n### WARNING: The --procs_per_node option must be an integer ≥1\n".
		"### SNPGenie terminated.\n\n";
} #else {
#	$param_file_contents .= "MINIMUM ALLELE FREQUENCY: $minfreq\n";
#}

#if(! $sliding_window_size) {
#	if($sliding_window_size != 0) { # Called as a flag, but given no value
#		$sliding_window_size = 10; # default behavior: 10-mer peptide
#	}
##	$param_file_contents .= "SLIDING WINDOW LENGTH: Default used; 9 codons\n";
#} elsif($sliding_window_size < 1) {
#	die "\n### WARNING: The --sliding_window_size option must be an integer ≥1\n".
#		"### SNPGenie terminated.\n\n";
#} #else {
##	$param_file_contents .= "SLIDING WINDOW LENGTH: None\n";
##}
#
#if(! $sliding_window_step) {
#	if($sliding_window_step != 0) { # Called as a flag, but given no value
#		$sliding_window_step = 1; # default behavior: one codon
#	}
#} elsif($sliding_window_step < 1) {
#	die "\n### WARNING: The --sliding_window_step option must be an integer ≥1\n".
#		"### SNPGenie terminated.\n\n";
#}

if(! $num_bootstraps) {
	if($num_bootstraps != 0) { # Called as a flag, but given no value
		$num_bootstraps = 1000; # default behavior: nonamer peptide
	}
#	$param_file_contents .= "SLIDING WINDOW LENGTH: Default used; 9 codons\n";
} elsif($num_bootstraps < 2) {
	die "\n### WARNING: The --num_bootstraps option must be an integer ≥2\n".
		"### SNPGenie terminated.\n\n";
} #else {
#	$param_file_contents .= "SLIDING WINDOW LENGTH: None\n";
#}


### Prepare for bootstrapping
##if(! $num_bootstraps) { # null or 0
##	$num_bootstraps = 0; # default behavior
##} elsif(! $num_bootstraps >= 10) {
##	die "\n### WARNING: The --num_bootstraps option must ≥10\n".
##		"### Script terminated.\n\n";
##}
##
### Optional BOOTSTRAP argument
##if($num_bootstraps >= 10) {
##	print "\n# Bootstrapping will be performed with $num_bootstraps replicates.\n";
##} else {
##	print "\n# No bootstrapping will occur.\n";
##}
##
### Set up parallelism
##if(! $procs_per_node) { # null or 0
##	$procs_per_node = 1; # default behavior
##} elsif($procs_per_node < 1) {
##	die "\n### WARNING: The --procs_per_node option must ≥1\n".
##		"### Script terminated.\n\n";
##}

print "\n# Number of parallel processes to be invoked is $procs_per_node\.\n";

print "\n################################################################################".
	"\n##                                                                            ##".
	"\n##                       Within-Group SNPGenie Initiated!                     ##".
	"\n##                                                                            ##".
	"\n################################################################################\n";


print "\nSNPGenie initiated at local time $local_time1\n";

#my $fasta_file_name = $ARGV[0];
#my $gtf_file_name = $ARGV[1];
#unless($fasta_file_name =~ /.fa/) { die "\n\n# FASTA file must contain .fa or .fasta extension. TERMINATED\n\n"; }
#unless($gtf_file_name =~ /.gtf/) { die "\n\n# GTF file must contain .gtf extension. TERMINATED\n\n"; }

# Store product information
my @product_names_arr = &get_product_names_from_gtf($gtf_file_name);
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

# Read in the group of sequences from the fasta file
my $seq = '';
my @seqs_arr;
#my $header = '';
#my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file_name") or die "Could not open file $fasta_file_name\n";

print "\nRecording coding sequence data for $fasta_file_name...\n";

while(<IN_FASTA>) {
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
			#print "\nseq: $seq\n";
			
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

close IN_FASTA;

$seq = uc($seq);
$seq =~ tr/U/T/;
push(@seqs_arr,$seq);
#push(@headers_arr,$header);

print "\n";

# Go through each product and DO DIS THANG
# HERE'S HOW to get the first array:
#print "@{$seqs_aa[0]}\n";

#my $counter=1;
#foreach(@{$seqs_aa[0]}) {
#	print "$counter: $_\n\n";
#	$counter++;
#}

my %products_seqs_hh;
# FORMAT:
#%products_seqs_hh = {
#	'ORF1a' => {
#		'seq_1' => 'ACGT',
#		'seq_2' => 'ACGT',
#		... ,
#		'seq_n' => 'ACGT',
#	},
#	
#	'ORF1b' => {
#		'seq_1' => 'ACGT',
#		'seq_2' => 'ACGT',
#		... ,
#		'seq_n' => 'ACGT',
#		}
#	}
#};

foreach(@product_names_arr) { # FOR EACH PRODUCT
	my $product_name = $_;
	my @product_coord_arr = @{$product_coordinates_ha{$product_name}->{product_coord_arr}};
	
	# New segments approach
	my %product_starts;
	my %product_stops;
	
	my $num_segments = (@product_coord_arr / 2);
	for(my $i=1; $i<=$num_segments; $i++) { # $i<=scalar(@product_coord_arr) # FOR EACH SEGMENT (2) # COME BACK COMEBACK
		$product_starts{$i} = $product_coord_arr[2*$i-2];
		$product_stops{$i} = $product_coord_arr[2*$i-1];
	}
	
	# Build an array of all the sequences, JUST ONE GROUP
	my $seq_num = 0;
	foreach my $sequence (@seqs_arr) { # for each seq in alignment (max coverage)
	
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
		$products_seqs_hh{$product_name}->{$seq_id} = $this_product_seq;
	}
}

#my @test_keys = keys %{$products_groups_seqs_hh{'ORF5'}->{'group_1'}};
#@test_keys = keys %products_groups_seqs_hh;
#print "\n\n @test_keys \n\n";
#my $test_seq = $products_groups_seqs_hh{'ORF5'}->{'group_1'}->{'seq_1'};
#print "\n\n $test_seq \n\n";

open(CODON_FILE,">>within\_group\_codon\_results\.txt");
print CODON_FILE "file\tproduct\tcodon\tvariability\tcomparisons\t".
	"N_sites\tS_sites\tN_diffs\tS_diffs\n";
close CODON_FILE;

open(VARIANT_FILE,">>within\_group\_variant\_results\.txt");
print VARIANT_FILE "file\tproduct_name\tcodon_num\tconsensus_codon\tnum_alleles\t".
	"codons\tamino_acids\tvariant_type\tcodon_counts\tcount_total\tcodon_freqs\t".
	"max_codon_count\tmin_codon_count\tnum_defined_seqs\n";
close VARIANT_FILE;

open(OUTFILE,">>within\_group\_product\_results\.txt");
print OUTFILE "file\tproduct\tN_sites\tS_sites\tN_diffs\tS_diffs\tdN\tdS\t";

if($num_bootstraps > 1) {
	print OUTFILE "SE_dN\tSE_dS\t";
}

print OUTFILE "dN_minus_dS\tdN_over_dS";
	
if($num_bootstraps > 1) {
	print OUTFILE "\tSE_dN_minus_dS\tZ_value\tsignificance\n";
} else {
	print OUTFILE "\n";
}

## Preparing sliding window output
#open(OUTFILE_WITHIN_SW,">>within\_group\_sw_" . $sliding_window_size . "codons\_results.txt");
#print OUTFILE_WITHIN_SW "analysis\tfamily\tpartition\tgroup_1\tgroup_2\twindow\tfirst_codon\tlast_codon\t".
#	"num_defined_codons_g1\tnum_defined_codons_g2\tmean_num_defined_codons\tmin_num_defined_codons\t".
#	"N_sites\tS_sites\t".
#	"N_diffs\tS_diffs\t".
#	"dN\tdS\t".
#	"dN\-dS\t".
#	"dN\/dS\tdN\>total_dS";
#	
#if($num_bootstraps > 1) {
#	print OUTFILE_WITHIN_SW "\tSE(dN-dS)\tZ_value\tsignificance\n";
#} else {
#	print OUTFILE_WITHIN_SW "\n";
#}
#
#close OUTFILE_WITHIN_SW;



#print "product\tN_sites\tS_sites\tN_diffs\tS_diffs\n";

# YES, it all works
#my %product_data_hh;
my $total_N_sites = 0;
my $total_S_sites = 0;
my $total_N_diffs = 0;
my $total_S_diffs = 0;

# We HAVE all the sequences stored by product in %products_seqs_hh as 
# $products_seqs_hh{$product_name}->{$seq_id}
# Next, let's simply build/analyze all codons for each as we go; so
# (1) loop products
# (2) loop codons
# (3) analyze as we go

#my $last_num_codons;

my %min_count2num_alleles_TOTALS;

foreach my $product_name (keys %products_seqs_hh) { # for each product

	# FOREACH SEQUENCE IN FASTA?
	# foreach my $sequence (@seqs_arr)

	print "\n################################################################################".
		"\nANALYZING $product_name:\n";
		
	# I want, for this product:
	#codonum_codon_aa[codon_index]->array of all codons at that site
	my @codonum_codon_aa;	
	my $last_num_codons;
	
	# STORE all codons here
	foreach my $seq_id (keys %{$products_seqs_hh{$product_name}}) { # for each sequence	
		
		my $product_seq = $products_seqs_hh{$product_name}->{$seq_id};
		my $product_length = length($product_seq);
		
		#print "\nProduct $product_name is $product_length nucleotides long:\n$product_seq\n\n";
		
		#foreach codon as in 1 + 3 + 3
		# Make sure it's a complete set of codons
		if(($product_length % 3) != 0) {
			die "\n\nDIE: A sequence in $product_name is not a multiple of 3 (complete codon set). TERMINATED.\n\n";
		}
		
		# Build all codons, put in @aa with @codon_num->@codons
		my $num_codons = $product_length / 3;
		
		if((! $last_num_codons) && ($last_num_codons ne '') && ($num_codons != $last_num_codons)) {
			die "\n\nDIE: In $product_name, there are sequences of different length ($last_num_codons and $num_codons). TERMINATED.\n\n";
		} else {
			$last_num_codons = $num_codons;
			#print "changed\n";
		}
		
		my $codon_index = 0; # so we're going to populate $codonum_codon_aa[0]->@
		for(my $i=1; $i<=$num_codons; $i++) { # for each codon
			my $array_index = $i - 1;
			#my $codon = substr($sequence,$codon_index,3);
			# More efficient to call substr at beginning and gobble
			my $codon = substr($product_seq,0,3,"");
			#print "$codon ";
			
			# Note that some may contain N's
			push(@{$codonum_codon_aa[$array_index]},$codon);
			
			#$codon_index+=3;
		}
		#print "\n";
	} # end sequences loop; end STORING codon information
	
	# ANALYZE all codons here; compare all codons to one another; all pairwise comparisons
	# We're going one product at a time, so we're within just one product now, 
	# having all codons stored in @codonum_codon_aa{INDEX IN SEQ}->@ARRAY OF CODONS AT INDEX FOR ALL SEQS
	# e.g., SIX (6) BIRDS, 15 comparisons


	my $num_seqs = scalar @{$codonum_codon_aa[0]};
	my @num_seqs_defined;
	#print "\nNum seqs: $num_seqs\n";
	my $num_codons_in_product = scalar @codonum_codon_aa;
	my @polymorphic_codons_arr; # same indices as codons; value 0 if invariant, 1 if poly

	
	# PARALLELIZE IT
	mkdir("$product_name\_polymorphic_codons_arr");
	chdir("$product_name\_polymorphic_codons_arr");
	my $pm_poly = Parallel::ForkManager->new($procs_per_node);
	
	# QUICKLY CHECK IF IT'S POLYMORPHIC to save time. It's faster just to check,
	# EVEN IF not polymorphic
	# FOR EACH CODON IN PRODUCT
	for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) {
		
		$pm_poly->start and next; # this is IT
		
		#print "Codon $codon_index\n";
		#my $codon_to_print = $groups_codons_aa[$i]->[$codon_index]->[1];
		#$| = 1;
		#print "$codon_to_print";
		my $between_s_comparisons = 0;
		my $this_codon_poly = 0;
		my $this_codon_num_seqs_defined = $num_seqs;
		
		OUTER: for (my $si_seq_index = 0; $si_seq_index < $num_seqs; $si_seq_index++) { # for each si sequence at this codon
			my $codon_si = $codonum_codon_aa[$codon_index]->[$si_seq_index];
			# Codon to compare against
			
			if(! ($codon_si =~ 'N') && ! ($codon_si =~ '-') ) {
				INNER: for (my $sj_seq_index = $si_seq_index+1; $sj_seq_index < $num_seqs; $sj_seq_index++) { # for each sj sequence at this codon
					#if ($codon_index == 0) { print "pw comparison\n"; }
					my $codon_sj = $codonum_codon_aa[$codon_index]->[$sj_seq_index];
					#print "Comparing codons $codon_si and $codon_sj\n";
					if($codon_si ne $codon_sj && ! ($codon_sj =~ 'N') && ! ($codon_sj =~ '-')) {
						$this_codon_poly = 1;
						last OUTER; # gotta break of out TWO loops here and go to next codon in product
					}
				} # INNER
			} else {
				$this_codon_num_seqs_defined -= 1;
			}
		} # OUTER
		
		open(THIS_POLY_TEMP_FILE,">>$codon_index");
		print THIS_POLY_TEMP_FILE "$this_codon_poly\t$this_codon_num_seqs_defined";
		close THIS_POLY_TEMP_FILE;
				
		$pm_poly->finish; # special name
		
#		push(@polymorphic_codons_arr,$this_codon_poly);
#		$num_seqs_defined[$codon_index] = $this_codon_num_seqs_defined;
	} # done checking if polymorphic
	
	
	$pm_poly->wait_all_children; # special name, methods within module
			
	my @polymorphic_codons_FILES_arr = glob "*";
	@polymorphic_codons_FILES_arr = sort {$a <=> $b} @polymorphic_codons_FILES_arr;
	
	#print "sorted polymorphic_codons_FILE_arr: @polymorphic_codons_FILES_arr\n";
	
	foreach(@polymorphic_codons_FILES_arr) { # file names, actually
		my $poly_file = $_;
		
		open(CURR_POLY_FILE, $poly_file) or die "\n## Cannot open $poly_file. TERMINATED.\n\n";
		while(<CURR_POLY_FILE>) {
			chomp;
			my @line_arr = split(/\t/,$_,-1);
			push(@polymorphic_codons_arr,$line_arr[0]);	
			$num_seqs_defined[$poly_file] = $line_arr[1];
		}
		close CURR_POLY_FILE;
		unlink $poly_file;
	}
	
	chdir("..");
	rmdir("$product_name\_polymorphic_codons_arr");
	# Finish parallelism for polymorphism detection
	
	
	
	
	
	
	
	#print "\nDone recording polymorphism in $product_name\nIt is: @polymorphic_codons_arr\n\n";
			
	print "\nAnalyzing polymorphic codons in $product_name\...\n";
		
	# Calculate number of total codons so that we can track % progress?
	# INTRO: skip if not poly, just do values for the first codon observed
	
	my $product_N_sites_sum = 0;
	my $product_S_sites_sum = 0;
	my $product_N_diffs_sum = 0;
	my $product_S_diffs_sum = 0;
	
	open(CODON_FILE,">>within\_group\_codon\_results\.txt");
	
	# COMEBACK: make a number of sites subroutine for the whole codon, not just the sites
	# OR, instead of calling the subroutine every time, just COUNT the comparisons!
	# $hh_comps{CGA}->{CGG}->256 YES!!	
	





	# PARALLEL
##	mkdir("$product_name\_codon\_analysis");
##	chdir("$product_name\_codon\_analysis");
##	my $pm_codons = Parallel::ForkManager->new($procs_per_node);
	# To do this, must write all output to files during the loop
	# after all codons finished, then we must read through and scoop results
	# while storing in the appropriate places


	### NEW POLYMORPHIC SITE ANALYSIS
	my %poly_site_data_hh; # $poly_site_data_hh{codon_index}->{unique_variants}
	my %codon_to_results_hh; # for bootstraps
	
	COMPL_DEL: for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) { # for each codon in product
##		$pm_codons->start and next; # this is IT
		
		my $codon_num = $codon_index+1;
		
		print CODON_FILE "$fasta_file_name\t$product_name\t" . $codon_num . " \t";
			
		if($polymorphic_codons_arr[$codon_index]) { # if codon is polymorphic, value is 1
			print CODON_FILE "polymorphic\t";
			
			my %comps_hh;
			for (my $si_seq_index = 0; $si_seq_index < $num_seqs; $si_seq_index++) { # for each si sequence at this codon
				my $codon_si = $codonum_codon_aa[$codon_index]->[$si_seq_index];
				
				# Sequence (codon) to compare against
				if(! ($codon_si =~ 'N') && ! ($codon_si =~ '-')) {
					for (my $sj_seq_index = $si_seq_index+1; $sj_seq_index < $num_seqs; $sj_seq_index++) { # for each sj sequence at this codon
						my $codon_sj = $codonum_codon_aa[$codon_index]->[$sj_seq_index];
						if(! ($codon_sj =~ 'N') && ! ($codon_sj =~ '-')) { # they don't have to be the same; syn still stored here
							$comps_hh{$codon_si}->{$codon_sj}+=1;
						}
					}
				}
			}
			
			# SUM UP STUFF HERE
			my $between_s_comparisons = 0;
			my $sum_N_sites = 0;
			my $sum_S_sites = 0;
			my $sum_N_diffs = 0;
			my $sum_S_diffs = 0;
			
			foreach my $codon_si (keys %comps_hh) {
				foreach my $codon_sj (keys %{$comps_hh{$codon_si}}) {
					
					if($codon_si =~ /-/ || $codon_sj =~ /-/) { # these will skip all additions
						print CODON_FILE "COMPLETE_DELETION\tNA\tNA\tNA\tNA\n";
						next COMPL_DEL;
					} else {						
						my $weight = $comps_hh{$codon_si}->{$codon_sj};
						
						### NEW POLYMORPHIC SITE ANALYSIS ADDED
						$poly_site_data_hh{$codon_index}->{$codon_si} += $weight;
						$poly_site_data_hh{$codon_index}->{$codon_sj} += $weight;
						### PARALLEL: NEED TO FIGURE THIS OUT
						
						print CODON_FILE "($codon_si\-$codon_sj)";
##						$this_codon_output .= "($codon_si\-$codon_sj)";
						
						# Sites codon gi
						my @codon_si_sites_1_arr = &get_number_of_sites($codon_si,1);
						my @codon_si_sites_2_arr = &get_number_of_sites($codon_si,2);
						my @codon_si_sites_3_arr = &get_number_of_sites($codon_si,3);
					
						my $codon_si_N_sites = ($codon_si_sites_1_arr[0] + $codon_si_sites_2_arr[0] + $codon_si_sites_3_arr[0]);
						my $codon_si_S_sites = ($codon_si_sites_1_arr[1] + $codon_si_sites_2_arr[1] + $codon_si_sites_3_arr[1]);							
						
						# Site codon gj
						my @codon_sj_sites_1_arr = &get_number_of_sites($codon_sj,1);
						my @codon_sj_sites_2_arr = &get_number_of_sites($codon_sj,2);
						my @codon_sj_sites_3_arr = &get_number_of_sites($codon_sj,3);
					
						my $codon_sj_N_sites = ($codon_sj_sites_1_arr[0] + $codon_sj_sites_2_arr[0] + $codon_sj_sites_3_arr[0]);
						my $codon_sj_S_sites = ($codon_sj_sites_1_arr[1] + $codon_sj_sites_2_arr[1] + $codon_sj_sites_3_arr[1]);
						
						#print "\nComparing codons $codon_si and $codon_sj:\n";
						
						my $mean_comp_N_sites = ($codon_si_N_sites + $codon_sj_N_sites) / 2;
						my $mean_comp_S_sites = ($codon_si_S_sites + $codon_sj_S_sites) / 2;
						
						$sum_N_sites += $weight * $mean_comp_N_sites;
						$sum_S_sites += $weight * $mean_comp_S_sites;
						
						# Differences
						my $N_diffs = 0;
						my $S_diffs = 0;
						
						if($codon_si ne $codon_sj) {
							my @diffs_arr = &return_avg_diffs($codon_si,$codon_sj);
							$N_diffs = $diffs_arr[0];
							$S_diffs = $diffs_arr[1];
						}
						
						$sum_N_diffs += $weight * $N_diffs;
						$sum_S_diffs += $weight * $S_diffs;							
	
						$between_s_comparisons += $weight;
					}
				}
			}
			print CODON_FILE "\t";
##			$this_codon_output .= "\t";
			
			#print "\nfor a total of $between_s_comparisons comparisons\n";
			
			my $mean_N_diffs = ($sum_N_diffs / $between_s_comparisons);
			my $mean_S_diffs = ($sum_S_diffs / $between_s_comparisons);
			my $mean_N_sites = ($sum_N_sites / $between_s_comparisons);
			my $mean_S_sites = ($sum_S_sites / $between_s_comparisons);
			
			
			print CODON_FILE "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";
##			$this_codon_output .= "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";
			
			$product_N_diffs_sum += $mean_N_diffs;
			$product_S_diffs_sum += $mean_S_diffs;
			$product_N_sites_sum += $mean_N_sites;
			$product_S_sites_sum += $mean_S_sites;
			
			
			#print ".. and it's all done.\n";
			
			# For bootstraps
			if($num_bootstraps > 1) {
				$codon_to_results_hh{$codon_num}->{N_diffs} = $mean_N_diffs;
				$codon_to_results_hh{$codon_num}->{S_diffs} = $mean_S_diffs;
				$codon_to_results_hh{$codon_num}->{N_sites} = $mean_N_sites;
				$codon_to_results_hh{$codon_num}->{S_sites} = $mean_S_sites;
			}
			
		} else { # not polymorphic; just use first codon from the first group because it's conserved
			my $conserved_codon = '';
			
			FIND_CONSERVED_CODON: for(my $k = 0; $k < $num_seqs; $k++) { # for each codon in product
				my $curr_codon_rep = $codonum_codon_aa[$codon_index]->[$k]; # grab first codon THAT DOESN'T CONTAIN N
				if($curr_codon_rep =~ "N" || $curr_codon_rep =~ "-") {
					next FIND_CONSERVED_CODON;
				} else {
					$conserved_codon = $curr_codon_rep;
					last FIND_CONSERVED_CODON;
				}
			}
			
			my @codon_sites_1_arr = &get_number_of_sites($conserved_codon,1);
			my @codon_sites_2_arr = &get_number_of_sites($conserved_codon,2);
			my @codon_sites_3_arr = &get_number_of_sites($conserved_codon,3);
			
			my $codon_N_sites = ($codon_sites_1_arr[0] + $codon_sites_2_arr[0] + $codon_sites_3_arr[0]);
			my $codon_S_sites = ($codon_sites_1_arr[1] + $codon_sites_2_arr[1] + $codon_sites_3_arr[1]);
			
			print CODON_FILE "conserved\t$conserved_codon\t".
				"$codon_N_sites\t$codon_S_sites\t0\t0\n";
			
##			$this_codon_output .= "conserved\t$conserved_codon\t".
##				"$codon_N_sites\t$codon_S_sites\t0\t0\n";
			
			$product_N_sites_sum += $codon_N_sites;
			$product_S_sites_sum += $codon_S_sites;
			
			
			# For bootstraps
			if($num_bootstraps > 1) {
				$codon_to_results_hh{$codon_num}->{N_diffs} = 0;
				$codon_to_results_hh{$codon_num}->{S_diffs} = 0;
				$codon_to_results_hh{$codon_num}->{N_sites} = $codon_N_sites;
				$codon_to_results_hh{$codon_num}->{S_sites} = $codon_S_sites;
			}
			
		} # end not polymorphic
		
##		open(THIS_CODON_TEMP_FILE,">>$codon_index");
##		print THIS_CODON_TEMP_FILE "$this_codon_output";
##		close THIS_CODON_TEMP_FILE;
				
##		$pm_codons->finish; # special name
	} # end all codons in product
	
	close CODON_FILE;

	
	
	
##	# PARALLEL
##	$pm_codons->wait_all_children; # special name, methods within module
##
##	# Vacuum up output and delete temp files	
##	my @codon_lines_arr;	
##	my @codon_analysis_FILES_arr = glob "*";
##	@codon_analysis_FILES_arr = sort {$a <=> $b} @codon_analysis_FILES_arr;
##	
##	#print "sorted polymorphic_codons_FILE_arr: @polymorphic_codons_FILES_arr\n";
##	
##	foreach(@codon_analysis_FILES_arr) { # file names, actually
##		my $codon_file = $_; # this is the INDEX
##		
##		open(CURR_CODON_FILE, $codon_file) or die "\n## Cannot open $codon_file. TERMINATED.\n\n";
##		while(<CURR_CODON_FILE>) {
##			chomp;
##			
##			push(@codon_lines_arr,$_);
##			
##			my @this_line_arr = split(/\t/,$_,-1);
##			#product / group_1 / group_2 / codon / variability / comparisons / N_sites / S_sites / N_diffs / S_diffs
##			my $poly = $this_line_arr[4];
##			
##			# IF POLY
##			if($poly eq 'polymorphic') {
##				$product_N_sites_sum += $this_line_arr[6];
##				$product_S_sites_sum += $this_line_arr[7];
##				$product_N_diffs_sum += $this_line_arr[8];
##				$product_S_diffs_sum += $this_line_arr[9];
##			
##			# IF CONSERVED
##			} elsif($poly eq 'conserved') {
##				$product_N_sites_sum += $this_line_arr[6];
##				$product_S_sites_sum += $this_line_arr[7];
##				
##				#STORE THIS CONSERVED CODON FOR LATER USE IN BOOTSTRAPPING
##				$codonIndex_conserved{$codon_file} = $this_line_arr[5]; # COMPARISONS column
##			# PROBLEM
##			} else {
##				die "\nPROBLEM: site neither conserved nor polymorphic!? DIE\n\n";
##			}
##			
##			
##			
##
##
##		}
##		close CURR_CODON_FILE;
##		unlink $codon_file;
##	}
##	
##	chdir("..");
##	rmdir("$product_name\_codon\_analysis");
##	# done sweeping output
##
##	# Concatenate and print parallelized codon results to a results file
##	open(CODON_FILE,">>between\_group\_codon\_results\.txt");
##	foreach(@codon_lines_arr) {
##		print CODON_FILE "$_\n";
##	}
##	close CODON_FILE;
	
	


	#################
	# BOOTSTRAPPING
	my $SE_piN_minus_piS = 0;
	my $product_boot_Z = 'NA';
	
	my $SE_dN;
	my $SE_dS;
	
	if($num_bootstraps > 1) {
		# MAKE BOOTSTRAP FILE
		open(BOOT_FILE_PRODUCT,">>$product_name\_bootstrap\_results\.txt");
		print BOOT_FILE_PRODUCT "file\tproduct_name\tbootstrap_num\tsim_product_N_sites_sum\t".
				"sim_product_S_sites_sum\tsim_product_N_diffs_sum\tsim_product_S_diffs_sum\n";
		
		#################
		# BOOTSTRAPPING TO CALCULATE STANDARD ERROR HERE?
		print "\nBootstrapping for standard error...\n";
		
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
			for(my $codon_index = 0; $codon_index < $num_codons_in_product; $codon_index++) {
				# Which codon site do we choose?
				my $random_codon_num = int(rand($num_codons_in_product + 1));
				
				my $mean_N_diffs = $codon_to_results_hh{$random_codon_num}->{N_diffs};
				my $mean_S_diffs = $codon_to_results_hh{$random_codon_num}->{S_diffs};
				my $mean_N_sites = $codon_to_results_hh{$random_codon_num}->{N_sites};
				my $mean_S_sites = $codon_to_results_hh{$random_codon_num}->{S_sites};
				
				$sim_product_N_diffs_sum += $mean_N_diffs;
				$sim_product_S_diffs_sum += $mean_S_diffs;
				$sim_product_N_sites_sum += $mean_N_sites;
				$sim_product_S_sites_sum += $mean_S_sites;
				
#				print CODON_FILE "$mean_N_sites\t$mean_S_sites\t$mean_N_diffs\t$mean_S_diffs\n";

			} # finished compiling all sampled codons (cols in alignment)
			
			# Print product totals for BOOTSTRAP
			my $out_line_boot = "$fasta_file_name\t$product_name\t$bootstrap_num\t$sim_product_N_sites_sum\t".
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
		
		my $actual_piN = '*'; 
		my $actual_piS = '*';
		my $actual_piN_minus_piS = '*';
		
		if($product_N_sites_sum > 0) {
			$actual_piN = $product_N_diffs_sum / $product_N_sites_sum;
		}
		
		if($product_S_sites_sum > 0) {
			$actual_piS = $product_S_diffs_sum / $product_S_sites_sum;
		}
		
		if($actual_piN >= 0 && $actual_piS >= 0) {
			$actual_piN_minus_piS = $actual_piN - $actual_piS;
		}
		

		# CALCULATE BOOTSTRAP STANDARD ERROR HERE! NEI & KUMAR (2000) EQUATIONS
		my @sim_piN_minus_piS;
		my @sim_dN;
		my @sim_dS;
		
		for(my $sim_num = 0; $sim_num < scalar(@sim_N_sites_arr); $sim_num++) {
			my $this_round_piN = '*';
			if($sim_N_sites_arr[$sim_num] > 0) {
				$this_round_piN = $sim_N_diffs_arr[$sim_num] / $sim_N_sites_arr[$sim_num];
			}
			
			my $this_round_piS = '*';
			if($sim_S_sites_arr[$sim_num] > 0) {
				$this_round_piS = $sim_S_diffs_arr[$sim_num] / $sim_S_sites_arr[$sim_num];
			}
			
			my $this_round_piN_minus_piS = '*';
			if($this_round_piN >= 0 && $this_round_piS >= 0) {
				$this_round_piN_minus_piS = $this_round_piN - $this_round_piS;
			}
			
			push(@sim_piN_minus_piS,$this_round_piN_minus_piS);
			push(@sim_dN, $this_round_piN);
			push(@sim_dS, $this_round_piS);
		}
		
		$SE_piN_minus_piS = &standard_deviation(@sim_piN_minus_piS);
		$SE_dN = &standard_deviation(@sim_dN);
		$SE_dS = &standard_deviation(@sim_dS);
		
		if($SE_piN_minus_piS > 0) {
			$product_boot_Z = $actual_piN_minus_piS / $SE_piN_minus_piS;
		}
		
		if($SE_piN_minus_piS == 0) {
			$SE_piN_minus_piS = 'NA';
		}
		
#		print "\nproduct SE(dN-dS) = $SE_piN_minus_piS\n".
#			"product Z-value = $product_boot_Z\n\n";
		
	} # END BOOTSTRAPS
	
	
	# Print product totals
	my $product_piN;
	my $product_piS;
	my $product_piN_over_piS;
	
	if($product_N_sites_sum > 0) {
		$product_piN = $product_N_diffs_sum / $product_N_sites_sum;
	} else {
		$product_piN = '*';
	}
	
	if($product_S_sites_sum > 0) {
		$product_piS = $product_S_diffs_sum / $product_S_sites_sum;
	} else {
		$product_piS = '*';
	}
	
	my $product_piN_minus_piS = $product_piN - $product_piS;
	
	if($product_piS > 0) {
		$product_piN_over_piS = $product_piN / $product_piS;
	} else {
		$product_piN_over_piS = '*';
	}
	
	my $out_line = "$fasta_file_name\t$product_name\t$product_N_sites_sum\t".
		"$product_S_sites_sum\t$product_N_diffs_sum\t$product_S_diffs_sum\t".
		"$product_piN\t$product_piS\t";
		
	if($num_bootstraps > 1) {
		$out_line .= "$SE_dN\t$SE_dS\t";
	}
		
	$out_line .= "$product_piN_minus_piS\t$product_piN_over_piS";
	
	my $rounded_dN = sprintf("%.5f",$product_piN);
	my $rounded_dS = sprintf("%.5f",$product_piS);
	my $rounded_dNdS = sprintf("%.5f",$product_piN_over_piS);
	
	if($product_piN_over_piS eq '*') {
		$rounded_dNdS = '*';
	}
	
	print "\ndN=$rounded_dN\ndS=$rounded_dS\ndN\/dS=$rounded_dNdS\n";
	
	if($num_bootstraps > 1) {
#		my $z_value;
		my $rounded_SE = sprintf("%.5f",$SE_piN_minus_piS);
		my $rounded_Z = sprintf("%.5f",$product_boot_Z);
		
		print "\nSE(dN-dS) = $rounded_SE\n".
			"Z-value = $rounded_Z\n\n";
		
		my $significance = '';
		
		if($SE_piN_minus_piS > 0) {
		
#			$z_value = $product_piN_minus_piS / $SE_piN_minus_piS; 
			
			if(abs($product_boot_Z) > 2.81) {
				$significance = '***';
			} elsif(abs($product_boot_Z) > 1.96) {
				$significance = '**';
			} elsif(abs($product_boot_Z) > 1.64) {
				$significance = '*';
			}
		
		}
		
		$out_line .= "\t$SE_piN_minus_piS\t$product_boot_Z\t$significance\n"; #####
		
		
	} else {
		$out_line .= "\n";
	}
	
	print OUTFILE "$out_line";
#	print "ACTUAL:\n$out_line\n\n";
	
	
	open(VARIANT_FILE,">>within\_group\_variant\_results\.txt");
	
	### NEW POLYMORPHIC SITE ANALYSIS ADDED
	#$poly_site_data_hh{$codon_index}->{$codon_si}; # it's a read depth count
	# LOOP POLY CODONS; CODON NUMBER WILL BE INDEX + 1
	
	my %min_count2num_alleles;
	
	foreach my $codon_index (sort {$a <=> $b} keys %poly_site_data_hh) {
		my $codon_num = $codon_index + 1;
		
		my $num_alleles = scalar(keys %{$poly_site_data_hh{$codon_index}});
		
		my $codons = '';
#		my $num_alleles = 0;
		my $codon_counts = '';
		my @codon_counts;
		my $count_total = 0;
		
		my $amino_acids = '';
		my @amino_acids_arr;
		
		my $max_codon_count = 0;
		my $consensus_codon = '';
		my $ambiguous_flag = 0;
		
		foreach my $codon_triplet (sort keys %{$poly_site_data_hh{$codon_index}}) {
			$codons .= "$codon_triplet\,";
#			$num_alleles ++;
			my $curr_codon_count = ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs_defined[$codon_index]-1));
			$codon_counts .= $curr_codon_count . "\,";
			push(@codon_counts,($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs_defined[$codon_index]-1)));
			$count_total += ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs_defined[$codon_index]-1));
			
			my $this_amino_acid = &get_amino_acid($codon_triplet);
			$amino_acids .= "$this_amino_acid\,";
			push(@amino_acids_arr,$this_amino_acid);
			
#			$codons .= "$codon_triplet\,";
#			$num_alleles ++;
#			my $curr_codon_count = ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs-1));
#			$codon_counts .= $curr_codon_count . "\,";
#			push(@codon_counts,($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs-1)));
#			$count_total += ($poly_site_data_hh{$codon_index}->{$codon_triplet} / ($num_seqs-1));
			
			# CONSENSUS CODON
			if($curr_codon_count > $max_codon_count) {
				$max_codon_count = $curr_codon_count;
				$consensus_codon = $codon_triplet;
			} elsif($curr_codon_count == $max_codon_count) {
				$ambiguous_flag = 1;
			}
		}
		
		my $nonsyn_comps_count = 0;
		my $syn_comps_count = 0;
		# Determine if variants nonsynonymous, synonymous, or ambiguous at this codon
		
		for(my $i = 0; $i < @amino_acids_arr; $i++) {
			my $amino_acid_1 = $amino_acids_arr[$i];
			
			for(my $j = $i+1; $j < @amino_acids_arr; $j++) {
				my $amino_acid_2 = $amino_acids_arr[$j];
				
				if($amino_acid_1 eq $amino_acid_2) {
					$syn_comps_count++;
				} else {
					$nonsyn_comps_count++;
				}
				
			}
		}
		
		my $variant_type;
		if($nonsyn_comps_count > 0 && $syn_comps_count > 0) {
			$variant_type = "Ambiguous";
		} elsif($nonsyn_comps_count > 0) {
			$variant_type = "Nonsynonymous";
		} else {
			$variant_type = "Synonymous";
		}
		
		#my @freqs_arr;
		my $codon_freqs = '';
		
		foreach (@codon_counts) {
			my $freq = ($_ / $count_total);
			#push(@freqs_arr, $freq);
			my $freq_rounded = sprintf("%.3f",$freq);
			$codon_freqs .= "$freq_rounded\,";
		}
		
		# GOBBLE / TRIM ENDS OF EACH COMPONENT
		chop($codons);
		chop($codon_counts);
		chop($codon_freqs);
		chop($amino_acids);
		
		# NUM SEQS MINOR ALLELE
		my $min_codon_count = $num_seqs_defined[$codon_index] - $max_codon_count;
#		my $min_codon_count = $num_seqs - $max_codon_count;
		
		my $out_variant_line = "$fasta_file_name\t$product_name\t$codon_num\t$consensus_codon\t$num_alleles\t".
			"$codons\t$amino_acids\t$variant_type\t".
			"$codon_counts\t$count_total\t$codon_freqs\t$max_codon_count\t$min_codon_count\t" . $num_seqs_defined[$codon_index];
		
		print VARIANT_FILE "$out_variant_line\n";
		
		$min_count2num_alleles{$min_codon_count}->{$num_alleles}++;
		$min_count2num_alleles_TOTALS{$min_codon_count}->{$num_alleles}++;
		
	} # end loop poly codons within this product
	
	close VARIANT_FILE;
	
	$total_N_sites += $product_N_sites_sum;
	$total_S_sites += $product_S_sites_sum;
	$total_N_diffs += $product_N_diffs_sum;
	$total_S_diffs += $product_S_diffs_sum;
	
	foreach my $min_codon_count (sort {$a <=> $b} keys %min_count2num_alleles) {
		
		foreach my $num_alleles (sort {$a <=> $b} keys %{$min_count2num_alleles{$min_codon_count}}) {
		
			#print "For codons with MAF of $min_codon_count, $num_alleles alleles: ".
			#	$min_count2num_alleles{$min_codon_count}->{$num_alleles} . " cases\n";
			
		}
	}
	
} # end product loop

close OUTFILE;








#########################################################################################
# NEW NONCODING DIVERSITY FILE
my @nc_A_count; # indices are alignment index (pos - 1)
my @nc_C_count;
my @nc_G_count;
my @nc_T_count;

my $aln_length = length($seqs_arr[0]);

# Initialize
for(my $seq_pos = 0; $seq_pos < $aln_length; $seq_pos++) {
	$nc_A_count[$seq_pos] = 0;
	$nc_C_count[$seq_pos] = 0;
	$nc_G_count[$seq_pos] = 0;
	$nc_T_count[$seq_pos] = 0;
}

my $num_seqs = scalar(@seqs_arr);

# Store nucleotide counts
for(my $seq_index = 0; $seq_index < $num_seqs; $seq_index++) {
	my $this_seq = $seqs_arr[$seq_index];
	
	for(my $seq_pos = 0; $seq_pos < $aln_length; $seq_pos++) {
		my $this_nt = substr($this_seq,$seq_pos,1);
		if($this_nt eq 'A') {
			$nc_A_count[$seq_pos]++;
		} elsif($this_nt eq 'C') {
			$nc_C_count[$seq_pos]++;
		} elsif($this_nt eq 'G') {
			$nc_G_count[$seq_pos]++;
		} elsif($this_nt eq 'T') {
			$nc_T_count[$seq_pos]++;
		}
	}
}

my @pi_by_site;

# Calculate and store non-coding (just plain) pi
open(OUTFILE_GENE_DIV,">>within\_group\_site\_results\.txt");

print OUTFILE_GENE_DIV "file\tsite\tmaj_nt\tmaj_nt_count\t".
	"pi\tcoverage\t".
	"A\tC\tG\tT\n";

for(my $seq_pos = 0; $seq_pos < $aln_length; $seq_pos++) {
	my $A = $nc_A_count[$seq_pos];
	my $C = $nc_C_count[$seq_pos];
	my $G = $nc_G_count[$seq_pos];
	my $T = $nc_T_count[$seq_pos];	
	
	my $maj_nt = 'A';
	my $maj_nt_count = $A;
	
	if($C > $maj_nt_count) {
		$maj_nt = 'C';
		$maj_nt_count = $C;
	}
	
	if($G > $maj_nt_count) {
		$maj_nt = 'G';
		$maj_nt_count = $G;
	}
	
	if($T > $maj_nt_count) {
		$maj_nt = 'T';
		$maj_nt_count = $T;
	}
	
	my $num_pw_diffs = ($A * $C) + ($A * $G) + ($A * $T) +
		($C * $G) + ($C * $T) + ($G * $T);
		
	my $num_pw_comps = ($num_seqs**2 - $num_seqs) / 2;
	
	my $mean_num_pw_diffs = $num_pw_diffs / $num_pw_comps; # this is pi
	
	$pi_by_site[$seq_pos] = $mean_num_pw_diffs;
	
	my $curr_site = $seq_pos + 1;
	
	print OUTFILE_GENE_DIV "$fasta_file_name\t$curr_site\t".
		"$maj_nt\t$maj_nt_count\t".
		"$mean_num_pw_diffs\t$num_seqs\t".
		"$A\t$C\t$G\t$T\n";
}
close OUTFILE_GENE_DIV;
########################################################################################



print "################################################################################";
print "\nTOTALS:\nN_sites: $total_N_sites\nS_sites: $total_S_sites\n".
	"N_diffs: $total_N_diffs\nS_diffs: $total_S_diffs\n";
print "################################################################################\n";

# FOR WNV, not for HPV
#foreach my $min_codon_count (sort {$a <=> $b} keys %min_count2num_alleles_TOTALS) {
#	foreach my $num_alleles (sort {$a <=> $b} keys %{$min_count2num_alleles_TOTALS{$min_codon_count}}) {
#		print "For codons with MAF of $min_codon_count, $num_alleles alleles: ".
#			$min_count2num_alleles_TOTALS{$min_codon_count}->{$num_alleles} . " cases\n";
#	}
#}

print "\n";

# Print a completion message to screen
&end_the_program;

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
				die "\n\n## WARNING: CDS annotation(s) in $gtf_file_name does not have a ".
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
# Get the amino acid (single-letter code) encoded by a given DNA or RNA codon
# Returns an array with:
#	returned[0] = number of nonsynonymous sites
#	returned[1] = number of synonymous sites
sub get_amino_acid {
	my ($codon) = @_;
	my $amino_acid;
	
	# Establish genetic code for use with synonymous sites; DNA or RNA
	my %code = (
		"AAA" => "K",
		"AAC" => "N",
		"AAG" => "K",
		"AAT" => "N",
		"AAU" => "N",
		"ACA" => "T",
		"ACC" => "T",
		"ACG" => "T",
		"ACT" => "T",
		"ACU" => "T",
		"AGA" => "R",
		"AGC" => "S",
		"AGG" => "R",
		"AGT" => "S",
		"AGU" => "S",
		"ATA" => "I",
		"ATC" => "I",
		"ATG" => "M",
		"ATT" => "I",
		"AUA" => "I",
		"AUC" => "I",
		"AUG" => "M",
		"AUU" => "I",
		"CAA" => "Q",
		"CAC" => "H",
		"CAG" => "Q",
		"CAT" => "H",
		"CAU" => "H",
		"CCA" => "P",
		"CCC" => "P",
		"CCG" => "P",
		"CCT" => "P",
		"CCU" => "P",
		"CGA" => "R",
		"CGC" => "R",
		"CGG" => "R",
		"CGT" => "R",
		"CGU" => "R",
		"CTA" => "L",
		"CTC" => "L",
		"CTG" => "L",
		"CTT" => "L",
		"CUA" => "L",
		"CUC" => "L",
		"CUG" => "L",
		"CUU" => "L",
		"GAA" => "E",
		"GAC" => "D",
		"GAG" => "E",
		"GAT" => "D",
		"GAU" => "D",
		"GCA" => "A",
		"GCC" => "A",
		"GCG" => "A",
		"GCT" => "A",
		"GCU" => "A",
		"GGA" => "G",
		"GGC" => "G",
		"GGG" => "G",
		"GGT" => "G",
		"GGU" => "G",
		"GTA" => "V",
		"GTC" => "V",
		"GTG" => "V",
		"GTT" => "V",
		"GUA" => "V",
		"GUC" => "V",
		"GUG" => "V",
		"GUU" => "V",
		"TAA" => "*",
		"TAC" => "Y",
		"TAG" => "*",
		"TAT" => "Y",
		"UAA" => "*",
		"UAC" => "Y",
		"UAG" => "*",
		"UAU" => "Y",
		"TCA" => "S",
		"TCC" => "S",
		"TCG" => "S",
		"TCT" => "S",
		"UCA" => "S",
		"UCC" => "S",
		"UCG" => "S",
		"UCU" => "S",
		"TGA" => "*",
		"TGC" => "C",
		"TGG" => "W",
		"TGT" => "C",
		"UGA" => "*",
		"UGC" => "C",
		"UGG" => "W",
		"UGU" => "C",
		"TTA" => "L",
		"TTC" => "F",
		"TTG" => "L",
		"TTT" => "F",
		"UUA" => "L",
		"UUC" => "F",
		"UUG" => "L",
		"UUU" => "F",
	);
	
	$amino_acid = $code{$codon};
	
	if($amino_acid eq '') {
		$amino_acid = '?';
	}
	
	return $amino_acid;
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
	
	open (CURRINFILE, $gtf_file_name);
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

	##### Numbers of nonsynonymous and synonymous sites for all three sites of every codon,
	# ignoring STOP codons per usual, stored in arrays by the codon's name in the format:
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