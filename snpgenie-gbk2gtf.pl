#! /usr/bin/perl

# Takes a Genbank file as its argument, and creates a GTF file for its products.

use strict;
#use warnings;
use IO::Handle;

my $genbank_file_name = $ARGV[0];

# Generate new file name prefix
my $new_file_prefix;
if($genbank_file_name =~/\.gbk/) { 
	$new_file_prefix = $`;
} else {
	$new_file_prefix = "inFile";
}

my $new_file_name = $new_file_prefix . ".gtf";

my $curr_ORF = '';
my $curr_start = '';
my $curr_stop = '';
my $curr_start2 = '';
my $curr_stop2 = '';

my %ORF_info_hh;
my @ORF_arr;

my $seen_label = 0;

open(GBK_FILE, $genbank_file_name) or die "Could not open file $genbank_file_name\n";

while(<GBK_FILE>) {
	#if($_ =~ /^\s*ORF\s+(\d+)\.\.(\d+)/) {
	if($_ =~ /^\s*CDS\s+(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	#} elsif($_ =~ /^\s*ORF\s+<(\d+)\.\.(\d+)/) {
	} elsif($_ =~ /^\s*CDS\s+<(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	#} elsif($_ =~ /^\s*ORF\s+join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)/) {
	} elsif($_ =~ /^\s*CDS\s+join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$curr_start2 = $3;
		$curr_stop2 = $4;
		$seen_label = 0;
	}
	
	if($_ =~ /^\s*ORF\s+(\d+)\.\.(\d+)/) {
	#if($_ =~ /^\s*CDS\s+(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	} elsif($_ =~ /^\s*ORF\s+<(\d+)\.\.(\d+)/) {
	#} elsif($_ =~ /^\s*CDS\s+<(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	} elsif($_ =~ /^\s*ORF\s+join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)/) {
	#} elsif($_ =~ /^\s*CDS\s+join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$curr_start2 = $3;
		$curr_stop2 = $4;
		$seen_label = 0;
	}
	
	if($_ =~ /^\s*misc_feature\s+(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	} elsif($_ =~ /^\s*misc_feature\s+<(\d+)\.\.(\d+)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$seen_label = 0;
	} elsif($_ =~ /^\s*misc_feature\s+join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)/) {
		$curr_start = $1;
		$curr_stop = $2;
		$curr_start2 = $3;
		$curr_stop2 = $4;
		$seen_label = 0;
	}
	
	#if($_ =~ /\s*\/label=(\w+)/) { # This was for the SIV and SHFVkrc1 and SHFVkrc2
	if($seen_label == 0) {
		if($_ =~ /\s*\/label="([\w\s']+ [\w\s']+)"/) { # try to for all non-whitespace names
			if($curr_start2 eq '') {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = '';
				$curr_stop = '';
				$curr_start2 = '';
				$curr_stop2 = '';
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label="([\w\s']+)"/) { # JUST ADDED THIS FOR MANESS DATA
			if($curr_start2 eq '') {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = '';
				$curr_stop = '';
				$curr_start2 = '';
				$curr_stop2 = '';
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/label=([\w\s']+)/) {
			if($curr_start2 eq '') {
				$curr_ORF = $1;
				#print "Here is \$curr_ORF: $curr_ORF";
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				#print "Here is \$curr_ORF: $curr_ORF";
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = '';
				$curr_stop = '';
				$curr_start2 = '';
				$curr_stop2 = '';
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product="([\w\s']+ [\w\s']+)"/) { # try to for all non-whitespace names
			if($curr_start2 eq '') {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = '';
				$curr_stop = '';
				$curr_start2 = '';
				$curr_stop2 = '';
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product="([\w\s']+)"/) {
			if($curr_start2 eq '') {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = '';
				$curr_stop = '';
				$curr_start2 = '';
				$curr_stop2 = '';
				
				$seen_label = 1;
			}
		} elsif($_ =~ /\s*\/product=([\w\s']+)/) {
			if($curr_start2 eq '') {
				$curr_ORF = $1;
				#print "Here is \$curr_ORF: $curr_ORF";
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				
				$seen_label = 1;
			} else {
				$curr_ORF = $1;
				#print "Here is \$curr_ORF: $curr_ORF";
				chomp($curr_ORF);
				push(@ORF_arr,$curr_ORF);
				$ORF_info_hh{$curr_ORF}->{start} = $curr_start;
				$ORF_info_hh{$curr_ORF}->{stop} = $curr_stop;
				$ORF_info_hh{$curr_ORF}->{start2} = $curr_start2;
				$ORF_info_hh{$curr_ORF}->{stop2} = $curr_stop2;
				
				$curr_ORF = '';
				$curr_start = '';
				$curr_stop = '';
				$curr_start2 = '';
				$curr_stop2 = '';
				
				$seen_label = 1;
			}
		}		
	}
}

close GBK_FILE;

my @ORF_arr_sorted = sort (@ORF_arr);

open(OUTFILE,">>$new_file_name");

foreach my $curr_product (@ORF_arr_sorted) {
	#print "$curr_product\n";

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



#print "\n\nProducts are:";

#foreach (@ORF_arr) {
#	print "\n\n" . $_ . "\nStart: " . $ORF_info_hh{$_}->{start} . "\nStop: " . $ORF_info_hh{$_}->{stop};
#	if (exists $ORF_info_hh{$_}->{start2}) {
#		print "\nStart2: " . $ORF_info_hh{$_}->{start2} . "\nStop2: " . $ORF_info_hh{$_}->{stop2};
#	}
#}


#	open (INFILE, $curr_snp_report_name);
#	while (<INFILE>) {
#		if($line == 0) {
#			$line++;
#		} else {
#			chomp;
#			
#			my @line_arr = split(/\t/,$_,-1);
#			my $type = $line_arr[$index_type];
#			my $min = $line_arr[$index_min];
#			my $max = $line_arr[$index_max];
#			my $product_name = $line_arr[$index_product];
#			
#			if($line_arr[$index_min] =~ /\<(\d+)/) {
#				$min = $1;
#			}
#			
#			# Store the CDS coordinates by the product name
#			if($type eq 'CDS' && $product_name ne '') { # &
#				$gtf_info_hh{$product_name}->{start} = $min;
#				$gtf_info_hh{$product_name}->{stop} = $max;
#			}
#			
#		}
#	}
#	close INFILE;

#if($line =~ /^\/\/\n/) { # If $line is the end of record, //\n
#	last;
#}

#STDOUT->autoflush(1);
#
#my @snp_report_file_names_arr = &get_snp_report_file_names;
#
#foreach my $curr_snp_report_name (@snp_report_file_names_arr) {
#	#print "The current report is: $curr_snp_report_name\n";
#	
#	# Generate new file name prefix
#	my $new_file_prefix;
#	if($curr_snp_report_name =~/\.txt/) { 
#		$new_file_prefix = $`;
#	} else {
#		$new_file_prefix = "inFile";
#	}
#	
#	my @header_names_arr = &get_header_names($curr_snp_report_name);
#	
#	my $index_min;
#	my $index_max;
#	my $index_type;
#	my $index_poly_type;
#	my $index_cds_position;
#	my $index_change;
#	my $index_percent;
#	my $index_cov;
#	my $index_product;
#	
#	# Determine the index of each column
#	for (my $i=0; $i<scalar(@header_names_arr); $i++) {
#		if ($header_names_arr[$i] eq 'Minimum') { # MUST MAKE SURE MIN = MAX
#			$index_min = $i;
#		} elsif ($header_names_arr[$i] eq 'Maximum') { # Will check the value BEGINS with 'SNP'
#			$index_max = $i;
#		} elsif ($header_names_arr[$i] eq 'Type') { # Will check the value BEGINS with 'SNP'
#			$index_type = $i;
#		} elsif ($header_names_arr[$i] eq 'Polymorphism Type') { # Will check the value BEGINS with 'SNP'
#			$index_poly_type = $i;
#		} elsif ($header_names_arr[$i] eq 'CDS Position') {
#			$index_cds_position = $i;
#		} elsif ($header_names_arr[$i] eq 'Change') {
#			$index_change = $i;
#		} elsif ($header_names_arr[$i] eq 'Variant Frequency') {
#			$index_percent = $i;
#		} elsif ($header_names_arr[$i] eq 'Coverage') {
#			$index_cov = $i;
#		} elsif ($header_names_arr[$i] eq 'product') {
#			$index_product = $i;
#		}
#	}
#	
#	my @product_names_arr = &get_product_names_geneious($curr_snp_report_name,$index_product,$index_type);
#	# N.B.: product names comes from both type eq "CDS" and type eq "Polymorphism", so
#	# we may need to add a check and print a warning if a product is known but does not
#	# have any CDS information listed
#	
#	my %gtf_info_hh;
#	
#	my $line = 0;
#	open (INFILE, $curr_snp_report_name);
#	while (<INFILE>) {
#		if($line == 0) {
#			$line++;
#		} else {
#			chomp;
#			
#			my @line_arr = split(/\t/,$_,-1);
#			my $type = $line_arr[$index_type];
#			my $min = $line_arr[$index_min];
#			my $max = $line_arr[$index_max];
#			my $product_name = $line_arr[$index_product];
#			
#			if($line_arr[$index_min] =~ /\<(\d+)/) {
#				$min = $1;
#			}
#			
#			# Store the CDS coordinates by the product name
#			if($type eq 'CDS' && $product_name ne '') { # &
#				$gtf_info_hh{$product_name}->{start} = $min;
#				$gtf_info_hh{$product_name}->{stop} = $max;
#			}
#			
#		}
#	}
#	close INFILE;
#	
#	my $out_gtf_file_name;
#
#	print "$curr_snp_report_name\n";
#
#	if($curr_snp_report_name =~/\.txt/) { 
#		$out_gtf_file_name = $` . '.gtf'; # $` is for the string BEFORE the match
#	} else {
#		$out_gtf_file_name = 'inFile.gtf';
#	}
#	
#	print "$out_gtf_file_name\n";
#	
#	open(OUTFILE,">>$out_gtf_file_name");
#	
#	my @current_products_arr = sort (keys %gtf_info_hh);
#	
#	foreach my $curr_product (@current_products_arr) {
#		print "$curr_product\n";
#		
#		my $this_min = $gtf_info_hh{$curr_product}->{start};
#		my $this_max = $gtf_info_hh{$curr_product}->{stop};
#		print OUTFILE "$curr_snp_report_name\tGeneious\tCDS\t$this_min\t$this_max\t\.\t\+\t0\tgene_id \"$curr_product\";\n";	
#	}	
#	close OUTFILE;
#	
#	print "\n";
##			PRINT, tab-delimited
##			REPORT NAME
##			CLC OR GENEIOUS
##			CDS
##			MIN
##			MAX
##			.
##			+
##			0
##			gene_id "";
#
#}
#
#
##########################################################################################
##########################################################################################
####################################### SUBROUTINES ######################################
##########################################################################################
##########################################################################################
#
#
##########################################################################################
## Obtains all file names in current directory ending in .txt
#sub get_snp_report_file_names { 
#	my @snp_report_file_names = glob "*.txt";
#	if (scalar (@snp_report_file_names) == 0) {
#		die "\n\n## WARNING:\n## There are no .txt SNP Reports.\n\n";
#	}
#	return 	@snp_report_file_names;
#}
#
##########################################################################################
#sub get_header_names {
#	my ($curr_snp_report_filename) = @_;
#	#print "\n$curr_snp_report_filename\n";
#	my $line = 0;
#	open (CURRINFILE, $curr_snp_report_filename);
#	while (<CURRINFILE>) {
#		if($line == 0) {
#			if(!($_ =~/\t/)) {
#				chdir('SNPGenie_Results');
#				open(ERROR_FILE,">>SNPGenie\_WARNINGS\.txt");
#				# FILE | PRODUCT | SITE | CODON | WARNING
#				print ERROR_FILE "$curr_snp_report_filename\tN/A\tN/A\t".
#					"File not TAB-delimited (\\t) or there is only one column. SNPGenie terminated\n";
#				close ERROR_FILE;
#				chdir('..');
#				
#				die "\n\n## WARNING:\n## The SNP Report $curr_snp_report_filename is ".
#					"not TAB-delimited (\\t), or there is only one column.\n\n";
#			}
#			chomp;
#			$line++;
#			my @line_arr = split(/\t/,$_);
#			return @line_arr;
#		}
#	}
#	close CURRINFILE;
#}
#
##########################################################################################
#sub get_product_names_geneious {
#	my ($curr_snp_report_filename,$index_product,$index_type) = @_;
#	my $line = 0;
#	my %products_hash;
#	open (CURRINFILE, $curr_snp_report_filename);
#	while (<CURRINFILE>) {
#		if ($line == 0) {
#			$line++;
#		} else {
#			chomp;
#			my @line_arr = split(/\t/,$_);
#			my $product = $line_arr[$index_product];
#			my $type = $line_arr[$index_type];
#			
#			if($product ne '') {
#				if ( ! exists $products_hash{$product} && 
#					(($type eq 'Polymorphism') || ($type eq 'CDS')) ) {
#						$products_hash{$product} = 1;
#				}
#			}
#		}
#	}
#	close CURRINFILE;
#	my @product_names = keys %products_hash;
#	return @product_names;
#}