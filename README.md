# snpgenie

Perl software for estimating evolutionary parameters from pooled next-generation sequencing single-nucleotide polymorphism data. 

# Introduction
SNPGenie is Perl-based program for calculating evolutionary parameters such as nonsynonymous and synonymous nucleotide diversities (πN and πS, respectively) from next-generation sequencing (NGS) SNP reports generated using CLC Genomics Workbench or Geneious. For the results to have biological meaning, the SNP report must contain SNPs called from the sequencing of a pooled nucleic acid sample that is representative of the population of interest. For example, if one is interested in determining the nucleotide diversity of a virus population within a single host, it would be appropriate to sequence the pooled nucleic acid content of the viruses in a blood sample from that host.

# Using snpgenie
SNPGenie version 1.2 is a command-line interface application written in Perl and also available as an executable (.exe) file. As such, it is limited only by the memory and processing capabilities of the local hardware. It accepts one or more reference sequence files in FASTA format (.fa/.fasta), one file with ORF information in Gene Transfer Format (.gtf), and one or more tab-delimited (.txt) SNP Reports in CLC format (http://www.geneious.com/). 
At minimum, the SNP report must include the following 8 default column selections, with the unaltered CLC labels: 
1. Reference Position, which refers to the start site of the polymorphism within the reference FASTA sequence;
2. Type, which refers to the nature of the record, usually the type of polymorphism, e.g., “SNV”;
3. Reference, the reference nucleotide(s) at that site(s);
4. Allele, the variant nucleotide(s) at that site(s);
5. Count, the number of reads containing the variant;
6. Coverage, the total number of sequencing reads at the site(s);
7. Frequency, the frequency of the variant as a percentage, e.g., “14.6” for 14.60%; and
8. Overlapping annotations, containing the name of the protein product or open reading frame (ORF), e.g., “CDS: ORF1”.

In addition to the aforementioned columns, the SNP report should also be free of thousand separators (,) in the Reference Position, Count, and Coverage columns (default format). The Variant Frequency must remain a percentage (default format). Finally, the user should verify that the reading frame in the CLC output is correct. SNPGenie will produce various errors to ensure that these things are so, e.g., by checking that all products begin with START and end with STOP codons, and checking for premature stop codons.
	The Gene Transfer Format (.gtf) file must include records for all ORFs present in your SNP Report(s). If a single ORF has multiple segments with different coordinates, simply enter one line for each segment, using the same product name. SNPGenie for CLC can currently handle 2 segments per ORF. For more information about GTF, please visit <http://mblab.wustl.edu/GTF22.html>.
	Given the appropriate files, SNPGenie calculates nucleotide diversities for nonsynonymous and synonymous sites in a protein-coding sequence. Nucleotide diversity may be defined as the average number of nucleotide variants per nucleotide site for all pairwise comparisons. To distinguish between nonsynonymous and synonymous differences and sites, it is necessary to consider the codon context of each nucleotide in a sequence. This is why the user must submit the starting and ending sites of the coding regions in the .gtf file, along with the reference FASTA sequence file, to allow an accurate estimation of the number of nonsynonymous and synonymous sites for each codon by the Nei-Gojobori (1986) method. SNPGenie first splits the coding sequence into codons, each of which contains 3 sites. The software then determines the number of these sites which are nonsynonymous and synonymous by testing all possible changes at each site of every codon in the sequence. Because different nucleotide variants at the same site may lead to both nonsynonymous and synonymous polymorphisms, fractional sites occur frequently (e.g., only 2 of 3 possible nucleotide substitutions at the third position of AGA cause an amino acid change; thus, that site is considered 2/3 nonsynonymous and 1/3 synonymous). Next, the SNP report is consulted for the presence of variants to produce a revised estimate. Variants are incorporated through averaging weighted by their frequency. Although it is a rare occurrence, high levels of sequence variation may alter the number of nonsynonymous and synonymous sites in a particular codon, contributing to an altered picture of natural selection. 
Next, SNPGenie calculates the number of nucleotide differences for each codon in each ORF specified in the .gtf file. Calculating nucleotide diversity codon-by-codon enables sliding-window analyses that may help to pinpoint important nucleotide regions subject to varying forms of natural selection. SNPGenie determines the average number of pairwise differences as follows: for every variant in the SNP Report, the number of variants is calculated as the product of the variant’s relative frequency and the coverage at that site. For each variant nucleotide (up to 3 non-reference nucleotides), the number of variants is stored, and their sum is subtracted from the coverage to yield the reference’s absolute frequency. Next, for each pairwise nucleotide comparison at the site, it is determined whether the comparison represents a nonsynonymous or synonymous change. If the former, the product of their absolute frequencies contributes to the number of nonsynonymous differences; if the latter, it contributes to the number of synonymous differences. When comparing codons with more than one nucleotide difference, all possible mutational pathways are considered, per the methods of Nei & Gojobori (1986). The sum of pairwise differences is divided by the total number of pairwise comparisons at the codon (i.e., cC2, where c = coverage) to yield the number of differences per site of each type. This is calculated separately for nonsynonymous and synonymous comparisons. For further background, see Nelson & Hughes (2015).
To run SNPGenie, first download the appropriate script and place it in your system’s PATH, or simply in your working directory. Next, place your SNP Report (.txt), FASTA (.fa/.fasta), and GTF (.gtf) files in your working directory (Figure 1). Open the command line prompt or terminal and navigate to the directory containing these files (Figure 2). Then simply execute SNPGenie (Figure 3).
SNPGenie for CLC automatically detects the number of FASTA files in your working directory. If there is one FASTA file, then SNPGenie enters its ONE-SEQUENCE MODE, where it is assumed that all SNP Reports refer to the same reference sequence, contained in the FASTA. This is the most common application. For example, one FASTA file might contain the reference (e.g., consensus) sequence for a virus population, while the GTF file contains the coordinates of the ORFs in the viral genome, and the SNP Reports contains SNPs called from different samples against the aforementioned reference.
There are some situations in which multiple FASTA sequences may be necessary, e.g., the segmented genome of an influenza virus. If there are multiple FASTA sequences in the working directory, SNPGenie for CLC will enter its MULTI-SEQUENCE MODE, and will look for a separate FASTA file for each ORF. These files must begin with the name of the respective ORF, followed by an underscore (“_”). For example, if one product was named ORF1, a FASTA file would need to exist by the name of ORF1_x.fa, where x is any string of alphanumeric characters.
Along with the aforementioned WARNINGS, SNPGenie will also alert you to such aberrations as missing FASTA files, incorrectly named ORFs, and reference sequences which do not match what is reported in the SNP Report. In most instances, SNPGenie will terminate. If this occurs, you may consult either the Terminal (command line) screen, or else the new SNPGenie_Results directory (created after running SNPGenie in your working directory) for the SNPGenie_WARNINGS.txt file. Also see section 5. Troubleshooting below.

Here's a list:

- Collapsing Toolbar
- FloatingActionButton
- View anchoring
- NavigationView
- Snackbar

Pre-requisites
--------------

- Android SDK v22
- Android Build Tools v22.0.1
- Android Support Repository v22.2

License
-------

Copyright 2014 The Android Open Source Project, Inc.

Licensed to the Apache Software Foundation (ASF) under one or more contributor
license agreements.  See the NOTICE file distributed with this work for
additional information regarding copyright ownership.  The ASF licenses this
file to you under the Apache License, Version 2.0 (the "License"); you may not
use this file except in compliance with the License.  You may obtain a copy of
the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
License for the specific language governing permissions and limitations under
the License.
