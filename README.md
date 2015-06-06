# snpgenie

Perl software for estimating evolutionary parameters from pooled next-generation sequencing single-nucleotide polymorphism data. 

# Introduction
SNPGenie is Perl-based program for calculating evolutionary parameters such as nonsynonymous and synonymous nucleotide diversities (πN and πS, respectively) from next-generation sequencing (NGS) SNP reports generated using CLC Genomics Workbench or Geneious. For the results to have biological meaning, the SNP report must contain SNPs called from the sequencing of a pooled nucleic acid sample that is representative of the population of interest. For example, if one is interested in determining the nucleotide diversity of a virus population within a single host, it would be appropriate to sequence the pooled nucleic acid content of the viruses in a blood sample from that host.

# Using snpgenie
SNPGenie version 1.2 is a command-line interface application written in Perl. As such, it is limited only by the memory and processing capabilities of the local hardware. As input, it accepts:

1. One or more **reference sequence** files in **FASTA** format (.fa/.fasta); 
2. One file with CDS information in **Gene Transfer Format** (.gtf); and 
3. One or more tab-delimited (.txt) **SNP Reports** in CLC format (http://www.geneious.com/). 

Further details on input are below.

## CLC Genomics Workbench Input
At minimum, a CLC SNP report must include the following 8 default column selections, with the unaltered CLC column headers: 
1. **Reference Position**, which refers to the start site of the polymorphism within the reference FASTA sequence;
2. **Type**, which refers to the nature of the record, usually the type of polymorphism, e.g., “SNV”;
3. **Reference**, the reference nucleotide(s) at that site(s);
4. **Allele**, the variant nucleotide(s) at that site(s);
5. **Count**, the number of reads containing the variant;
6. **Coverage**, the total number of sequencing reads at the site(s);
7. **Frequency**, the frequency of the variant as a percentage, e.g., “14.6” for 14.60%; and
8. **Overlapping annotations**, containing the name of the protein product or open reading frame (ORF), e.g., “CDS: ORF1”.

In addition to the aforementioned columns, the SNP report should also be free of thousand separators (,) in the Reference Position, Count, and Coverage columns (default format). The Variant Frequency must remain a percentage (default format). Finally, the user should verify that the reading frame in the CLC output is correct. SNPGenie will produce various errors to ensure that these things are so, e.g., by checking that all products begin with START and end with STOP codons, and checking for premature stop codons.

## Geneious Input

## Gene Transfer Input
	The Gene Transfer Format (.gtf) file must include records for all ORFs present in your SNP Report(s). If a single ORF has multiple segments with different coordinates, simply enter one line for each segment, using the same product name. SNPGenie for CLC can currently handle 2 segments per ORF. For more information about GTF, please visit <http://mblab.wustl.edu/GTF22.html>.
	Given the appropriate files, SNPGenie calculates nucleotide diversities for nonsynonymous and synonymous sites in a protein-coding sequence. Nucleotide diversity may be defined as the average number of nucleotide variants per nucleotide site for all pairwise comparisons. To distinguish between nonsynonymous and synonymous differences and sites, it is necessary to consider the codon context of each nucleotide in a sequence. This is why the user must submit the starting and ending sites of the coding regions in the .gtf file, along with the reference FASTA sequence file, to allow an accurate estimation of the number of nonsynonymous and synonymous sites for each codon by the Nei-Gojobori (1986) method. SNPGenie first splits the coding sequence into codons, each of which contains 3 sites. The software then determines the number of these sites which are nonsynonymous and synonymous by testing all possible changes at each site of every codon in the sequence. Because different nucleotide variants at the same site may lead to both nonsynonymous and synonymous polymorphisms, fractional sites occur frequently (e.g., only 2 of 3 possible nucleotide substitutions at the third position of AGA cause an amino acid change; thus, that site is considered 2/3 nonsynonymous and 1/3 synonymous). Next, the SNP report is consulted for the presence of variants to produce a revised estimate. Variants are incorporated through averaging weighted by their frequency. Although it is a rare occurrence, high levels of sequence variation may alter the number of nonsynonymous and synonymous sites in a particular codon, contributing to an altered picture of natural selection. 
Next, SNPGenie calculates the number of nucleotide differences for each codon in each ORF specified in the .gtf file. Calculating nucleotide diversity codon-by-codon enables sliding-window analyses that may help to pinpoint important nucleotide regions subject to varying forms of natural selection. SNPGenie determines the average number of pairwise differences as follows: for every variant in the SNP Report, the number of variants is calculated as the product of the variant’s relative frequency and the coverage at that site. For each variant nucleotide (up to 3 non-reference nucleotides), the number of variants is stored, and their sum is subtracted from the coverage to yield the reference’s absolute frequency. Next, for each pairwise nucleotide comparison at the site, it is determined whether the comparison represents a nonsynonymous or synonymous change. If the former, the product of their absolute frequencies contributes to the number of nonsynonymous differences; if the latter, it contributes to the number of synonymous differences. When comparing codons with more than one nucleotide difference, all possible mutational pathways are considered, per the methods of Nei & Gojobori (1986). The sum of pairwise differences is divided by the total number of pairwise comparisons at the codon (i.e., cC2, where c = coverage) to yield the number of differences per site of each type. This is calculated separately for nonsynonymous and synonymous comparisons. For further background, see Nelson & Hughes (2015).
To run SNPGenie, first download the appropriate script and place it in your system’s PATH, or simply in your working directory. Next, place your SNP Report (.txt), FASTA (.fa/.fasta), and GTF (.gtf) files in your working directory (Figure 1). Open the command line prompt or terminal and navigate to the directory containing these files (Figure 2). Then simply execute SNPGenie (Figure 3).
SNPGenie for CLC automatically detects the number of FASTA files in your working directory. If there is one FASTA file, then SNPGenie enters its ONE-SEQUENCE MODE, where it is assumed that all SNP Reports refer to the same reference sequence, contained in the FASTA. This is the most common application. For example, one FASTA file might contain the reference (e.g., consensus) sequence for a virus population, while the GTF file contains the coordinates of the ORFs in the viral genome, and the SNP Reports contains SNPs called from different samples against the aforementioned reference.
There are some situations in which multiple FASTA sequences may be necessary, e.g., the segmented genome of an influenza virus. If there are multiple FASTA sequences in the working directory, SNPGenie for CLC will enter its MULTI-SEQUENCE MODE, and will look for a separate FASTA file for each ORF. These files must begin with the name of the respective ORF, followed by an underscore (“_”). For example, if one product was named ORF1, a FASTA file would need to exist by the name of ORF1_x.fa, where x is any string of alphanumeric characters.
Along with the aforementioned WARNINGS, SNPGenie will also alert you to such aberrations as missing FASTA files, incorrectly named ORFs, and reference sequences which do not match what is reported in the SNP Report. In most instances, SNPGenie will terminate. If this occurs, you may consult either the Terminal (command line) screen, or else the new SNPGenie_Results directory (created after running SNPGenie in your working directory) for the SNPGenie_WARNINGS.txt file. Also see section 5. Troubleshooting below.

# Output
SNPGenie creates a new folder called SNPGenie_Results within the working directory. This contains 3 TAB-delimited files:
	WARNINGS.txt. This file alerts you to peculiarities observed in your data or any indicators that your annotations might be incorrect. All warnings are also printed to the Terminal (Unix) window.

	Nucleotide_diversity_results.txt. This file contains the nucleotide diversity results for all SNP Reports processed. The columns are:
	File. The SNP Report analyzed.
	Product. The protein product or ORF.
	Site. The site of the first nucleotide in this row’s codon, relative to the reference sequence.
	Codon. The codon to which this row applies.
	Num Overlap ORF Nts. The number of nucleotides within this codon (ranging from 0 to 3) which overlap another ORF. For example, if a codon is contained within ORF1 and ORF2, this value will be 3.
	Nonsyn Diffs. For the codons observed at this site, determined using variant and coverage data, this column gives the mean number of pairwise differences per site which are nonsynonymous (amino acid-altering). For codon comparisons in which ≥2 nucleotides differ, all possible mutational pathways are considered via an adaptation of the Nei & Gojobori (1986) method (Nelson & Hughes 2015).
	Syn Diffs. For all possible codons observed at this site, determined using variant and coverage data, this column gives the mean number of pairwise differences per site which are synonymous (not amino acid-altering). For codon comparisons in which ≥2 nucleotides differ, all possible mutational pathways are considered via an adaptation of the Nei & Gojobori (1986) method (Nelson & Hughes 2015).
	Nonsyn Sites. The number of nucleotide sites in the current codon which are nonsynonymous, determined using an adaptation of the Nei & Gojobori (1986) method (Nelson & Hughes 2015) which considers codons frequencies. The sum of this value and the value of Syn Sites for a codon sum to 3, as there are 3 nucleotide sites per codon.
	Syn Sites. The number of nucleotide sites in the current codon which are synonymous, determined using an adaptation of the Nei & Gojobori (1986) method (Nelson & Hughes 2015) which considers codons frequencies. The sum of this value and the value of Nonsyn Sites for a codon sum to 3, as there are 3 nucleotide sites per codon. 
	Nonsyn Sites (Ref). The number of nucleotide sites in the reference sequence of the current codon which are nonsynonymous., determined using the methods of Nei & Gojobori (1986). The sum of this value and the value of Syn Sites (Ref) will be 3.
	Syn Sites (Ref). The number of nucleotide sites in the reference sequence of the current codon which are synonymous, determined using the methods of Nei & Gojobori (1986). The sum of this value and the value of Nonsyn Sites (Ref) will be 3.
	πN. The mean number of nonsynonymous differences per nonsynonymous site in the population of sequences.
	πS. The mean number of nonsynonymous differences per nonsynonymous site in the population of sequences.
	Mean dN vs. Ref. Each individual sequence read can be compared to the reference sequence to yield the number of nonsynonymous differences. The Jukes-Cantor correction is applied to account for multiple mutations at the same site, yielding dN. The mean of all such comparisons is the value given in this column.
	Mean dS vs. Ref. Each individual sequence read can be compared to the reference sequence to yield the number of synonymous differences. The Jukes-Cantor correction is applied to account for multiple mutations at the same site, yielding dS. The mean of all such comparisons is the value given in this column.

	Gene_diversity_results.txt. This file contains the gene diversity results for all SNP Reports processed. The columns are:
	File. The SNP Report analyzed.
	Product. The protein product or ORF.
	Site. The site of the nucleotide in the reference sequence.
	Nucleotide. The nucleotide to which this row applies.
	Position in Codon. The position (1, 2, or 3) of the nucleotide within its codon. 
	Overlapping ORFs. The number of other ORFs (besides the one for which the nucleotide is reported) which overlap this nucleotide site. For example, if ORF1, ORF2, and ORF3 all overlap at a site, the value of this column will be 2 each time it is reported (once for ORF1, once for ORF2, and once for ORF3).
	Codon Start Site. The position of the first nucleotide in the relevant codon within the reference sequence.
	Codon. The codon in which the nucleotide exists.
	Gene Diversity. The value of gene diversity (expected heterozygosity, HO) at this site, defined as 1-∑_(i=1)^4▒x_i^2 , where xi is the frequency of nucleotide i at the site (1=A, 2=C, 3=G, 4=T).
	Site Classification vs. Ref. Possible categories for polymorphic sites are 1=synonymous, 2=nonsynonymous, 3=ambiguous (Knapp et al. 2011) as compared to the reference.
	Site Classification vs. All Codons. Possible categories for polymorphic sites are 1=synonymous, 2=nonsynonymous, 3=ambiguous (Knapp et al. 2011). All codons present at the site, given the variant data in the SNP Report, are compared to determine this classification.
	Site Classification vs. Ref. Cha.
	
Last paragraph.

# Geneious
SNPGenie SNPGenie version 1.2 is a command-line interface application written in Perl and also available as an executable (.exe) file. As such, it is limited only by the memory and processing capabilities of the local hardware. It accepts a reference sequence in FASTA format and a SNP Report in a comma-separated Geneious format (http://www.geneious.com/). At minimum, the SNP report must include the following default column selections, with the unaltered Geneious labels: 
(1) Type, which refers to the nature of the record entry, e.g., “Polymorphism”; 
(2) Minimum and (3) Maximum, which refer to the start and end sites of the polymorphism within the reference FASTA sequence, and will hold the same value for SNP records; 
(4) product, containing the name of the protein product or open reading frame, e.g., ORF1; 
(5) Amino Acid Change, which holds the single letter representations of the reference and variant amino acids, separated by spaces and a dash-greater-than (->) symbol, e.g., I -> V; N.B.: in the case of a synonymous SNP, this record is empty; 
(6) Codon Change and (7) Change, which contain the reference and variant nucleotides as full codons and single nucleotides, e.g., AUU -> GUU and A -> G, respectively, and are always populated for SNP records;
(8) Coverage, containing the number of sequencing reads that include the site; 
(9) Protein Effect, containing the word “Substitution” in the event of a nonsynonymous SNP, and “None” in the event of a synonymous SNP; and 
 (10) Variant Frequency, which contains the frequency of the nucleotide variant as a percentage, e.g., 14.60%.

In addition to the aforementioned columns, the SNP report must also be free of thousand separators (,) in the Minimum, Maximum, and Coverage columns. The Variant Frequency must be a percentage (default format). The line endings must be line feed (\n) or carriage return-line feed (\r\n) format. If an inappropriate line ending format is used, the program simply returns nucleotide diversities of zero. Mac OSX users may ensure an appropriate format by saving the SNP reports as Windows Comma Separated (.csv) files in Excel, or else manually changing the line ending format in an application such as TextWrangler. Finally, it is critically important for the user to verify that the reading frame in the Geneious output is correct by manually checking the START and STOP codons, as well as the intervening sequences for premature stop codons that may have been introduced by inaccurate data annotations.
	Given the appropriate files, SNPGenie calculates average nucleotide diversities at nonsynonymous and synonymous sites in a protein-coding sequence using its Sites Counter. Nucleotide diversity may be defined as the average number of nucleotide variants per nucleotide site for all pairwise comparisons. To distinguish between nonsynonymous and synonymous sites, it is necessary to consider the codon context of each nucleotide in a sequence. Thus, the user must submit the starting and ending sites of the coding region, along with the reference FASTA sequence, to allow an accurate estimation of the number of nonsynonymous and synonymous sites for each codon by the Nei-Gojobori (1986) method. SNPGenie first splits the coding sequence into codons, and then determines the number of possible changes at each site in every codon that would result in either nonsynonymous or synonymous polymorphisms. Each codon contains a total of 3 sites. Because different nucleotide variants at each site may lead to both nonsynonymous and synonymous polymorphisms, fractional sites occur frequently (e.g., only 2 of 3 possible nucleotide substitutions at the third position of AGA cause an amino acid change; thus, that site is considered 2/3 nonsynonymous and 1/3 synonymous). Next, the SNP report is consulted for the presence of variants to produce a revised estimate. Variant are incorporated through averaging weighted by their frequency. Although it is a rare occurrence, high levels of sequence variation may alter the number of nonsynonymous and synonymous sites in a particular codon, contributing to an altered picture of natural selection. 
To determine the number of sites, the user must execute the file SNPGenie_SiteCounter from the command line with the following arguments, separated by spaces: (1) reference sequence file; (2) SNP report file; (3) start site; and (4) end site, as follows:

SNPGenie_SiteCounter.exe MyRefSeq.fasta MySNPReport.csv start end

SNPGenie requires that the STOP codon be included in this range (i.e., the end site must be the final site of the STOP codon). The reference and SNP report files must be in the same directory as the executable file. Results are written to two output files, also in the same directory: one containing the number of nonsynonymous and synonymous sites by codon (file name suffix “_SitesByCodon”), and another containing the sums of such sites for the full sequence (file name suffix “_SumSites”). Because STOP codons (TAA, TAG, and TGA) are considered to contain no sites, they do not contribute to this estimate.
Next, SNPGenie calculates the number of nucleotide differences by codon and as a sum for the entire SNP report using its Differences Counter. Calculating nucleotide diversity codon-by-codon enables sliding-window analyses that may help to pinpoint important nucleotide regions subject to varying forms of natural selection. SNPGenie determines the average number of pairwise differences as follows: for every variant in the SNP Report, the number of variants is calculated as the product of the variant’s relative frequency and the coverage at that site. For each variant nucleotide (up to 3 non-reference nucleotides), the number of variants is stored, and their sum is subtracted from the coverage to yield the reference frequency. Next, for each pairwise nucleotide comparison at the site, it is determined whether the comparison represents a nonsynonymous or synonymous change. If the former, the product of their absolute frequencies contributes to the number of nonsynonymous differences; if the latter, it contributes to the number of synonymous differences. This sum is divided by the total number of pairwise comparisons at the site (i.e., cC2, where c = coverage) to yield the number of differences per site of each type:

((A*C)+(A*G)+(A*T)+(C*G)+(C*T)+(G*T))/(〖(coverage〗^2-coverage)/2)

where the numerator represents the total number of pairwise differences, and coverage is obtained from the SNP Report. Note that this is calculated separately for nonsynonymous and synonymous comparisons. Thus, if the change A->G is synonymous, it will be included when calculating synonymous diversity, but not nonsynonymous nucleotide diversity. Most terms in the numerator are 0 unless more than one SNP is present at a site.
To determine the number of differences, the user executes the file SNPGenie_DiffsCounter using the same arguments as the Sites Counter:

SNPGenie_DiffsCounter.exe MyRefSeq.fasta MySNPReport.csv start end

Two files are produced: one containing the nonsynonymous and synonymous differences annotated by codon (file name suffix “_DiffsByCodon”), and another containing the sums of such differences for all codons (file name suffix “_SumDiffs”). After this is completed, the results from the Difference Counter must be divided by the number of sites from the Site Counter to obtain the average number of pairwise differences per site of each type (yielding πN and πS). In other words, divide the output of SNPGenie_DiffsCounter.exe by the corresponding output of SNPGenie_SitesCounter.exe.

# Example
Suppose a blood sample is taken from an animal infected with an known influenza virus inoculum. This inoculum reference sequence is present in the file fluRef.fasta. Suppose Illumina sequencing has been performed on the pooled viral sample present in the blood sample, and Geneious has been used to generate a SNP report, fluPopSNPs.csv. Suppose further that we are interested in examining πN and πS in a protein product that is encoded by nucleotides 112 through 949. After consulting the SNP report to ensure it meets all format requirements, the following operations are performed:

SNPGenie_SitesCounter.exe fluRef.fasta fluPopSNPs.csv 112 949
SNPGenie_DiffsCounter.exe fluRef.fasta fluPopSNPs.csv 112 949

Finally, the first output is divided by the appropriate second output to obtain πN and πS.

# Citation
When using this software, please refer to and cite:

	Nelson CW, Hughes AL (2015) Within-host nucleotide diversity of virus populations: Insights from next-generation 		sequencing. Infection, Genetics and Evolution 30:1-7. 
	doi: 10.1016/j.meegid.2014.11.026

# Troubleshooting
•	Are the files tab-delimited? (All SNP Report and FASTA files must end in .txt with columns separated by tabs [\t])

•	Are (end-of-line) newline characters in either Unix LF (\n; preferable) or Windows CRLF (\r\n) format?

•	Are there no commas in the integer fields (CDS Position, Coverage, Minimum, and Maximum)? This can be readily fixed by formatting only those columns as type “General” in Excel.

•	Is the Variant Frequency column in a Percentage format—NOT a decimal (e.g., 0.1% rather than 0.001)?

•	Is the FRAME correct (i.e., do the codons in the SNP report begin with ATG and end with TAA, TAG, or TGA)?

# References
Knapp EW, Irausquin SJ, Friedman R, Hughes AL (2011) PolyAna: analyzing synonymous and nonsynonymous polymorphic sites. Conserv Genet Resour 3:429-431.
Nei M, Gojobori T (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Mol Biol Evol 3:418-26.
Nelson CW, Hughes AL (2015) Within-host nucleotide diversity of virus populations: Insights from next-generation sequencing. Infection, Genetics and Evolution 30:1-7. 
doi: 10.1016/j.meegid.2014.11.026


Copyright 2015
