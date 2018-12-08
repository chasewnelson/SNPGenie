<img src="https://github.com/chasewnelson/snpgenie/blob/master/logo1a.jpg?raw=true" alt="SNPGenie logo" width="450" height="175" align="middle">

SNPGenie is a collection of Perl scripts for estimating *π*<sub>N</sub>/*π*<sub>S</sub>, *d*<sub>N</sub>/*d*<sub>S</sub>, and other evolutionary parameters from next-generation sequencing (NGS) single-nucleotide polymorphism (SNP) variant data. Each type of analysis uses its own script:

1. **WITHIN-POOL ANALYSIS**. Use **snpgenie.pl**, the original SNPGenie. This will perform an analysis of within-sample *π*<sub>N</sub>/*π*<sub>S</sub> from **pooled NGS SNP data**. SNP reports (acceptable in a variety of formats) must each correspond to a single population, with variants called relative to a single reference sequence (one sequence in one FASTA file). For example, call:

        snpgenie.pl --vcfformat=4 --snpreport=<variants>.vcf --fastafile=<ref_seq>.fa --gtffile=<CDS_annotations>.gtf

2. **WITHIN-GROUP ANALYSIS**. Use **snpgenie\_within\_group.pl**. This will perform a traditional within-group analysis of *d*<sub>N</sub>/*d*<sub>S</sub> among aligned sequences in a FASTA file, such as possible using the <a target="_blank" href="http://www.megasoftware.net/">MEGA software</a>, with automation, sliding windows, parallelism, and overlapping gene features. For example, call:

        snpgenie_within_group.pl --fasta_file_name=<aligned_seqs>.fa --gtf_file_name=<CDS_annotations>.gtf --num_bootstraps=10000 --procs_per_node=8

3. **BETWEEN-GROUP ANALYSIS**. For a traditional between-group analysis comparing different groups of aligned sequences in two or more FASTA files, such as possible using the <a target="_blank" href="http://www.megasoftware.net/">MEGA software</a>, with automation, sliding windows, parallelism, and overlapping gene features. For example, call:

        snpgenie_between_group.pl

**UPDATE FOR VCF INPUT:** given myriad VCF formats in use, it is **NECESSARY TO SPECIFY THE FORMAT OF THE VCF** SNP report input using the **--vcfformat** argument. Format 4 is the most common. For details, see the section on [VCF](#vcf).

## Contents

* [Introduction](#introduction)
* [SNPGenie](#snpgenie)
	* [SNPGenie Input](#snpgenie-input)
		* **Input 1:** [Reference Sequence](#ref-seq)
		* **Input 2:** [Gene Transfer Format](#gtf)
		* **Input 3:** [SNP Report(s)](#SNP-Reports)
			* [CLC](#clc)
			* [Geneious](#geneious)
			* [VCF](#vcf)
		* [A Note on Reverse Complement ('—' Strand) Records](#revcom)
	* [Options](#options)
	* [Output](#output)
* [SNPGenie Within-Group](#snpgenie-within)
* [SNPGenie Between-Group](#snpgenie-between)
* [How SNPGenie Works](#how-snpgenie-works)
* [Additional Scripts](#additional-scripts)
* [Troubleshooting](#troubleshooting)
* [Citation](#citation)
* [Studies Using SNPGenie](#studies-using-snpgenie)
* [Contact](#contact)
* [References](#references)


## <a name="introduction"></a>Introduction

New applications of next-generation sequencing (NGS) use pooled samples containing DNA from multiple individuals to perform population genetic analyses. SNPGenie is a Perl program which analyzes the variant results of a single-nucleotide polymorphism (SNP) caller to calculate evolutionary parameters, such as nucleotide diversity (including its nonsynonymous and synonymous partitions, *π*<sub>N</sub> and *π*<sub>S</sub>) and gene diversity. Variant calls are typically present in annotation tables and assume that the pooled DNA sample is representative of some population of interest. For example, if one is interested in determining the nucleotide diversity of a virus population within a single host, it might be appropriate to sequence the pooled nucleic acid content of the virus in a blood sample from that host. Comparing *π*<sub>N</sub> and *π*<sub>S</sub> for a gene product, or comparing gene diversity at polymorphic sites of different types, may help to identify genes or gene regions subject to positive (Darwinian) selection, negative (purifying) selection, and random genetic drift. SNPGenie also includes such features as minimum allele frequency trimming (see [Options](#options)), and can be combined with upstream applications such as maximum-likelihood SNP calling techniques (*e.g.*, see Lynch *et al.* 2014). For additional background, see Nelson & Hughes (2015) in the [References](#references).

## <a name="snpgenie"></a>SNPGenie

### <a name="snpgenie-input"></a>SNPGenie Input

SNPGenie is a command-line application written in Perl, with no additional dependencies beyond the standard Perl package (*i.e.*, simply download and run). Input includes:

1. One [**Reference Sequence**](#ref-seq) file, containing one reference sequence, in **FASTA** format (.fa/.fasta); 
2. One file with CDS annotations in [**Gene Transfer Format**](#gtf) (.gtf); and 
3. One or more tab-delimited (.txt) [**SNP Reports**](#SNP-Reports), each corresponding to a single pooled-sequencing run (*i.e.*, sample from a population) in CLC, VCF, or Geneious formats, with variants called relative to the reference sequence.

**All input files must be placed in the same directory. Then simply run the appropriate SNPGenie script from the command line** (see [Options](#options) if you wish for more control). To do this, first download the **snpgenie.pl** script and place it in your system’s $PATH, or simply in your working directory. Next, place your SNP report(s), FASTA (.fa/.fasta), and GTF (.gtf) files in your working directory. Open the command line prompt (or Terminal) and navigate to the directory containing these files using the "cd" command in your shell. Finally, execute SNPGenie by typing the name of the script and pressing the \<RETURN\> (\<ENTER\>) key:

    snpgenie.pl

Further details on [input](#ref-seq) and [options](#options) are below.

#### <a name="ref-seq"></a>INPUT (1) One Reference Sequence File
Only one reference sequence must be provided in a single **FASTA** (.fa/.fasta) file (*i.e.*, one file containing one sequence). Thus, all SNP positions in the SNP reports are called relative to the single reference sequence. Because of this one-sequence stipulation, a script has been provided to split a multi-sequence FASTA file into its constituent sequences for multiple or parallel analyses; see [Additional Scripts](#additional-scripts) below.

#### <a name="gtf"></a>INPUT (2) One Gene Transfer Format File
The **Gene Transfer Format** (.gtf) file is tab (\t)-delimited, and must include *non-redundant* records for all CDS elements (*i.e.*, open reading frames, or ORFs) present in your SNP report(s). Thus, if multiple records exist with the same gene name (*e.g.*, different transcripts), ONLY ONE may be included in your analysis at a time. The GTF should also include any ORFs which do not contain any variants (if they exist). **SNPGenie expects every coding element to be labeled as type "CDS", and for its product name to follow a "gene\_id" tag.** In the case of CLC and Geneious SNP reports, this name must match that present in the SNP report. If a single coding element has multiple segments (*e.g.*, exons) with different coordinates, simply enter one line for each segment, using the same product name; they will be concatenated. Finally, for sequences with genes on the reverse '–' strand, SNPGenie must be run twice, once for each strand, with the minus strand's own set of input files (*i.e.*, the '–' strand FASTA, GTF, and SNP report); see [A Note on Reverse Complement ('–' Strand) Records](#revcom) below. For more information about GTF, please visit <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">The Brent Lab</a>. To convert a GFF file to GTF format, use a tool like <a target="_blank" href="https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread">gffread</a>. A short GTF example follows:

	reference.gbk	CLC	CDS	5694	8369	.	+	0	gene_id "ORF1";
	reference.gbk	CLC	CDS	8203	8772	.	+	0	gene_id "ORF2";
	reference.gbk	CLC	CDS	1465	4485	.	+	0	gene_id "ORF3";
	reference.gbk	CLC	CDS	5621	5687	.	+	0	gene_id "ORF4";
	reference.gbk	CLC	CDS	7920	8167	.	+	0	gene_id "ORF4";
	reference.gbk	CLC	CDS	5395	5687	.	+	0	gene_id "ORF5";
	reference.gbk	CLC	CDS	7920	8016	.	+	0	gene_id "ORF5";
	reference.gbk	CLC	CDS	4439	5080	.	+	0	gene_id "ORF6";
	reference.gbk	CLC	CDS	5247	5549	.	+	0	gene_id "ORF7";
	reference.gbk	CLC	CDS	4911	5246	.	+	0	gene_id "ORF8";

#### <a name="SNP-Reports"></a>INPUT (3) SNP Report File(s)

Each SNP report should contain variant calls for a single pooled-sequencing run (*i.e.*, sample from population), or alternatively a summary of multiple individual sequences. Its site numbers should refer to the exact reference sequence present in the FASTA input. SNP reports must be provided in one of the following formats:

##### <a name="clc"></a>CLC Genomics Workbench
A <a target="_blank" href="http://www.clcbio.com/products/clc-genomics-workbench/">CLC Genomics Workbench</a> SNP report must include the following columns to be used with SNPGenie, with the unaltered CLC column headers: 

* **Reference Position**, which refers to the start site of the polymorphism within the reference FASTA sequence;
* **Type**, which refers to the nature of the record, usually the type of polymorphism, *e.g.*, "SNV” for single-nucleotide variants;
* **Reference**, the reference nucleotide(s) at that site(s);
* **Allele**, the variant nucleotide(s) at that site(s);
* **Count**, the number of reads containing the variant;
* **Coverage**, the total number of sequencing reads at the site(s);
* **Frequency**, the frequency of the variant as a percentage, *e.g.*, “14.6” for 14.60%; and
* **Overlapping annotations**, containing the name of the protein product or open reading frame (ORF), *e.g.*, “CDS: ORF1”.

In addition to containing the aforementioned columns, the SNP report should ideally be free of thousand separators (,) in the Reference Position, Count, and Coverage columns (default .txt output). The Frequency must remain a percentage (default .txt output). Finally, the user should verify that the reading frame in the CLC input file is correct.

##### <a name="geneious"></a>Geneious
A <a target="_blank" href="http://www.geneious.com/">Geneious</a> SNP report must include the following columns to be used with SNPGenie, with the unaltered Geneious column headers:

* **Minimum** and **Maximum**, which refer to the start and end sites of the polymorphism within the reference FASTA sequence (for SNPs, both will be the same site number);
* **CDS Position**, the nucleotide position relative to the start site of the relevant CDS annotation;
* **Type**, which refers to the nature of the record entry, *e.g.*, “Polymorphism”;
* **Polymorphism Type**, which gives the type of polymorphism;
* **product**, containing the name of the protein product or open reading frame, e.g., ORF1; 
* **Change**, which contains the reference and variant nucleotides, *e.g.*, "A -> G", and is always populated for SNP records;
* **Coverage**, containing the number of sequencing reads that overlap the site; and
* **Variant Frequency**, which contains the frequency of the nucleotide variant as a percentage, *e.g.*, 14.60%.

In addition to containing the aforementioned columns, the SNP report should ideally be free of extraneous characters such as thousand separators (,). The Variant Frequency must remain a percentage (default .txt output). Finally, the user should verify that the reading frame in the Geneious output is correct.

##### <a name="vcf"></a>Variant Call Format (VCF)
A <a target="_blank" href="https://github.com/samtools/hts-specs">VCF</a> SNP report must include the following columns, with the unaltered VCF column headers:

* **CHROM**, the name of the reference genome;
* **POS**, which refers to the start site of the polymorphism within the reference FASTA sequence;
* **REF**, the reference nucleotide(s) at that site(s);
* **ALT**, the variant nucleotide(s) at that site(s), with multiple variants at the same site separated by commas (,);
* **QUAL**, the Phred quality score for the variant;
* **FILTER**, the filter status, based on such metrics as minimum frequencies and minimum quality scores;
* **INFO**, information regarding read depth and/or allele frequencies; 
* **FORMAT**, sometimes an alternative to the **INFO** column for read depth and/or allele frequency data; 
* **\<SAMPLE\>**, column(s) named by sample, a common alternative to the **INFO** column for read depth and/or allele frequency data.

Because VCF files vary widely in format, **SNPGenie now requires** users to specify exactly which VCF format is being used with the **--vcfformat** argument. New formats are being added on a case-by-case basis; users should [contact the author](#contact) to have new formats incorporated. **FORMAT 4 is the most common, and necessary if your file contains multiple \<SAMPLE\> columns**. All supported formats are:

1. **FORMAT (1): --vcfformat=1**. Multiple individual genomes have been sequenced separately, with SNPs summarized in the VCF file. SNPGenie will require the following in the **INFO** column:
	* **AN**, the number of alleles sequenced (*e.g.*, "AN=60"). If absent, will use **NS**, the number of samples (*i.e.*, individual sequencing experiments) being summarized (*e.g.*, "NS=30");
	* **AF**, the fractional (decimal) allele frequency(-ies) for the variant alleles in the same order as listed in the ALT column (*e.g.*, "AF=0.200"). If multiple variants exist at the same site, their frequencies must be separated by commans (*e.g.*, "AF=0.200,0.087").

2. **FORMAT (2):  --vcfformat=2**. Variants have been called from a pooled deep-sequencing sample (genomes from multiple individuals sequenced together). SNPGenie will require the following in the **INFO** column:
	* **DP**, the coverage or total read depth (*e.g.*, "DP=4249");
	* **AF**, the fractional (decimal) allele frequency(-ies) of the variant alleles in the same order as listed in the ALT column (*e.g.*, "AF=0.01247"; "AF=0.01247,0.08956" for two variants; etc.).

3. **FORMAT (3): --vcfformat=3**. Like format 2, variants have been called from a pooled deep-sequencing sample (genomes from multiple individuals sequenced together). Unlike format 2, SNPGenie will require the following in the **INFO** column:
	* **DP4**, containing the number of reference and variant reads on the forward and reverse strands in the format *DP4=\<num. fw ref reads\>,\<num. rev ref reads\>,\<num. fw var reads\>,\<num. rev var reads\>* (*e.g.*, "DP4=11,9,219,38").
	* **N.B.**: if multiple single nucleotide variants exist at the same site, this VCF format does not specify the number of reads for each variant. SNPGenie thus approximates by dividing the number of non-reference reads evenly amongst the variant alleles.

4. **FORMAT (4):  --vcfformat=4**. **THE MOST COMMON FORMAT.** Like formats 2 and 3, variants have been called from a pooled deep-sequencing sample (genomes from multiple individuals sequenced together). For this format, SNPGenie will require AD and DP data in the **FORMAT** and **\<SAMPLE\>** columns, where FORMAT is the final column header before the **\<SAMPLE\>** column(s) begins. The order of the data keys in the FORMAT column and the data values in the **\<SAMPLE\>** columns must match.
	* For reference allele depth, include the **AD** tag in the **FORMAT** column, which refers to the read depth value for each allele in the **\<SAMPLE\>** column(s). Values for variant allele(s) must be in the same order as listed in the **ALT** column (*e.g.*, "AD" in the FORMAT column and "75,77" in the \<SAMPLE\> column);
	* For coverage (total read depth), include **DP** in the **FORMAT** column, which refers to the total read depth in the **\<SAMPLE\>** column(s) (*e.g.*, "DP" in the FORMAT column and "152" in the \<SAMPLE\> column). If the value of **DP** does not equal the sum of reads for all variant nucleotides (sum of **AD** values), the latter will be used.

As usual, you will want to make sure to maintain the VCF file's features, such as TAB(\t)-delimited columns. Unlike some other formats, the allele frequency in VCF is a decimal.

#### <a name="revcom"></a>A Note on Reverse Complement ('–' Strand) Records
Many large genomes have coding products on both strands. In this case, SNPGenie must be run twice: once for the '+' strand, and once for the '—' strand. This requires FASTA, GTF, and SNP report input for the '–' strand. I provide a script for converting input files to their reverse complement strand in the [Additional Scripts](#additional-scripts) below. Note that, regardless of the original SNP report format, the reverse complement SNP report is in a CLC-like format that SNPGenie will recognize. For both runs, the GTF should include all products for both strands, with products on the strand being analyzed labeled '+' and coordinates defined with respect to the beginning of that strand's FASTA sequence. Also note that a GTF file containing *only* '—' strand records will not run; SNPGenie does calculations only for the products on the current + strand, using the '—' strand products only to determine the number of overlapping reading frames for each variant.

### <a name="options"></a>Options

To alter how SNPGenie works, the following options may be used:

* **--minfreq**: the minimum allele (SNP) frequency to include. Enter as a proportion/decimal (*e.g.*, 0.01), **not** as a percentage (*e.g.*, 1.0%). Default: 0.
* **--snpreport**: the (one) SNP report to analyze. Default: auto-detect all .txt and .csv file(s).
* **--vcfformat**: the exact format of VCF input. See [VCF](#vcf).
* **--fastafile**: the (one) reference sequence. Default: auto-detect a .fa or .fasta file.
* **--gtffile**: the (one) file with CDS annotations. Default: auto-detect the .gtf file.
* **--slidingwindow**: the length of the (optional) sliding (codon) window used in the analysis. Default: 9 codons.

For example, if you wanted to specify a minimum allele frequency of 1% and specify all three input files, you could run the following:

	snpgenie.pl --minfreq=0.01 --snpreport=mySNPreport.txt --fastafile=myFASTA.fa --gtffile=myGTF.gtf



### <a name="output"></a>Output

SNPGenie creates a new folder called **SNPGenie_Results** within the working directory. This contains the following TAB-delimited results files:

1. **SNPGenie\_parameters.txt**, containing the input parameters and file names.

2. **SNPGenie\_LOG.txt**, documenting any peculiarities or errors encountered. Warnings are also printed to the Terminal (shell) window.

3. **site\_results.txt**, providing results for all polymorphic sites. Note that, if the population is genetically homogenous at a site, even if it differs from the reference or ancestral sequence, it will not be considered polymorphic. Also keep in mind that columns are sorted by product first, then site number, with noncoding sites at the end of the file. Columns are:
	* *file*. The SNP report analyzed.
	* *product*. The CDS annotation to which the site belongs; "noncoding" if none on this strand.
	* *site*. The site coordinate of the nucleotide in the reference sequence.
	* *ref\_nt*. The identity of the nucleotide in the reference sequence.
	* *maj\_nt*. The most common nucleotide in the population at this site.
	* *position\_in\_codon*. If present in a CDS annotation, the position of this site within its codon (1, 2, or 3).
	* *overlapping\_ORFs*. The number of additional CDS annotations overlapping this site. For example, if the site is part of only one open reading frame, the value will be 0. If the site is part of two open reading frames, the value will be 1.
	* *codon\_start\_site*. The site coordinate of the relevant codon's first nucleotide in the reference sequence.
	* *codon*. The identity of the relevant codon.
	* *pi*. Nucleotide diversity at this site.
	* *gdiv*. Gene diversity at this site.
	* *class\_vs\_ref*. This site's classification, as compared to the reference sequence. For example, if the site contains only one SNP, and that SNP is synonymous, the site will be classified as Synonymous. Nonsynonymous or Synonymous.
	* *class*. This site's classification as compared to all sequences present in the population. For example, if the population contains both A and G residues at the third site of a GAA (reference) codon, then the site will be Synonymous, because both GAA and GAG encode Glu. On the other hand, if the population also contains a C at this site, the site will be Ambiguous, because GAC encodes Asp, meaning both nonsynonymous and synonymous polymorphisms exist at the site. Nonsynonymous, Synonymous, or Ambiguous.
	* *coverage*. The NGS read depth at the site.
	* *A*. The number of reads containing an A (adenine) nucleotide at this site. N.B.: may be fractional if the coverage and variant frequency given in the SNP report do not imply a whole number.
	* *C*. For C (cytosine), as for A.
	* *G*. For G (guanine), as for A.
	* *T*. For T (thymine), as for A.

4. **codon\_results.txt**, providing results for all polymorphic sites. Columns are:
	* *file*. The SNP report analyzed.
	* *product*. The CDS annotation to which the site belongs; "noncoding" if none.
	* *site*. The site coordinate of the nucleotide in the reference sequence.
	* *codon*. The identity of the relevant codon.
	* *num\_overlap\_ORF\_nts*. The number of nucleotides in this codon (up to 3) which overlap other ORFs (in addition to the current "product" annotation).
	* *N\_diffs*. The mean number of pairwise nucleotide comparisons in this codon which are nonsynonymous (i.e., amino acid-altering) in the pooled sequence sample. The numerator of *π*<sub>N</sub>.
	* *S\_diffs*. The mean number of pairwise nucleotide comparisons in this codon which are synonymous (i.e., amino acid-conserving) in the pooled sequence sample. The numerator of *π*<sub>S</sub>.
	* *N\_sites*. The mean number of sites in this codon which are nonsynonymous, given all sequences in the pooled sample. The denominator of *π*<sub>N</sub> and mean *d*<sub>N</sub> versus the reference.
	* *S\_sites*. The mean number of sites in this codon which are synonymous, given all sequences in the pooled sample. The denominator of *π*<sub>S</sub> mean *d*<sub>S</sub> versus the reference.
	* *N\_sites\_ref*. The number of sites in this codon which are nonsynonymous in the reference sequence.
	* *S\_sites\_ref*. The number of sites in this codon which are synonymous in the reference sequence.
	* *N\_diffs\_vs\_ref*. This codon's mean number of nonsynonymous nucleotide differences from the reference sequence in the pooled sequence sample. The numerator of mean *d*<sub>N</sub> versus the reference.
	* *S\_diffs\_vs\_ref*. This codon's mean number of synonymous nucleotide differences from the reference sequence in the pooled sequence sample. The numerator of mean *d*<sub>S</sub> versus the reference.
	* *gdiv*. Mean gene diversity (observed heterozygosity) for this codon's nucleotide sites.
	* *N\_gdiv*. Mean gene diversity for this codon's nonsynonymous polymorphic sites.
	* *S\_gdiv*. Mean gene diversity for this codon's synonymous polymorphic sites.

5. **product\_results.txt**, providing results for all CDS elements present in the GTF file for the '+' strand. Columns are:
	* *file*. The SNP report analyzed.
	* *product*. The CDS annotation to which the site belongs; "noncoding" if none.
	* *N\_diffs*. The sum over all codons in this product of the mean number of pairwise nucleotide comparisons which are nonsynonymous (i.e., amino acid-altering) in the pooled sequence sample. The numerator of *π*<sub>N</sub>.
	* *S\_diffs*. The sum over all codons in this product of the mean number of pairwise nucleotide comparisons which are synonymous (i.e., amino acid-conserving) in the pooled sequence sample. The numerator of *π*<sub>S</sub>.
	* *N\_diffs\_vs\_ref*. The sum over all codons in this product of the mean number of nonsynonymous nucleotide differences from the reference sequence in the pooled sequence sample. The numerator of mean *π*<sub>N</sub> versus the reference.
	* *S\_diffs\_vs\_ref*. The sum over all codons in this product of the mean number of synonymous nucleotide differences from the reference sequence in the pooled sequence sample. The numerator of mean *π*<sub>S</sub> versus the reference.
	* *N\_sites*. The mean number of sites in this product which are nonsynonymous, given all sequences in the pooled sample. The denominator of *π*<sub>N</sub> and mean *d*<sub>N</sub> versus the reference.
	* *S\_sites*. The mean number of sites in this product which are synonymous, given all sequences in the pooled sample. The denominator of *π*<sub>S</sub> and mean *d*<sub>S</sub> versus the reference.
	* *piN*. (*π*<sub>N</sub>.) The mean number of pairwise nonsynonymous differences per nonsynonymous site in this product.
	* *piS*. (*π*<sub>S</sub>.) The mean number of pairwise synonymous differences per synonymous site in this product.
	* *mean\_dN\_vs\_ref*. The mean number of nonsynonymous differences from the reference per nonsynonymous site in this product.
	* *mean\_dS\_vs\_ref*. The mean number of synonymous differences from the reference per synonymous site in this product.
	* *mean\_gdiv\_polymorphic*. Mean gene diversity (observed heterozygosity) at all polymorphic nucleotide sites in this product.
	* *mean\_N\_gdiv*. Mean gene diversity at all nonsynonymous polymorphic nucleotide sites in this product.
	* *mean\_S\_gdiv*. Mean gene diversity at all synonymous polymorphic nucleotide sites in this product.

6. **population\_summary.txt**, providing summary results for each population's sample (SNP report) with respect to the '+' strand. Columns are:
	* *file*. The SNP report analyzed.
	* *sites*. Total number of sites in the reference genome.
	* *sites\_coding*. Total number of sites in the reference genome which code for a protein product on the analyzed '+' strand, given the CDS annotations in the GTF file.
	* *sites\_noncoding*. Total number of sites in the reference genome which do not code for a protein product, given the CDS annotations in the GTF file.
	* *pi*. Mean number of pairwise differences per site in the pooled sample across the whole genome.
	* *pi\_coding*. Mean number of pairwise differences per site in the pooled sample across all coding sites in the genome.
	* *pi\_noncoding*. Mean number of pairwise differences per site in the pooled sample across all noncoding sites in the genome.
	* *N\_sites*. The mean number of sites in the genome which are nonsynonymous, given all sequences in the pooled sample. The denominator of *π*<sub>N</sub> and mean *d*<sub>N</sub> versus the reference.
	* *S\_sites*. The mean number of sites in the genome which are synonymous, given all sequences in the pooled sample. The denominator of *π*<sub>S</sub> and mean *d*<sub>S</sub> versus the reference.
	* *piN*. The mean number of pairwise nonsynonymous differences per nonsynonymous site across the genome of the pooled sample.
	* *piS*. The mean number of pairwise synonymous differences per synonymous site across the genome of the pooled sample.
	* *mean\_dN\_vs\_ref*. The mean number of nonsynonymous differences from the reference per nonsynonymous site across the genome of the pooled sample.
	* *mean\_dS\_vs\_ref*. The mean number of synonymous differences from the reference per synonymous site across the genome of the pooled sample.
	* *mean\_gdiv\_polymorphic*. Mean gene diversity (observed heterozygosity) at all polymorphic nucleotide sites in the genome of the pooled sample.
	* *mean\_N\_gdiv*. Mean gene diversity at all nonsynonymous polymorphic nucleotide sites in the genome of the pooled sample.
	* *mean\_S\_gdiv*. Mean gene diversity at all synonymous polymorphic nucleotide sites in the genome of the pooled sample.
	* *mean\_gdiv*. Mean gene diversity at all nucleotide sites in the genome of the pooled sample.
	* *sites\_polymorphic*. The number of sites in the genome of the pooled sample which are polymorphic.
	* *mean\_gdiv\_coding\_poly*. Mean gene diversity at all polymorphic nucleotide sites in the genome of the pooled sample which code for a protein product, given the CDS annotations in the GTF file.
	* *sites\_coding\_poly*. The number of sites in the genome of the pooled sample which are polymorphic and code for a protein product, given the CDS annotations in the GTF file.
	* *mean\_gdiv\_noncoding\_poly*. Mean gene diversity at all polymorphic nucleotide sites in the genome of the pooled sample which do not code for a protein product, given the CDS annotations in the GTF file.
	* *sites\_noncoding\_poly*. The number of sites in the genome of the pooled sample which are polymorphic and do not code for a protein product, given the CDS annotations in the GTF file.

7. **sliding\_window\_length\<Length\>\_results.txt**, containing codon-based results over a sliding window, with a default length of 9 codons.

## <a name="snpgenie-within"></a>SNPGenie Within-Group

The script **snpgenie\_within\_group.pl** estimates mean *d*<sub>N</sub> and *d*<sub>S</sub> within a group of aligned sequences in a single FASTA file. Designed for use with sequence data that outsizes what can be handled by other software platforms, this script invokes parallelism for codons and bootstrapping, and as a result has one dependency: the <a target="_blank" href="http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm">**Parallel::ForkManager**</a> module. If this isn't already installed, please use CPAN to install it before running the script. <a target="_blank" href="http://www.cpan.org/modules/INSTALL.html">Click here</a> to learn how.

Note that within-group mean *d*<sub>N</sub> and *d*<sub>S</sub> are mathematically equivalent to *π*<sub>N</sub> and *π*<sub>S</sub> when groups contain sequences sampled from a single population. 

Provide the script with the following arguments: 

* **--fasta\_file\_name**: the name or path of the (one) FASTA file containing aligned nucleotide sequences.
* **--gtf\_file\_name**: the name or path of the (one) [GTF file](#gtf) containing CDS information.
* **--num\_bootstraps**: OPTIONAL. the number of bootstrap replicates to be performed for calculating the standard error of *d*<sub>N</sub>-*d*<sub>S</sub>. DEFAULT: 0.
* **--procs\_per\_node**: OPTIONAL. The number of parallel processes to be invoked for processing codons and bootstraps (speeds up SNPGenie if high performance computing resources are available). DEFAULT: 1.

The format for using this script is:

        snpgenie_within_group.pl --fasta_file_name=<aligned_sequences>.fa --gtf_file_name=<coding_annotations>.gtf --num_bootstraps=<bootstraps> --procs_per_node=<procs_per_node>

For example, on a node with access to 64 processors in a high performance computing environment, you might submit a job including the following (**DON'T** run this from the head node!):

        snpgenie_within_group.pl --fasta_file_name=my_virus_genomes.fa --gtf_file_name=my_virus_genes.gtf --num_bootstraps=10000 --procs_per_node=64

Results files will be generated in the working directory. For a detailed explanation of the results, consult the [output descriptions](#output) above (however, note that not all files and columns generated for the within-pool analysis will be present for the within-group analysis).

To execute bootstrapping and/or sliding window analyses, run the processing script **snpgenie\_within\_group\_processor.pl** with the following arguments:

* **--codon\_file**: the name or path of the **within\_group\_codon\_results.txt** file produced using snpgenie\_within\_group.pl.
* **--sliding\_window\_size**: the size of the sliding windows to be calculated, in codons.
* **--sliding\_window\_step**: the number of codons to advance forward for each subsequent sliding window.
* **--num\_bootstraps**: OPTIONAL. the number of bootstrap replicates to be performed for calculating the standard error of *d*<sub>N</sub>-*d*<sub>S</sub> for genes and sliding windows. DEFAULT: 0.

The format for using the script is:

	snpgenie_within_group_processor.pl --codon_file=within_group_codon_results.txt --sliding_window_size=10 --sliding_window_step=1 --num_bootstraps=10000

A typical workflow thus includes one run of **snpgenie\_within\_group.pl** to obtain codon results, followed by one or more runs of **snpgenie\_within\_group\_processor.pl** to obtain sliding window results with bootstraps.

## <a name="snpgenie-between"></a>SNPGenie Between-Group

The script **snpgenie\_between\_group.pl** estimates mean *d*<sub>N</sub> and *d*<sub>S</sub> between two or more groups of sequences in FASTA format. Designed for use with sequence data that outsizes what can be handled by other software platforms, this script invokes parallelism for codons and bootstrapping, and as a result has one dependency: the <a target="_blank" href="http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm">**Parallel::ForkManager**</a> module. If this isn't already installed, please use CPAN to install it before running the script. <a target="_blank" href="http://www.cpan.org/modules/INSTALL.html">Click here</a> to learn how.

Just run the script in a directory containing **two or more FASTA** files, each containing a group of aligned sequences which are also aligned to one another across groups (files), ending with a .fa, .fas, or .fasta extension. In other words, each FASTA file contains one group. For example, one might wish to compare sequences of the West Nile Virus 2K gene derived from three different mosquito vector species. In such a case, all sequences would first be aligned. Next, all sequences from mosquito species 1 would be placed together in one FASTA file; all sequences from mosquito species 2 would be placed in a second FASTA file; and all sequences from mosquito species 3 would be placed in a third FASTA file.

SNPGenie automatically detects all FASTA files in the working directory. The GTF file (see [Gene Transfer Format](#gtf) for details) must be specified at the command line, along with an optional the number of bootstraps and processes:

	snpgenie_between_group.pl --gtf_file=<CDS_annotations>.gtf --num_bootstraps=<integer> --procs_per_node=<integer>

Results files will be generated in the working directory. For a detailed explanation of the results, consult the [output descriptions](#output) above (however, note that not all files and columns generated for the within-pool analysis will be present for the within-group analysis).

To execute bootstrapping and/or sliding window analyses, run the processing script **snpgenie\_between\_group\_processor.pl** with the following arguments:

* **--between\_group\_codon\_file**: the name or path of the **between\_group\_codon\_results.txt** file produced using snpgenie\_between\_group.pl.
* **--sliding\_window\_size**: the size of the sliding windows to be calculated, in codons.
* **--sliding\_window\_step**: the number of codons to advance forward for each subsequent sliding window.
* **--codon\_min\_taxa\_total**: the minimum number of total distinct sequences or taxa *across both groups* having defined (non-gap) codons at a codon position. If this criterion is not met for the groups being compared, the codon position is ignored.
* **--codon\_min\_taxa\_group**: the minimum number of distinct sequences or taxa *in a single groups* having defined (non-gap) codons at a codon position. If this criterion is not met for either of the groups being compared, the codon position is ignored.
* **--num\_bootstraps**: OPTIONAL. the number of bootstrap replicates to be performed for calculating the standard error of *d*<sub>N</sub>-*d*<sub>S</sub> for genes and sliding windows. DEFAULT: 0.
* **--procs\_per\_node**: OPTIONAL. The number of parallel processes to be invoked for processing codons and bootstraps (speeds up SNPGenie if high performance computing resources are available). DEFAULT: 1.

The format for using the script is:

	snpgenie_between_group_processor.pl --between_group_codon_file=<between_group_codon_results>.txt --sliding_window_size=10 --sliding_window_step=1 --codon_min_taxa_total=5 --codon_min_taxa_group=2 --num_bootstraps=10000 --procs_per_node=8 

A typical workflow thus includes one run of **snpgenie\_between\_group.pl** to obtain codon results, followed by one or more runs of **snpgenie\_between\_group\_processor.pl** to obtain sliding window and/or whole-gene results with various codon criteria, numbers of bootstraps, and sliding window sizes.

## <a name="how-snpgenie-works"></a>How SNPGenie Works

SNPGenie calculates gene and nucleotide diversities for sites in a protein-coding sequence. Nucleotide diversity, *π*, is the average number of nucleotide variants per nucleotide site for all pairwise comparisons. To distinguish between nonsynonymous and synonymous differences and sites, it is necessary to consider the codon context of each nucleotide in a sequence. This is why the user must submit the starting and ending sites of the coding regions in the .gtf file, along with the reference FASTA sequence file, so that the numbers of nonsynonymous and synonymous sites for each codon may be accurately estimated by a derivation of the Nei-Gojobori (1986) method. 

SNPGenie first splits the coding sequence into codons, each of which contains 3 sites. The software then determines which are nonsynonymous and synonymous by testing all possible polymorphisms present at each site of every codon in the sequence. Because different nucleotide variants at the same site may lead to both nonsynonymous and synonymous polymorphisms, fractional sites may occur (*e.g.*, only 2 of 3 possible nucleotide substitutions at the third position of AGA cause an amino acid change; thus, that site is considered 2/3 nonsynonymous and 1/3 synonymous). Next, the SNP report is consulted for the presence of variants to produce a revised estimate of the number of sites based on the presence of variant alleles. Variants are incorporated through averaging, weighted by their frequency. Although it is relatively rare, high levels of sequence variation may alter the number of nonsynonymous and synonymous sites in a particular codon. 

SNPGenie calculates the number of nucleotide differences for each codon in each ORF as follows. For every variant in the SNP Report, the number of variants is calculated as the product of the variant’s relative frequency and the coverage at that site. For each variant nucleotide (up to 3 non-reference nucleotides), the number of variants is stored, and their sum is subtracted from the coverage to yield the reference’s absolute frequency. Next, for each pairwise nucleotide comparison at the site, it is determined whether the comparison represents a nonsynonymous or synonymous change. If the former, the product of their absolute frequencies contributes to the number of nonsynonymous pairwise differences; if the latter, it contributes to the number of synonymous pairwise differences. When comparing codons with more than one nucleotide difference, all possible mutational pathways are considered, per the method of Nei & Gojobori (1986). The sum of pairwise differences is divided by the total number of pairwise comparisons at the codon (<sub>*n*</sub>C<sub>2</sub>, where *n* = coverage) to yield the mean number of differences per site of each type. This is calculated separately for nonsynonymous and synonymous comparisons. Calculating nucleotide diversity codon-by-codon enables sliding-window analyses that may help to pinpoint important nucleotide regions subject to varying forms of natural selection.

For further background, see Nelson & Hughes (2015) in the [References](#references).

## <a name="additional-scripts"></a>Additional Scripts

Additional scripts, provided to help in the preparation of SNPGenie input, have moved to the <a target="_blank" href="https://github.com/chasewnelson/CHASeq">CHASeq</a> repository. There you can find gb2gtf.pl, gff2gtf.pl, split_fasta.pl, vcf2revcom.pl, and lots of other useful stuff. Enjoy!

## <a name="troubleshooting"></a>Troubleshooting

* **"coverage is 0" WARNING.** SNPGenie will run even when some sites lack sequencing coverage. Sometimes this presents no problems, *e.g.*, if all sequences lack coverage in the same regions or at terminal sites. **However**, if coverage differs between samples, this will affect diversity estimates and will need to be accounted for: one sample might have relatively low diversity simply because numerous sites lack coverage. Such a difference is not biological, but rather a result of incomplete data. As always, it is up to the user to understand their own data, and to only supply samples that are biologically meaningful.

* Using **Windows**? Unfortunately, SNPGenie was written for Unix systems, including Mac, which have Perl installed by default. Windows sometimes doesn't, but getting Perl installed is as simple as following these <a target="_blank" href="http://learn.perl.org/installing/windows.html">three-minute download instructions</a>, and with any luck the program will work! Just open the Windows Command Prompt, and remember to type "perl" first when you run SNPGenie, *i.e.*, type "perl snpgenie.pl". Unfortunately, if this does not work, we are unable to provide any further Windows support.

* Try preceding the whole command line with "perl" to make sure SNPGenie is being treated as a script. For example:

        perl snpgenie.pl --minfreq=0.01 --snpreport=mySNPreport.txt --fastafile=myFASTA.fa --gtffile=myGTF.gtf
    
* Make sure that the first line of the script contains the correct path to your machine's copy of Perl. Right now SNPGenie assumes it's located at **/usr/bin/perl**. You can find out if this matches your machine by typing **which perl** at the command line. If it's a different path than the one above, copy and paste your machine's path after the **#!** at the beginning of the SNPGenie script.

* Make sure that the script is executable at the command line, as follows:

        chmod +x snpgenie.pl	

* Make sure end-of-line characters are in Unix LF (\n) format. Although an attempt has been made to ensure SNPGenie can accept Windows CRLF (\r\n) or Mac CR (\r) formats, these can sometimes introduce unexpected problems causing SNPGenie to crash or return all 0 values. This is often a problem when saving Excel spreadsheets on Mac, which defaults to a Mac newline format. Trying changing the newline character to Unix LF using a free program such as <a target="_blank" href="http://www.barebones.com/products/textwrangler/">TextWrangler</a>.

* Make sure CLC files are tab (\t)-delimited.

* Make sure Geneious files are comma-separated.

* Make sure the SNP reports contain all necessary columns. (See sections on Input above.)

* Make sure the Frequency (CLC) or Variant Frequency (Geneious) SNP report columns contain a percentage, not a decimal (e.g., 11.0% rather than 0.11).

* Make sure the SNP calling frame correct (e.g., do the codons for a product in the SNP report begin with ATG and end with TAA, TAG, or TGA?).

* Make sure the product coordinates in the GTF file are correct. (You might use a free program such as <a target="_blank" href="http://www.megasoftware.net/">MEGA</a> to check that the CDS coordinates begin with ATG and end with TAA, TAG, or TGA, if that's what you expected.)

## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>Nelson CW, Moncla LH, Hughes AL (2015) <a target="_blank" href="http://bioinformatics.oxfordjournals.org/content/31/22/3709.long">SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data</a>. *Bioinformatics* **31**(22):3709-11. doi: 10.1093/bioinformatics/btv449.

and this page:

>https://github.com/chasewnelson/snpgenie

## <a name="studies-using-snpgenie"></a>Studies Using SNPGenie

* Schmitt K, *et al.* (2017) <a target="_blank" href="https://www.sciencedirect.com/science/article/pii/S0042682217302167">Modeling the evolution of SIV sooty mangabey progenitor virus towards HIV-2 using humanized mice</a>. *Virology* **510**: 175-84.
* Huang SW, *et al.* (2017) <a target="_blank" href="http://jvi.asm.org/content/91/23/e01062-17.long">A selective bottleneck shapes the evolutionary mutant spectra of enterovirus A71 during viral dissemination in humans</a>. *Journal of Virology* **91**(23): e01062-17.
* Wilson JR, *et al.* (2017) <a target="_blank" href="https://www.sciencedirect.com/science/article/pii/S0042682217302799?via%3Dihub">An influenza A virus (H7N9) anti-neuraminidase monoclonal antibody protects mice from morbidity without interfering with the development of protective immunity to subsequent homologous challenge</a>. *Virology* **511**: 214-21.
* Mirabello L, *et al.* (2017) <a target="_blank" href="http://www.cell.com/cell/abstract/S0092-8674(17)30889-9">HPV16 E7 genetic conservation is critical to carcinogenesis</a>. *Cell* **170**(6): 1164-74.
* Blanc-Mathieu R, *et al.* (2017). <a target="_blank" href="http://advances.sciencemag.org/lookup/doi/10.1126/sciadv.1700239">Population genomics of picophytoplankton unveils novel chromosome hypervariability</a>. *Science Advances* **3**(7):e1700239.
* Kutnjak D, Elena SF, and Ravnikar M (2017).  <a target="_blank" href="http://jvi.asm.org/content/91/16/e00690-17.long">Time-Sampled Population Sequencing Reveals the Interplay of Selection and Genetic Drift in Experimental Evolution of *Potato Virus Y*</a>. *Journal of Virology* **91**(16):e00690-17.
* Stam R, Scheikl D, and Tellier A (2016). <a target="_blank" href="https://academic.oup.com/gbe/article-lookup/doi/10.1093/gbe/evw094">Pooled Enrichment Sequencing Identifies Diversity and Evolutionary Pressures at NLR Resistance Genes within a Wild Tomato Population</a>. *Genome Biology and Evolution* **8**(5):1501–1515.
* Bailey AL, *et al.* (2016) <a target="_blank" href="http://jvi.asm.org/content/90/15/6724.abstract">Arteriviruses, pegiviruses, and lentiviruses are common among wild African monkeys</a>. *Journal of Virology* **90**(15):6724-6737.
* Moncla LH, *et al.* (2016) <a target="_blank" href="http://www.cell.com/cell-host-microbe/abstract/S1931-3128(16)30010-5">Selective bottlenecks shape evolutionary pathways taken during mammalian adaptation of a 1918-like avian influenza virus</a>. *Cell Host & Microbe* **19**(2):169-80.
* Gellerup D, *et al.* (2015) <a target="_blank" href="http://jvi.asm.org/content/early/2015/10/16/JVI.02587-15.abstract">Conditional immune escape during chronic SIV infection</a>. *Journal of Virology* **90**(1):545-52.
* Nelson CW and Hughes AL (2015) <a target="_blank" href="http://www.sciencedirect.com/science/article/pii/S1567134814004468">Within-host nucleotide diversity of virus populations: Insights from next-generation sequencing</a>. *Infection, Genetics and Evolution* **30**:1-7.
* Bailey AL, *et al.* (2014) <a target="_blank" href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090714">High genetic diversity and adaptive potential of two simian hemorrhagic fever viruses in a wild primate population</a>. *PLoS ONE* **9**(3):e90714.
* Wilker PR, *et al.* (2013) <a target="_blank" href="http://www.nature.com/ncomms/2013/131023/ncomms3636/abs/ncomms3636.html?message-global=remove">Selection on haemagglutinin imposes a bottleneck during mammalian transmission of reassortant H5N1 influenza viruses</a>. *Nature Communications* **4**:2636.

## <a name="contact"></a>Contact
If you have questions about using SNPGenie, please click on the "<a target="_blank" href="https://github.com/chasewnelson/snpgenie/issues">Issues</a>" tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Correspondance should be addressed to Chase W. Nelson at the American Museum of Natural History: 

* cnelson <**AT**> amnh <**DOT**> org

## <a name="references"></a>References

* Knapp EW, Irausquin SJ, Friedman R, Hughes AL (2011) <a target="_blank" href="http://link.springer.com/article/10.1007%2Fs12686-010-9372-5">PolyAna: analyzing synonymous and nonsynonymous polymorphic sites</a>. *Conserv Genet Resour* **3**:429-431.
* Lynch M, Bost D, Wilson S, Maruki T, Harrison S (2014) <a target="_blank" href="http://gbe.oxfordjournals.org/content/early/2014/04/30/gbe.evu085">Population-genetic inference from pooled-sequencing data</a>. *Genome Biol. Evol.* **6**(5):1210-1218.
* Nei M, Gojobori T (1986) <a target="_blank" href="http://mbe.oxfordjournals.org/content/3/5/418.short?rss=1&ssource=mfc">Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions</a>. *Mol Biol Evol* **3**:418-26.
* Nelson CW, Hughes AL (2015) <a target="_blank" href="http://www.sciencedirect.com/science/article/pii/S1567134814004468">Within-host nucleotide diversity of virus populations: Insights from next-generation sequencing</a>. *Infection, Genetics and Evolution* **30**:1-7.

Image Copyright 2015 Elizabeth Ogle
