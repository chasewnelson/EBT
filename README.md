# CHASeq
Perl scripts providing **C**omputational **H**elp for the **A**nalysis of **Seq**uence data.

## <a name="contents"></a>Contents
* **[Before You Begin](#before-you-begin)**
* **[Citation](#citation)**
* **[Troubleshooting](#troubleshooting)**
* **[Scripts](#scripts)**
	* **[aligned\_fasta\_group\_diffs.pl](#aligned-fasta-group-diffs)**. You have two or more groups of sequences, aligned to one another but in separate FASTA files, and want to identify the sites at which the groups exhibit differences.
	* **[aligned\_fasta2site\_nt\_freqs.pl](#aligned-fasta-2-site-nt-freqs)**. You want to tabulate the number (and proportion) of each nucleotide at each variable position across an aligned nucleotide FASTA file.
	* **[calculate\_p\_distance.pl](#calculate-p-distance)**. You want to calculate a *p*-distance between two nucleotide sequences.
	* **[determine\_consensus\_IUPAC\_seqs.pl](#determine-consensus-IUPAC-seqs)**. You want to determine the consensus and/or IUPAC sequence(s) for an aligned nucleotide FASTA file. 
	* **[extract\_codon\_from\_ANN\_VCFs.pl](#extract-codon-from-ANN-VCFs)**. You want to pull codon variants out of annotated VCF files (i.e., files that have been annotated using the <a target="_blank" href="http://snpeff.sourceforge.net/SnpSift.html">snpeff -formatEff</a> option), so that reference and variant codons can be compared. 
	* **[extract\_fasta\_by\_sites.pl](#extract-fasta-by-sites)**. You want to create separate FASTA files for segments (e.g., each of the genes) in an aligned sequence file.
	* **[gb2gtf.pl](#gb-to-gtf)**. You want to create a <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file from a GenBank file, compatible with SNPGenie input.
	* **[get\_random\_integers.pl](#get-random-integers)**. You want to get a certain number of random integer values in a range.
	* **[get\_unique\_values.pl](#get-unique-values)**. You have a table with a column containing multiple rows with the same value(s) (i.e., duplicates), and you want to extract just the unique ones (i.e., remove duplicate values) and generate summary statistics for the duplicates.
	* **[gff2gtf.pl](#gff-to-gtf)**. You want to convert a GFF file to a simpler <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, compatible with SNPGenie input.
	* **[LinkGe\_all\_site\_pairs](#LinkGe-all-site-pairs)**. You have an aligned BAM file and a list of sites with their SNPs, and want to use Gabriel Starrett's <a target="_blank" href="https://github.com/gstarrett/LinkGe">LinkGe</a> program to check for linkage of SNPs in the same reads.
	* **[remove\_seqs\_with\_stops.pl](#remove-seqs-with-stops)**. You have a FASTA file where each sequence is an in-frame coding sequence, and want to remove all sequences containing mid-sequence STOP codons.
	* **[split\_fasta.pl](#split-fasta)**. You have a FASTA file, and want to create an individual FASTA file for each sequence inside.
	* **[store\_fasta\_by\_ID.pl](#store-fasta-by-ID)**. You want to create a new FASTA file containing only a certain subset of another FASTA.
	* **[vcf2revcom.pl](#vcf-to-revcom)**. You want to convert a VCF SNP report file, along with its accompanying FASTA file and <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, to the reverse complement strand.

## <a name="before-you-begin"></a>Before You Begin

I provide these scripts for anyone to use freely for help in performing routine bioinformatics tasks with nucleotide sequence and related data. Please just [cite this GitHub page](#citation). The scripts are meant for use on Unix/Mac machines; no support is offered for Windows. Readers of code, be warned: very little effort has been made to 'clean up' my coding comments and alternative operations. Think of them as pseudogenes!

### How to Use

Most scripts are designed to take (unnamed) arguments for simple data file manipulations, in the following format:

        CHASeq_script.pl <argument1> <argument2>

where \<argument1\> is ommitted and replaced with the desired input value. Some more complicated scripts will contain named arguments in Perl's long form, in the following format:

        CHASeq_script.pl --argument-name=<value>
        
where \<value\> is ommitted and replaced with the desired input value.

## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>CHASeq software package, https://github.com/chasewnelson/CHASeq

## <a name="troubleshooting"></a>Troubleshooting

If a script isn't working, try working through the following checklist:

* Remember that you must open the Terminal and navigate to the directory which contains the script before you're able to call it. More advanced users may wish to place the script(s) in a directory included in your computer's PATH so that it can be called from any directory.
* Did you place the script in the same directory as your data? If not, you could move the script to the data-containing directory, or *vice versa*. Another option is to provide the full path of the input files to the script, e.g., **/Users/cwnelson88/Desktop/my\_project/my\_data.fasta**
* These scripts assume that your computer's copy of Perl is located at **/usr/bin/perl**. You can check whether this is the case by opening the Terminal and typing **which perl** at the command line. If this does not return **/usr/bin/perl**, then (1) copy the path it provides; (2) open the CHASeq script; (3) replace **#! /usr/bin/perl** at the top of the script with your own computer's path.
* Scripts must be made executable before use. You can check whether this is the case my opening the Terminal, navigating to the directory containing the script, and typing **ls -l**, which will return something like the following:

        -rwxr-xr-x@  1  name  staff  155  Sep 25  2013  script.pl

	If the string of letters at the beginning of the line containing your script (here, it's **-rwxr-xr-x**) *does not* contain an 'x', it means it is not yet executable. You can add executable status by typing **chmod +x \<script.pl\>**, where \<script.pl\> is replaced with the script's name.

## <a name="scripts"></a>Scripts

* <a name="aligned-fasta-group-diffs"></a>**aligned\_fasta\_group\_diffs.pl**. You have two or more groups of sequences, aligned to one another but in separate FASTA files, and want to identify the sites at which the groups exhibit differences. First, make sure all sequences are aligned to one another (even across groups), place each group of sequences in a separate FASTA file, and create a directory containing all the group FASTA files (but no others). Then, at the command line, call this script. A results file will be placed in the working directory describing positions at which the groups differ in their major nucleotide. For more control use the following options:
	* *--min\_variant\_maj\_nt_freq* to specify a minimum frequency cutoff for assigning a group's majority (consensus) nucleotide. Must be a decimal. Default=0.9 (90%).
	* *--min\_site\_coverage* to specify a minimum number of defined sequences required for calculating a group's frequencies. Must be an integer. For example, some sites in some groups may be mostly gaps (-) or undetermined (N), in which case we may not want to consider frequency calculations reliable. Default=5.

	Here is an example using defaults:

        aligned_fasta_group_diffs.pl

	Here is an example using using a frequency cutoff of 90% with a minimum of 8 defined nucleotides per site for each group:

        aligned_fasta_group_diffs.pl --min_variant_maj_nt_freq=.9 --min_site_coverage=8

* <a name="aligned-fasta-2-site-nt-freqs"></a>**aligned\_fasta2site\_nt\_freqs.pl**. You want to tabulate the number (and proportion) of each nucleotide at each variable position across an aligned nucleotide FASTA file. At the command line, call this script with one argument, an aligned FASTA file. One TAB-delimited results file will be placed in the working directory (\*\_site_summary.txt) and brief summary statistics will be printed to the Terminal. Here's an example:

        aligned_fasta2site_nt_freqs.pl <aligned_seqs.fasta>

* <a name="calculate-p-distance"></a>**calculate\_p\_distance.pl**. You want to calculate a *p*-distance between two nucleotide sequences. At the command line, provide this script with two arguments: two FASTA (.fa or .fasta) files, each containing one sequence, which are aligned to each other. This script will exclude positions which are gaps (-) or undetermined (N) in both sequences and return a *p*-distance. Here's an example:

        calculate_p_distance.pl <aligned_seq_1.fasta> <aligned_seq_2.fasta>

* <a name="determine-consensus-IUPAC-seqs"></a>**determine\_consensus\_IUPAC\_seqs.pl**. You want to determine the consensus and/or IUPAC sequence(s) for an aligned nucleotide FASTA file. At the command line, call this script with one argument, an aligned FASTA file. Two results files will be placed in the working directory, one with a consensus sequence (\*\_consensus.fa) and one with an IUPAC sequence (\*\_IUPAC.fa). Brief summary statistics will be printed to the Terminal. Here's an example:

        determine_consensus_IUPAC_seqs.pl <aligned_seqs.fasta>

* <a name="extract-codon-from-ANN-VCFs"></a>**extract\_codon\_from\_ANN\_VCFs.pl**. You want to pull codon variants out of annotated VCF files (i.e., files that have been annotated using the <a target="_blank" href="http://snpeff.sourceforge.net/SnpSift.html">snpeff -formatEff</a> option), so that reference and variant codons can be compared. At the command line, call this script in a directory containing one or more files ending in **.ann.vcf**. It will extract the information in the "CODON: Codon change" field (e.g. **ggT/ggG**) for each record and output to the files ref\_codons\_seq.txt, alt\_codons\_seq.txt, and codon\_metadata.txt. Here's an example:

        extract_codon_from_ANN_VCFs.pl

* <a name="extract-fasta-by-sites"></a>**extract\_fasta\_by\_sites.pl**. You want to create separate FASTA files for segments (e.g., each of the genes) in an aligned sequence file. At the command line, provide this script with two arguments: (1) one FASTA file containing one or more aligned sequences; and (2) a <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file containing the CDS products to extract. It will extract the sites for each coding element (CDS) annotation in the Gene Transfer Format, and output one FASTA file for each. Here's an example:

        extract_fasta_by_sites.pl <multiple_aligned_seqs.fasta> <gene_coordinates_to_extract.gtf>

* <a name="gb-to-gtf"></a>**gb2gtf.pl**. (*Helpful for preparing **<a target="_blank" href="https://github.com/chasewnelson/snpgenie">SNPGenie</a>** input!*) You want to create a <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file from a GenBank file, compatible with SNPGenie input. At the command line, provide this script with one argument: a GenBank (.gbk) file. It will extract the coding element (CDS) annotations to produce a Gene Transfer Format (.gtf) file ready for SNPGenie. Not working? Let us know, and we'll improve it! Here's an example:

        gb2gtf.pl <my_genbank_file.gbk>
        
* <a name="get-random-integers"></a>**get\_random\_integers.pl**. You want to get a certain number of random integer values in a range. At the command line, provide this script with three arguments: (1) number of random integer values wanted; (2) bottom of range (inclusive); (3) top of range (inclusive). Returns the specified number of random integers in the desired range. Here's an example to return 10 values in the range [3,999]:

        get_random_integers.pl 10 3 999

* <a name="get-unique-values"></a>**get\_unique\_values.pl**. You have a table with a column containing multiple rows with the same value(s) (i.e., duplicates), and you want to extract just the unique ones (i.e., remove duplicate values) and generate summary statistics for the duplicates. At the command line, provide this script with two arguments: (1) the name of the input tab-delimited file; and (2) the column number to analyze. The script outputs the unique values in the column, i.e., duplicate values eliminated, including the header. Also creates a summary statistics file in the working directory, by the name *\_unique\_values\_summary.txt. Here's an example:

        get_unique_values.pl <my_table.txt> <column_number>

* <a name="gff-to-gtf"></a>**gff2gtf.pl**. (*Helpful for preparing **<a target="_blank" href="https://github.com/chasewnelson/snpgenie">SNPGenie</a>** input!*) You want to convert a GFF file to a simpler <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, compatible with SNPGenie input. At the command line, provide this script with one argument: a General Feature Format (.gff) file. It will extract the coding element annotations to produce a Gene Transfer Format (.gtf) file ready for SNPGenie, with "gene_id" annotations identified using the GFF "ID" tag. Not working, or need a different tag? Let us know, and we'll improve it! Here's an example:

        gff2gtf.pl <my_gff_file.gff>

* <a name="remove-seqs-with-stops"></a>**remove\_seqs\_with\_stops.pl**. You have a FASTA file where each sequence is an in-frame coding sequence, and want to remove all sequences containing mid-sequence STOP codons. At the command line, provide this script with one argument: a FASTA (.fa or .fasta) file containing multiple coding sequences that all begin at the first position of the first codon (they need not be aligned to each other). This script will create a new FASTA file in the working directory, containing just those sequences lacking an in-frame mid-sequence STOP codon. Here's an example:

        remove_seqs_with_stops.pl <my_fasta_file.fasta>

* <a name="LinkGe-all-site-pairs"></a>**LinkGe\_all\_site\_pairs.pl**. You have an aligned BAM file and a list of sites with their SNPs, and want to use Gabriel Starrett's <a target="_blank" href="https://github.com/gstarrett/LinkGe">LinkGe</a> program to check for linkage of SNPs in the same reads. For example, you may wish to test if high-frequency variant nucleotides are present in linked haplotypes in a deep-sequenced viral sample from a single host. At the command line, provide this script with the named arguments:
	* *--snp\_file*: REQUIRED. One CLC-style "SNP report" file with, at minimum, the following three columns (with headers):
		* **"Reference Position"**: the site number of the single nucleotide variant within the aligned BAM reads.
		* **"Reference"**: the reference (majority) nucleotide.
		* **"Allele"**: the allele (variant or mutant) nucleotide.
	* *--BAM\_file*: REQUIRED. Standard format; must be aligned; reads of any length.
	* *--max\_dist*: the maximum distance at which to search for pairs of sites covered by the same read. DEFAULT: 1000.

	This script will then determine all pairs of sites provided in the **snp\_file** which fall within **max\_dist** of one another; run LinkGe for all identified pairs; and output read linkage information for all individual pairs (*LinkGe_\** file generated in working directory) and total summary data for linkage of reference and variant nucleotides (printed to screen). Here's an example:
	
		LinkGe_all_site_pairs.pl --snp_file=<my_SNPs>.txt --BAM_file=<aligned_BAM>.bam --max_dist=500

	Note that this script requires you to <a target="_blank" href="https://github.com/gstarrett/LinkGe">install LinkGe</a>, and to <a target="_blank" href="http://osxdaily.com/2014/08/14/add-new-path-to-path-command-line/">add it to your computer's PATH</a>. Three dependencies are also required: Bio::DB::Sam, BioPerl, and Samtools version 0.1.15. See the <a target="_blank" href="https://github.com/gstarrett/LinkGe">LinkGe page</a> for details.

* <a name="split-fasta"></a>**split\_fasta.pl**. You want to create an individual FASTA file for each sequence in another FASTA file. At the command line, provide this script with one argument: a FASTA (.fa or .fasta) file containing multiple sequences. This script will create multiple files in the working directory, each containing one of the sequences. Here's an example:

        split_fasta.pl <my_fasta_file.fasta>

* <a name="store-fasta-by-ID"></a>**store\_fasta\_by\_ID.pl**. You want to create a new FASTA file containing only a certain subset of another FASTA. At the command line, provide this script with two arguments: (1) a FASTA (.fa or .fasta) file containing multiple sequences; and (2) a .txt file containing one column with the FASTA IDs you want in your new FASTA file — there should be NO HEADER. This script will create a new FASTA file in the working directory, containing the sequences whose headers begin with the IDs in argument (2). Redirect (>) the output to create a file. Here's an example:

        store_fasta_by_ID.pl <all_seqs.fasta> <wanted_seqs_headers.txt> > <just_wanted_seqs.fasta>

* <a name="vcf-to-revcom"></a>**vcf2revcom.pl**. (*Helpful for preparing **<a target="_blank" href="https://github.com/chasewnelson/snpgenie">SNPGenie</a>** input!*) You want to convert a VCF SNP report file, along with its accompanying FASTA file and <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, to the reverse complement strand. This script automates the creation of such reverse complement files, compatible with SNPGenie input. Note that the resulting SNP report will not be a VCF file, but rather a CLC Genomics Workbench format file. At the command line, provide this script with three arguments, in the following order: 
	1. A '+' strand FASTA (.fa or .fasta) file containing the reference sequence (ALL UPPERCASE) against which SNPs were called;
	2. A '+' strand <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file containing both '+' and '–' strand products from the '+' strand point of view; and 
	3. A '+' strand SNP report in VCF format.

	This script will then create a '-' strand (reverse complement) version of each file in the working directory, with "_revcom" concatenated to the original file name. Here's an example:

        vcf2revcom.pl <my_reference_sequence.fasta> <my_cds_file.gtf> <my_snp_report.vcf>