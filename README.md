# CHASeq
Perl scripts providing **C**omputational **H**elp for the **A**nalysis of **Seq**uence data

## <a name="contents"></a>Contents
* **[Before You Begin](#before-you-begin)**
* **[Citation](#citation)**
* **[Troubleshooting](#troubleshooting)**
* **[Scripts](#scripts)**
	* **[calculate\_p\_distance.pl](#calculate-p-distance)**. You want to calculate a *p*-distance between two nucleotide sequences.
	* **[extract\_codon\_from\_ANN\_VCFs.pl](#extract-codon-from-ANN-VCFs)**. You want to pull codon variants out of annotated VCF files (i.e., files that have been annotated using the <a target="_blank" href="http://snpeff.sourceforge.net/SnpSift.html">snpeff -formatEff</a> option), so that reference and variant codons can be compared. 
	* **[extract\_fasta\_by\_sites.pl](#extract-fasta-by-sites)**. You want to create separate FASTA files for segments (e.g., each of the genes) in an aligned sequence file.
	* **[gb2gtf.pl](#gb-to-gtf)**. You want to create a <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file from a GenBank file, compatible with SNPGenie input.
	* **[get\_random\_integers.pl](#get-random-integers)**. You want to get a certain number of random integer values in a range.
	* **[gff2gtf.pl](#gff-to-gtf)**. You want to convert a GFF file to a simpler <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, compatible with SNPGenie input.
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

* <a name="calculate-p-distance"></a>**calculate\_p\_distance.pl**. You want to calculate a *p*-distance between two nucleotide sequences. At the command line, provide this script with two arguments: two FASTA (.fa or .fasta) files, each containing one sequence, which are aligned to each other. This script will exclude positions which are gaps (-) or undetermined (N) in both sequences and return a *p*-distance. Here's an example:

        calculate_p_distance.pl <aligned_seq_1.fasta> <aligned_seq_2.fasta>

* <a name="extract-codon-from-ANN-VCFs"></a>**extract\_codon\_from\_ANN\_VCFs.pl**. You want to pull codon variants out of annotated VCF files (i.e., files that have been annotated using the <a target="_blank" href="http://snpeff.sourceforge.net/SnpSift.html">snpeff -formatEff</a> option), so that reference and variant codons can be compared. At the command line, call this script in a directory containing one or more files ending in **.ann.vcf**. It will extract the information in the "CODON: Codon change" field (e.g. **ggT/ggG**) for each record and output to the files ref\_codons\_seq.txt, alt\_codons\_seq.txt, and codon\_metadata.txt. Here's an example:

        extract_codon_from_ANN_VCFs.pl

* <a name="extract-fasta-by-sites"></a>**extract\_fasta\_by\_sites.pl**. You want to create separate FASTA files for segments (e.g., each of the genes) in an aligned sequence file. At the command line, provide this script with two arguments: (1) one FASTA file containing one or more aligned sequences; and (2) a <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file containing the CDS products to extract. It will extract the sites for each coding element (CDS) annotation in the Gene Transfer Format, and output one FASTA file for each. Here's an example:

        extract_fasta_by_sites.pl <multiple_aligned_seqs.fasta> <gene_coordinates_to_extract.gtf>

* <a name="gb-to-gtf"></a>**gb2gtf.pl**. (*Helpful for preparing **<a target="_blank" href="https://github.com/chasewnelson/snpgenie">SNPGenie</a>** input!*) You want to create a <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file from a GenBank file, compatible with SNPGenie input. At the command line, provide this script with one argument: a GenBank (.gbk) file. It will extract the coding element (CDS) annotations to produce a Gene Transfer Format (.gtf) file ready for SNPGenie. Not working? Let us know, and we'll improve it! Here's an example:

        gb2gtf.pl <my_genbank_file.gbk>
        
* <a name="get-random-integers"></a>**get\_random\_integers.pl**. You want to get a certain number of random integer values in a range. At the command line, provide this script with three arguments: (1) number of random integer values wanted; (2) bottom of range (inclusive); (3) top of range (inclusive). Returns the specified number of random integers in the desired range. Here's an example to return 10 values in the range [3,999]:

        get_random_integers.pl 10 3 999
        
* <a name="gff-to-gtf"></a>**gff2gtf.pl**. (*Helpful for preparing **<a target="_blank" href="https://github.com/chasewnelson/snpgenie">SNPGenie</a>** input!*) You want to convert a GFF file to a simpler <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, compatible with SNPGenie input. At the command line, provide this script with one argument: a General Feature Format (.gff) file. It will extract the coding element annotations to produce a Gene Transfer Format (.gtf) file ready for SNPGenie, with "gene_id" annotations identified using the GFF "ID" tag. Not working, or need a different tag? Let us know, and we'll improve it! Here's an example:

        gff2gtf.pl <my_gff_file.gff>

* <a name="split-fasta"></a>**split\_fasta.pl**. You want to create an individual FASTA file for each sequence in another FASTA file. At the command line, provide this script with one argument: a FASTA (.fa or .fasta) file containing multiple sequences. This script will create multiple files in the working directory, each containing one of the sequences. Here's an example:

        split_fasta.pl <my_multi_fasta_file.fasta>

* <a name="store-fasta-by-ID"></a>**store\_fasta\_by\_ID.pl**. You want to create a new FASTA file containing only a certain subset of another FASTA. At the command line, provide this script with two arguments: (1) a FASTA (.fa or .fasta) file containing multiple sequences; and (2) a .txt file containing one column with the FASTA IDs you want in your new FASTA file — there should be NO HEADER. This script will create a new FASTA file in the working directory, containing the sequences whose headers begin with the IDs in argument (2). Redirect (>) the output to create a file. Here's an example:

        store_fasta_by_ID.pl <all_seqs.fasta> <wanted_seqs_headers.txt> > <just_wanted_seqs.fasta>

* <a name="vcf-to-revcom"></a>**vcf2revcom.pl**. (*Helpful for preparing **<a target="_blank" href="https://github.com/chasewnelson/snpgenie">SNPGenie</a>** input!*) You want to convert a VCF SNP report file, along with its accompanying FASTA file and <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file, to the reverse complement strand. This script automates the creation of such reverse complement files, compatible with SNPGenie input. Note that the resulting SNP report will not be a VCF file, but rather a CLC Genomics Workbench format file. At the command line, provide this script with three arguments, in the following order: 
	1. A '+' strand FASTA (.fa or .fasta) file containing the reference sequence against which SNPs were called;
	2. A '+' strand <a target="_blank" href="http://mblab.wustl.edu/GTF22.html">GTF</a> file containing both '+' and '–' strand products from the '+' strand point of view; and 
	3. A '+' strand SNP report in VCF format.

	This script will then create a '-' strand (reverse complement) version of each file in the working directory, with "_revcom" concatenated to the original file name. Here's an example:

        vcf2revcom.pl <my_snp_report.vcf> <my_reference_sequence.fasta> <my_cds_file.gtf>