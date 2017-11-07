# CHASeq
Perl scripts providing Computational Help for the Analysis of SEQuence data

## <a name="before-you-begin"></a>Before You Begin

I provide these scripts for anyone to use freely. Please just cite this GitHub page (below). They are meant for use on Unix/Mac machines; no support is offered for Windows.

Most scripts are designed to take unnamed arguments for simple data file manipulations, in the following format:

        CHASeq_script.pl <argument1> <argument2>

where \<argument1\> is ommitted and replaced with the desired input value. Some more complicated scripts will contain named arguments in Perl's long form, in the following format:

        CHASeq_script.pl --argument-name=<value>
        
where \<value\> is ommitted and replaced with the desired input value.

### Troubleshooting

If a script isn't working, try working through the following checklist:

* Remember that you must open the Terminal and navigate to the directory which contains the script before you're able to call it. More advanced users may wish to place the script(s) in a directory included in your computer's PATH so that it can be called from any directory.
* Did you place the script in the same directory as your data? If not, you could move the script to the data-containing directory, or *vice versa*. Another option is to provide the full path of the input files to the script, e.g., **/Users/cwnelson88/Desktop/my\_project/my\_data.fasta**
* These scripts assume that your computer's copy of Perl is located at **/usr/bin/perl**. You can check whether this is the case by opening the Terminal and typing **which perl** at the command line. If this does not return **/usr/bin/perl**, then (1) copy the path it provides; (2) open the CHASeq script; (3) replace **#! /usr/bin/perl** at the top of the script with your own computer's path.
* Scripts must be made executable before use. You can check whether this is the case my opening the Terminal, navigating to the directory containing the script, and typing **ls -l**, which will return something like the following:

        -rwxr-xr-x@  1  name  staff  155  Sep 25  2013  script.pl

	If the string of letters at the beginning of the line containing your script (here, it's **-rwxr-xr-x**) *does not* contain an 'x', it means it is not yet executable. You can add executable status by typing **chmod +x \<script.pl\>**, where \<script.pl\> is replaced with the script's name.

## <a name="citation"></a>Citation

When using this software, please refer to an cite:

>CHASeq software package, https://github.com/chasewnelson/CHASeq

## <a name="script-descriptions"></a>Script Descriptions

* **gb2gtf.pl**. At the command line, provide this script with one argument: a GenBank (.gbk) file. It will extract the coding element (CDS) annotations to produce a Gene Transfer Format (.gtf) file ready for SNPGenie. Not working? Let us know, and we'll improve it! Here's an example:

        gb2gtf.pl my_genbank_file.gbk
        
* **gff2gtf.pl**. At the command line, provide this script with one argument: a General Feature Format (.gff) file. It will extract the coding element annotations to produce a Gene Transfer Format (.gtf) file ready for SNPGenie, with "gene_id" annotations identified using the GFF "ID" tag. Not working, or need a different tag? Let us know, and we'll improve it! Here's an example:

        gff2gtf.pl my_gff_file.gff

* **split_fasta.pl**. At the command line, provide this script with one argument: a FASTA (.fa or .fasta) file containing multiple sequences. This script will create multiple files in the working directory, each containing one of the sequences. Here's an example:

        split_fasta.pl my_multi_fasta_file.fasta

* **vcf2revcom.pl**. This script automates the creation of the reverse complement input files. At the command line, provide this script with three arguments, in the following order: 
	1. A '+' strand FASTA (.fa or .fasta) file containing the reference sequence against which SNPs were called;
	2. A '+' strand GTF file containing both '+' and '–' strand products from the '+' strand point of view; and 
	3. A '+' strand SNP report in VCF format.

	This script will then create a '-' strand (reverse complement) version of each file in the working directory, with "_revcom" concatenated to the original file name. Here's an example:

        vcf2revcom.pl my_snp_report.vcf my_reference_sequence.fasta my_cds_file.gtf