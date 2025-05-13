# Getting started with Bowtie 2: Lambda phage example
Bowtie 2 is available as a module on Explorer. The files you will be working with are the lambda phage genome and a set of paired-end reads that you downloaded for the previous tutorial. (If you do not have them, download them from the "examples" folder at [this GitHub site: ]([https://pages.github.com/)(https://github.com/BenLangmead/bowtie2)). In this exercise, adapted from the Bowtie2 manual, you will learn to index a genome and align the reads to it using a simple phage genome to start.

Use `head` to peek at each of the three files. Using UNIX commands, determine the length of the lambda phage genome, and the number of sequencing reads in each of the two fastq files. Show your commands and the output, with hash-comments, for the entire tutorial.

Before proceeding further, request a computing node and load the bowtie2 module in Explorer.

### Indexing a reference genome
The command to index a genome includes the index command (bowtie2-build) followed by the path to the genome sequence to be indexed followed by the prefix that will be placed at the start of all indexed files. The command as written below assumes I'm running the command one level below where the lambda genome file lives. You may need to change the path!

`bowtie2-build ../lambda_virus.fa lambda_virus`

The command should print many lines of output then quit. When the command completes, the current directory will contain six new files that all start with lambda_virus and end with .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. These files contain the indexed version of the lambda phage genome, in a compressed format. 

You can use bowtie2-build to create an index for a set of FASTA files obtained from any source, including sites such as UCSC, NCBI, and Ensembl. Later, you'll index a human chromosome, and you'll also work with an indexed version of the entire human genome!

### Aligning example reads
Stay in the directory that contains the lambda_virus index files. Next, run the command to align the reads to the indexed genome. `bowtie2` is the "align" command, "-x" indicates that the next component is the prefix of the indexed files, "-U" means the reads you're aligning are unpaired and appear immediately after, "-S" is the format of the output file, which you're calling "eg1.sam". Again, pay attention to where your files are and alter the path from the example as necessary.

`bowtie2 -x lambda_virus -U ../reads_1.fq -S eg1.sam`

The alignment results in SAM format are written to the file eg1.sam, and a short alignment summary is written to the console. Use the `head` command to see the first few lines of the SAM output, which will look something like this, but with different values:

![SAM file](https://github.com/jamiehenzy/Genome_analysis_SE_R/blob/main/SAM.png)

The first few lines (beginning with @) are SAM header lines, and the rest of the lines are SAM alignments, one line per read or mate. See the Bowtie 2 manual section on SAM output and the SAM specification for details about how to interpret the SAM file format.

### Paired-end example
To align paired-end reads included with Bowtie 2, stay in the same directory and run:

`bowtie2 -x lambda_virus -1 ../reads_1.fq -2 ../reads_2.fq -S eg2.sam`

This aligns a set of paired-end reads to the reference genome, with results written to the file eg2.sam. If you had used the filename eg1.sam here, your new results would have overwritten your previous results, so be careful!

### Using SAMtools/BCFtools downstream
SAMtools is a collection of tools for manipulating and analyzing SAM and BAM alignment files. BCFtools is a collection of tools for calling variants and manipulating VCF and BCF files, and it is typically distributed with SAMtools. Both of these programs are available as modules in Explorer. Using these tools together allows you to get from alignments in SAM format to variant calls in VCF format. 

Load the modules for samtools and bcftools into your environment (keeping bowtie2 active, as well) before proceeding.

Use samtools `view` to convert the SAM file for the paired-end alignment you performed into a BAM file. BAM is the binary format corresponding to the SAM text format. Run:

`samtools view -bS eg2.sam > eg2.bam`

Nothing will print to the console but you'll see a new file ending in .bam in your directory.

Use samtools `sort` to convert the BAM file to a sorted BAM file that can be used for calling variants.

`samtools sort eg2.bam -o eg2.sorted.bam`

We now have a sorted BAM file called eg2.sorted.bam. Sorted BAM is a useful format because the alignments are (a) compressed, which is convenient for long-term storage, and (b) sorted, which is convenient for variant discovery. To generate variant calls in a format known as VCF, run:

`bcftools mpileup -f ../lambda_virus.fa eg2.sorted.bam | bcftools view -Ov - > eg2.raw.bcf`

Then to view the variants, run the command below, directing output to a file that you can name whatever you want:

`bcftools view eg2.raw.bcf > variants.txt`

View the file!

## Submission
+ Directory and subdirectories as needed with all files related to the tutorial
+ User-friendly file and directory naming and arrangement
+ Saved session showing commands and output, hashed with brief statements of what the command will do
  
