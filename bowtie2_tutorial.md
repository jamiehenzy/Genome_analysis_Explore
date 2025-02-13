# Getting started with Bowtie 2: Lambda phage example
Bowtie 2 is available as a module on Explorer. The files you will be working with are the lambda phage genome, a set of paired-end reads, and a set of long reads available in our course data directory (if you have not previously downloaded them). In this exercise, adapted from the [Bowtie2 Manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), you will learn to index a genome and align the reads to it using a simple phage genome to start.

For this exercise, have your lambda sequences in one directory. In that same directory, make a new temporary directory. Remember to request a computing node, and load the bowtie2 module.

### Indexing a reference genome
Change into the temporary directory you made and run:

`bowtie2-build ../lambda_virus.fa lambda_virus`

The command should print many lines of output then quit. When the command completes, the current directory will contain six new files that all start with lambda_virus and end with .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. These files constitute the index - you're done!

You can use bowtie2-build to create an index for a set of FASTA files obtained from any source, including sites such as UCSC, NCBI, and Ensembl. When indexing multiple FASTA files, specify all the files using commas to separate file names. For more details on how to create an index with bowtie2-build, see the manual section on index building. You may also want to bypass this process by obtaining a pre-built index. See using a pre-built index below for an example.

### Aligning example reads
Stay in the directory created in the previous step, which now contains the lambda_virus index files. Next, run:

`bowtie2 -x lambda_virus -U ../reads_1.fq -S eg1.sam`

This runs the Bowtie 2 aligner, which aligns a set of unpaired reads to the Lambda phage reference genome using the index generated in the previous step. The alignment results in SAM format are written to the file eg1.sam, and a short alignment summary is written to the console. (Actually, the summary is written to the "standard error" or "stderr" filehandle, which is typically printed to the console.)

To see the first few lines of the SAM output, run:

`head eg1.sam`

You will see something like this:

![SAM file](https://github.com/jamiehenzy/Genome_analysis_SE_R/blob/main/SAM.png)

The first few lines (beginning with @) are SAM header lines, and the rest of the lines are SAM alignments, one line per read or mate. See the Bowtie 2 manual section on SAM output and the SAM specification for details about how to interpret the SAM file format.

### Paired-end example
To align paired-end reads included with Bowtie 2, stay in the same directory and run:

`bowtie2 -x lambda_virus -1 ../reads_1.fq -2 ../reads_2.fq -S eg2.sam`

This aligns a set of paired-end reads to the reference genome, with results written to the file eg2.sam. If you had used the filename eg1.sam here, your new results would have overwritten your previous results, so be careful!

### Local alignment example
To use local alignment to align some longer reads included with Bowtie 2, stay in the same directory and run:

`bowtie2 --local -x lambda_virus -U ../longreads.fq -S eg3.sam`

This aligns the long reads to the reference genome using local alignment, with results written to the file eg3.sam.

### Using SAMtools/BCFtools downstream
SAMtools is a collection of tools for manipulating and analyzing SAM and BAM alignment files. BCFtools is a collection of tools for calling variants and manipulating VCF and BCF files, and it is typically distributed with SAMtools. Both of these programs are available as module in Explorer. Using these tools together allows you to get from alignments in SAM format to variant calls in VCF format. Load the modules for samtools and bcftools into your environment (keeping bowtie2 active, as well) before proceeding.

Use samtools view to convert the SAM file for the paired-end alignment you performed into a BAM file. BAM is the binary format corresponding to the SAM text format. Run:

`samtools view -bS eg2.sam > eg2.bam`

Nothing will print to the console but you'll see a new file ending in .bam in your directory.

Use samtools sort to convert the BAM file to a sorted BAM file that can be used for calling variants.

`samtools sort eg2.bam -o eg2.sorted.bam`

We now have a sorted BAM file called eg2.sorted.bam. Sorted BAM is a useful format because the alignments are (a) compressed, which is convenient for long-term storage, and (b) sorted, which is conveneint for variant discovery. To generate variant calls in a format known as VCF, run:

`bcftools mpileup -f ../lambda_virus.fa eg2.sorted.bam | bcftools view -Ov - > eg2.raw.bcf`

Then to view the variants, run the command below, directing output to a file that you can name whatever you want:

`bcftools view eg2.raw.bcf > variants.txt`

View the file!
