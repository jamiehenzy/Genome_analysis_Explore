# What is a Genetic Variant?

The files we'll work with are in the Sea_cuke directory in the course "data" folder. You should see 6 files with a `.fastq.gz` extension and 1 genome file with a `.fna.gz` extension. **Copy** these into your own folder before working with them **(remember to request a computing node before copying – these are large files)**.

The data are from an excellent marine genomics study on sea cucumber population genetics ([Xuereb et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14589)). For simplification, we're using only a subsample of the reference genome and raw reads from only 6 individuals of the species _Parastichopus californicus_. 

The programs we'll use are the following:

+ samtools: allows us to filter and view our mapped data
+ bowtie2: to map our reads to the reference genome
+ cutadapt: will trim adaptor sequences from the reads
+ fastqc: used to view the quality of the read files
+ angsd: will find variants from the aligned reads

Some of these tools are available as modules. Others are installed in conda environments that are available in the course "shared" folder.

After you have loaded and/or activated the various programs, we're ready to get going. It's always good to first have a look at our data to make sure we know what (and where) everything is. First unzip the files, then use UNIX commands to determine:

+ The header information of the genome sequence
+ What the read files look like
+ How many reads are contained in the fastq files 

## Raw read quality control

Remember how in the first assignment, we got a feel for the quality of reads by using `grep` to search for occurrences of multiple N's? Well as you may have guessed, there ARE more sophisticated tools to determine the overall quality of a set of reads! The program `fastqc` can do this, with this command:

```html
$ fastqc SRR6805880.tiny.fastq

```
* Readout will say: 
  + Started analysis for SRR6805880.tiny.fastq
  + Analysis complete for SRR6805880.tiny.fastq


Let's check that it worked
```html
$ ls

Ppar_tinygenome.fna       SRR6805880.tiny_fastqc.zip  SRR6805883.tiny.fastq
SRR6805880.tiny.fastq        SRR6805881.tiny.fastq    SRR6805884.tiny.fastq
SRR6805880.tiny_fastqc.html  SRR6805882.tiny.fastq    SRR6805885.tiny.fastq


```

Looks good! Fastqc generated two new output files with extensions `.html` and a `.zip`. Notice also that the command worked on compressed files. It will work on either compressed OR uncompressed files.

You may have noticed that fastqc just prefixed the file name to different extensions for the two output files. We can take advantage of that by running fastqc on all our datafiles with the wildcard `*`. Let's run fastqc on the remaining files, then look at the output. 

```html
$ fastqc SRR680588*

```
You initially get an error message because fastqc doesn't see the .fastq file extension on some of our files. It simply skips these and moves on the the next file. 

To view the output of fastqc, we'll go to the OOD portal of Explorer and view our files from there. On OOD, use "change directory" to type in the file path of the directory you've been working in for this exercise. When you see your files displayed, you can use the dropdown menu next to any html file that was produced and either hit "view" or "download". 

Voilà! You should see a beautiful graphic output displaying various aspects of your data quality!

## Trimming to remove adapters

Use `less` or `cat` to scan one of the read files. You should notice that the first few bases of each read are exactly the same. Why would that be?

Sequencing machines sometimes read into part of the adaptor for various reasons, resulting in some bases from the adaptors at the ends of the reads. These are not part of the sequence-of-interest and need to be ignored in any analysis. 

There are many programs that can be used to "trim" such sequences from the reads. Keep in mind they are only trimming from the read data, and not from any actual molecules of DNA! The program "cutadapt" is relatively easy to run once we have identified our adaptor sequence, and takes the general form given below. Based on your peek into the sequencing reads, what do you think will replace SEQUENCETOTRIM in the command?


```html
$ cutadapt -g SEQUENCETOTRIM name_of_input_file -o name_of_output_file 

```

Let's do this on one of our files to test it out.

```html
cutadapt -g TGCAG SRR6805880.tiny.fastq -o SRR6805880.tiny_trimmed.fastq 

```
This works for a single file, but if we want to do it for all our read files we need to either do them all individually (slow and error prone) or use a for loop.

```html

for filename in *.tiny.fastq
do

  base=$(basename $filename .tiny.fastq)
  echo ${base}

  cutadapt -g TGCAG ${base}.tiny.fastq -o ${base}.tiny_trimmed.fastq 

done

```

Yay! You should see a little report for each of these files showing how many reads were trimmed and some other info (how long are the reads, etc.).

You can check if the trimmed files are there with:
```html
ls *trimmed*
```
And you can use `less` to check that the starts of the reads in a file are no longer identical, since the portion representing adapter sequence has been removed.

Our reads are now ready to be mapped to the genome.

## Building an index of our genome

First we have to index our genome. We'll do that with the bowtie2-build command, just as we did for lambda phage in the previous tutorial. 

Recall that we give bowtie2-build two things – the name of our genome, and a general name (prefix) to label the output files. A good practice is to keep the name of the output files the same as the original genome file (without the .fna.gz or .fna extension) to avoid confusion.

```html

bowtie2-build Ppar_tinygenome.fna Ppar_tinygenome

```
This should produce several output files with extensions including: .bt2 and rev.1.bt2 etc (six files in total). I like to make a new directory (e.g. "index") and move these files into it. 

## Map reads to the genome

Let's map those reads using a for loop. I've made a directory `index` that includes my .bt2 files.

```html
for filename in *.tiny_trimmed.fastq
do

  base=$(basename $filename .tiny_trimmed.fastq)
  echo ${base}

  bowtie2 -x index/Ppar_tinygenome -U ${base}.tiny_trimmed.fastq -S ${base}.sam

done

```

You should see a bunch of text telling you how well the reads mapped to the genome. For this example we're getting a low percentage (20-30%) because of how the genome and reads were subset for this exercise. The full genome and full read files have a much higher mapping rate (70-80%) than our subset. 

Use `less` to look at one of the SAM files, just as you did in the lambda phage exercise. See how much of the information you can decipher.


## sam to bam file conversion

The next step is to convert our sam file to a bam (Binary Alignment Map file). This gets our file ready to be read by angsd – the program we're going to use to call SNPs (type of variants known as "single nucleotide polymorphisms").

```html
for filename in *.sam
do

  base=$(basename $filename .sam)
  echo ${base}
  
  samtools view -bhS ${base}.sam | samtools sort -o ${base}.bam

done

```

## Genotype likelihoods

There are many ways and many programs that call genotypes. The program that we will use calculates genotype likelihoods, which account for uncertainty due to sequencing errors and/or mapping errors and is included in the package ANGSD. The purpose of this class is not to discuss which program is the "best", but to teach you to use some commonly used programs.

angsd needs a text file with the `.bam` file names listed. We can make this file by running the command:

```html

ls *.bam > bam.filelist

```

Look at the list:
```html
cat bam.filelist
```

We want to run the command `angsd` in the angsd program to calculate genotype likelihoods. Take a look at the files in the angsd_env folder to see if you can find the `angsd` command. Hint: commands are often listed in a folder called "bin"! 

Now alter the command below to reflect the path for the `angsd` command. Write a **bash script** that loads anaconda, activates the angsd conda environment, and runs the command. Think carefully about where the script is running from, and whether your script as written will allow the program to find the `angsd` command and your "bam.filelist" file.

```html

<path_to>/angsd -bam bam.filelist -GL 1 -out genotype_likelihoods -doMaf 2 -SNP_pval 1e-2 -doMajorMinor 1

```

If the script worked, you'll see two new files. The file with a .arg extension contains a record of the script we ran to generate the output,  and a .maf file contains the **m**inor **a**llele **f**requencies and is the main output file. We'll go over the components of the .maf file in class.


### In-class exercises
1. Map the untrimmed files to the genome. How do the alignments compare?

2. Use cutadapt to trim the sequences to 70 bp like they did in the Xuereb et al. 2018 paper. Write the output of cutadapt to a new file (named such that you know what it contains) and then map these new trimmed reads to the genome. How do they compare to our first trimmed reads?

3. Change the parameters of the angsd genotype likelihoods command. How many more/less SNPs do we recover if we lower or raise the SNP p-value? 

4. Run fastqc on our .trimmed reads and compare the html with the untrimmed files. 

## Submission
+ Saved session with hash-comments indicating each step performed
+ Session should show commands and output, with any extraneous output cleaned up
+ Text file with answers to in-class exercises (no need to show commands and output)
+ For question 3 (in-class) indicate two other p-values you tried and how many SNPs each value produced
