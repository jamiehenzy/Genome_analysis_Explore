# What is a Genetic Variant?

The files we'll work with are in the course "data" folder. You should see 6 files with a `.fastq.gz` extension and 1 tiny genome file with a `.fna.gz` extension. Copy these into your own folder before working with them (remember to request a computing node before copying – these are large files).

The data are from an excellent marine genomics study on sea cucumber population genetics ([Xuereb et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14589)). For simplification, we're using only a subsample of the reference genome and raw reads from only 5 individuals of the species _Parastichopus californicus_. 

The programs we'll use are the following:

+ samtools: allows us to filter and view our mapped data
+ bowtie2: to map our reads to the reference genome
+ cutadapt: will trim adaptor sequences from the reads
+ fastqc: used to view the quality of the read files
+ angsd: will find variants from the aligned reads

Some of these tools are available as modules. Others are installed in conda environments that are available in the course "shared" folder.
What additional module do you need to load in order to activate a conda environment?

After you have loaded and/or activated the various programs, we're ready to get going. It's always good to first have a look at our data to make sure we know what (and where) everything is. Use UNIX commands to determine:

+ The header information of the genome sequence
+ What the read files look like

How could you use simple UNIX commands to determine how many reads are contained in the fastq files? 

## Raw read quality control

Remember how in the first assignment, we got a feel for the quality of reads by using `grep` to search for occurrences of multiple N's? Well as you may have guessed, there ARE more sophisticated tools to determine the overall quality of a set of reads! The program `fastqc` can do this, with this command:

```html
$ fastqc SRR6805880.tiny.fastq.gz

```
* Readout will say: 
  + Started analysis for SRR6805880.tiny.fastq.gz
  + Analysis complete for SRR6805880.tiny.fastq.gz


Let's check that it worked
```html
$ ls

Ppar_tinygenome.fna.gz       SRR6805880.tiny_fastqc.zip  SRR6805883.tiny.fastq.gz
SRR6805880.tiny.fastq.gz     SRR6805881.tiny.fastq.gz    SRR6805884.tiny.fastq.gz
SRR6805880.tiny_fastqc.html  SRR6805882.tiny.fastq.gz    SRR6805885.tiny.fastq.gz


```

Looks good! Fastqc generated two outputs for us – a `.html` and a `.zip` directory. Notice also that the command worked on compressed files. It will work on either compressed OR uncompressed files.

Let's run fastqc on the remaining files, then look at the output. You may have noticed fastqc just used the same file name to produce our output with different extensions. We can take advantage of that by running fastqc on all our datafiles with the wildcard `*`.

```html
$ fastqc SRR680588*

```
You initially get an error message because fastqc doesn't see the .fastq file extension on some of our files. It simply skips these and moves on the the next file. 

To view the output of fastqc, we'll go to the OOD portal of Explorer and view our files from there. Use "change directory" to type in the file path of the directory you've been working in for this exercise. When you see your files displayed on OOD, you can use the dropdown menu next to any html file that was produced and either hit "view" or "download". 

Voilá! You should see a beautiful graphic output displaying various aspects of your data quality!

## Trimming to remove adapters

Use `less` or `cat` to scan one of the read files. You should notice that the first few bases of each read are exactly the same. Why would that be?

Sequencing machines sometimes read into part of the adaptor for various reasons, resulting in some bases from the adaptors at the ends of the reads. These are not part of the sequence-of-interest and need to be ignored in any analysis. 

There are many programs that can be used to "trim" such sequences from the reads. Keep in mind they are only trimming from the read data, and not from any actual molecules of DNA. The program "cutadapt" is relatively easy to run with the code below, once we have identified our adaptor sequence, and takes the general form below.


```html
$ cutadapt -g SEQUENCETOTRIM -o name_of_input_file name_of_output_file 

```

Let's do this on one of our files to test it out.

```html
cutadapt -g TGCAG SRR6805880.tiny.fastq.gz -o SRR6805880.tiny_trimmed.fastq.gz 

```
This works for a single file, but if we want to do it for all our read files we need to either do them all individually (slow and error prone) or use a for loop.

```html

for filename in *.tiny.fastq.gz
do

  base=$(basename $filename .tiny.fastq.gz)
  echo ${base}

  cutadapt -g TGCAG ${base}.tiny.fastq.gz -o ${base}.tiny_trimmed.fastq.gz 

done

```

Yay! You should see a little report for each of these files showing how many reads were trimmed and some other info (how long are the reads, etc.).

You can check if the trimmed files are there with:
```html
ls *trimmed*
```

Our reads are now ready to be mapped to the genome.

## Building an index of our genome

First we have to index our genome. We'll do that with the bowtie2-build command. This will generate a lot of files that describe different aspects of our genome

We give bowtie2-build two things, the name of our genome, and a general name to label the output files. I always keep the name of the output files the same as the original genome file (without the .fna.gz extension) to avoid confusion (what's this file for?).

```html

bowtie2-build Ppar_tinygenome.fna.gz Ppar_tinygenome

```
This should produce several output files with extensions including: .bt2 and rev.1.bt2 etc (six files in total)

## Map reads to the genome

Let's map those reads using a for loop

```html
for filename in *.tiny_trimmed.fastq.gz
do

  base=$(basename $filename .tiny_trimmed.fastq.gz)
  echo ${base}

  bowtie2 -x Ppar_tinygenome -U ${base}.tiny_trimmed.fastq.gz -S ${base}.sam

done

```

You should see a bunch of text telling you all about how well our reads mapped to the genome. For this example we're getting a low percentage (20-30%) because of how the genome and reads were subset for this exercise. The full genome and full read files have a much higher mapping rate (70-80%) than our subset. 

You'll also notice that we have made a bunch of .sam files. .sam stands for Sequence Alignment Map file. Let's use `less` to look at one of these files using `less`

There are several columns of data in a sam file

## sam to bam file conversion

The next step is to convert our sam file to a bam (Binary Alignment Map file). This gets our file ready to be read by angsd the program we're going to use to call SNPs.

```html
for filename in *.sam
do

  base=$(basename $filename .sam)
  echo ${base}
  
  samtools view -bhS ${base}.sam | samtools sort -o ${base}.bam

done

```

## Genotype likelihoods

There are many ways and many programs that call genotypes. The program that we will use calculates genotype likelihoods, which account for uncertainty due to sequencing errors and/or mapping errors and is one of several programs in the package ANGSD. The purpose of this class is not to discuss which program is the "best", but to teach you to use some commonly used programs.

angsd needs a text file with the `.bam` file names listed. We can make that by running the command below

```html

ls *.bam > bam.filelist

```

Look at the list:
```html
cat bam.filelist
```

Run the following code to calculate genotype likelihoods

```html

../../angsd/angsd -bam bam.filelist -GL 1 -out genotype_likelihoods -doMaf 2 -SNP_pval 1e-2 -doMajorMinor 1

```

This will generate two files, one with a .arg extension, this has a record of the script we ran to generate the output, and a .maf file that will give you the minor allele frequencies and is the main output file. If you see these two files, Yay!! We did it!



### Suggested Exercises
> For our coding session you can re-run through the above code as it is written. You can also do the below suggestions to extend or modify what we did in Tuesdays class. 

> A possible answer is located beneath each activity, but it's possible you will correctly perform the suggestion in a different way. 


> 1. map the untrimmed files to the genome. How do the alignments compare?

<details><summary><span style="color: purple;">Solution</span></summary>
<p>

> 1. As a for loop:
> `for filename in *tiny.fastq.gz; do base=$(basename $filename .tiny.fastq.gz); echo=${base}; bowtie2 -x Ppar_tinygenome -U ${base}.tiny.fastq.gz -S ${base}.nottrimmed.sam; done`

> You should see something that by trimming the adapters off we get a higher overall mapping rate:

> ![results](./figs/week4/trimm_notrimm_thatisthequestion.jpg)

</p>
</details>
&nbsp;



> 2. Run the mapping for loop as a shell script using bash (i.e., store the for loop in a text editor (NANOs or other) and execute the .sh script with bash)

<details><summary><span style="color: purple;">Solution</span></summary>
<p>

> 2. this can be done by copying and pasting the for loop in a text editor that you save as for example `map_samples_bowtie2.sh`. This script is then executed by `bash map_samples_bowtie2.sh`

</p>
</details>
&nbsp;


> 3. use cutadapt to trim the sequences to 70 bp like they did in the Xuereb et al. 2018 paper. Write the output of cutadapt to an .70bp.trimmed.fastq.gz and then map these 70bp, trimmed reads to the genome. How do they compare to our .trimmed reads?

<details><summary><span style="color: purple;">Solution</span></summary>
<p>
> 3. to find the parameter for maximum read length in cutadapt: `cutadapt - help` There are a few ways to do this.
> `cutadapt -g TGCAG ${base}.tiny.fastq.gz -u 70 -o ${base}.tiny_70bp_trimmed.fastq.gz` 

</p>
</details>
&nbsp;


> 4. change the parameters of the angsd genotype likelihoods command. How many more/less SNPs do we recover if we lower or raise the SNP p-value? To see what the other parameters do run `../../angsd/angsd -h

<details><summary><span style="color: purple;">Solution</span></summary>
<p>

> 4. If we remove the `-SNP_pval` command entirely we get ~72000 sites retained! Wow! That seems like a lot given our ~20% maping rate. If you instead increase the p-value threshold to 1e-3 we find 3 SNPs.

</p>
</details>
&nbsp;

> 5. Run fastqc on our .trimmed reads and compare the html with the untrimmed files. 

<details><summary><span style="color: purple;">Solution</span></summary>
<p>

> 5. We should no longer see the red error flag for the per base sequence quality or base pairs conten. code: fastqc *trimmed.fastq.gz 

</p>
</details>
&nbsp;
