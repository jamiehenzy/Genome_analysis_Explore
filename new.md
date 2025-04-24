This last component of the "Explore" assignment pulls together what you've learned from the previous components while introducing you to another powerful tool for indexing genomes and aligning reads – STAR (Spliced Transcripts Alignment to a Reference). When reads are from RNA transcripts, a given read may span two exons which can be separated by thousands of base pairs! As the name implies, STAR is especially powerful for aligning spliced reads to genomes. 

For this assignment, you are given less guidance because you are expected to use knowledge you have gained in previous assignments and class. However, please always remember to request a computing node before beginning any work.

# Day lab RNA-seq data

The lab of NU researcher Tovah Day generated some RNA-seq data from an experiment you'll learn more about for the differential gene expression assignment. You'll find the files in our course data folder. For now, we want to:

+ assess the quality of the reads and do any necessary modifications
+ index the human genome (with STAR)
+ align the reads (with STAR)

## Gather the files

You'll need the **human genome fasta file** and **gtf file**, both of which you downloaded in the Genomes_tutorial.

You'll also need to write a script (bash or BATCH) to copy the **Day_lab_reads** from the course "data" folder into your own directory.

## Raw read quality control

What's the program you can use to generate those snazzy quality control reports? In a previous tutorial, you used a wildcard in the filenames to run the program on all the files. Run it and look for any red flags in the reports.

```html
$ <program_name_here> <file_name>*

```
Remember that to view the reports, you need to access the html files through the OOD portal of Explorer. 

*Are there any red flags? Use chat or claude to help you figure out if you should be concerned.
*Do you need to trim the sequences?

## Building an index of our genome

We'll use STAR to index our genome. Write a BATCH script to do this, and note that it can take several hours for the task to complete. By default, task time has a limit of four hours, but we can increase it by adding a line to our header:

```html
#!/bin/bash
#SBATCH -J indxJH
#SBATCH -N 2
#SBATCH -n 16
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -t 24:00:00

#Load the STAR module
module load star/2.7.11b

#Run the indexing command (you'll need to customize to reflect the locations of your files)
STAR --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf --runThreadN 16

![image](https://github.com/user-attachments/assets/4ea405d9-6674-4227-b3df-5bf475be37fc)

```
If it works, you'll see fourteen files added to your directory.

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

We want to run the command `angsd` in angsd program to calculate genotype likelihoods. Take a look at the files in the angsd_env folder to see if you can find the `angsd` command. Hint: commands are often listed in a folder called "bin"! 

Now alter the command given below to reflect the path for the `angsd` command. Write a bash script that loads anaconda3, activates the angsd conda environment, and runs the command. Think carefully about where the script is running from, and whether the program will be able to find the `angsd` command and your "bam.filelist" file.

```html

<path_to>/angsd -bam bam.filelist -GL 1 -out genotype_likelihoods -doMaf 2 -SNP_pval 1e-2 -doMajorMinor 1

```

If the script worked, you'll see two new files. The file with a .arg extension contains a record of the script we ran to generate the output,  and a .maf file contains the minor allele frequencies and is the main output file. We'll go over the components of the .maf file in class.


### In-class exercises
1. Map the untrimmed files to the genome. How do the alignments compare?

2. Use cutadapt to trim the sequences to 70 bp like they did in the Xuereb et al. 2018 paper. Write the output of cutadapt to an .70bp.trimmed.fastq.gz and then map these 70bp, trimmed reads to the genome. How do they compare to our .trimmed reads?

3. Change the parameters of the angsd genotype likelihoods command. How many more/less SNPs do we recover if we lower or raise the SNP p-value? 

4. Run fastqc on our .trimmed reads and compare the html with the untrimmed files. 
