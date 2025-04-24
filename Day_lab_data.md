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

We'll use STAR to index our genome. Although star is listed as a module on Explorer, their version doesn't work so we'll use the one I installed in a conda environment in the course "shared" directory. 

In the STAR indexing command, you'll need to customize to reflect the location of your files. The --genome-Dir option tells the program to look for the reference genome sequence in a directory named "ref". You may want to create such a directory for your reference and use that option. If your reference sequence is in the same folder as your gtf sequence and the script, you can remove the "--genome-Dir /ref" portion from the command. The --genomeFastaFiles option specifies that the reference genome is in FASTA format and should be left intact. 

Write a BATCH script for indexing, and note that it can take several hours for the task to complete. By default, task time has a limit of four hours, but we can increase it by adding another line to our header to increase the task time to 24 hours:

```html
#!/bin/bash
#SBATCH -J indxJH
#SBATCH -N 2
#SBATCH -n 16
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -t 24:00:00

#Activate the STAR conda environment (two steps!)


#Run the indexing command (you'll need to customize this to reflect the locations of your files)
STAR --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf --runThreadN 16

```

If it works, you'll see fourteen files added to your directory.

## Map reads to the genome

Let's align the reads to the indexed genome. Once again, think hard about how you'll need to customize the command to reflect the location of your files!

```html
#!/bin/bash
#SBATCH -J alnJH
#SBATCH -N 2
#SBATCH -n 12
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -t 24:00:00

#Activate the STAR conda environment (two steps!)

# Run a for loop with the align command 
for i in *R1.fastq; 
do
	name=$(basename ${i} _R1.fastq);
	STAR --runMode alignReads --genomeDir ../ref/ --outSAMtype BAM SortedByCoordinate --readFilesIn ${name}_R1.fastq ${name}_R2.fastq --outFileNamePrefix ../mapped/${name} --runThreadN 12;
done

```
This assignment is the first in which you have used a full-scale data set, giving you a taste of the time these operations can take and some of the snags that can trip you up.

Keep troubleshooting until it works! For help, use these resources:

*Copy your error message into chatgpt or claude

*Consult with your classmates on the Discussion forum on Canvas

