This last component of the "Explore" assignment pulls together what you've learned from the previous components while introducing you to another powerful tool for indexing genomes and aligning reads – STAR (Spliced Transcripts Alignment to a Reference). When reads are from RNA transcripts, a given read may span two exons which can be separated by thousands of base pairs! As the name implies, STAR is especially powerful for aligning spliced reads to genomes. 

For this assignment, you are given less guidance because you are expected to use knowledge you have gained in previous assignments and class. However, please always remember to request a computing node before beginning any work.

# RNA-seq data

In the Genomes tutorial, you learned to download genome sequences from various sources. You also downloaded a couple of files of sequencing reads from GitHub. However, a major repository of sequencing reads is the GEO (Gene Expression Omnibus) database. You can use a tool called sratoolkit, available as a module on Explorer, to download a set of reads from an experiment involving human cell lines. Then you'll align (map) those reads to the indexed human genome sequence.

+ assess the quality of the reads and do any necessary modifications
+ index the human genome (with STAR)
+ align the reads (with STAR)

## Gather the files

Genome: We'll need the **human genome fasta file** and **gtf file**, both of which you downloaded in the Genomes_tutorial.

Reads to map: From this [GEO dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778), we want the following accession numbers: 
SRR1039516
SRR1039517
SRR1039508
SRR1039509
SRR1039512
SRR1039513
SRR1039520
SRR1039521

These records correspond to reads from the untreated and treated cells of four different cell lines. Each record has millions of reads. We first have to "prefetch" the records, then convert them into fastq format. To fetch all the records in one run, we need to make a text file with all eight accession numbers and name it "accessions.txt". Create this file using nano (or "echo" into a new file).

Next we'll write a batch script to run the job on Explorer. Be sure to customize the script to suit your needs.

```html
#!/bin/bash
#SBATCH -J pfJH
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -t 24:00:00

# load sratoolkit
module load sratoolkit/12Dec2024

# carry out prefetch of the list of accessions
prefetch --option-file accessions.txt

```

Check the new output and error files in your directory. If everything is running smoothly, your output file will continue to update progress. It takes a few minutes or more for all eight files to be processed. After it runs, there will be a new directory, "sra_files", with eight files that are in a format that you cannot read. We need to convert these to fastq format. Again, be sure to make changes as needed to the script:

```html
#!/bin/bash
#SBATCH -J fqJH
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o output_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -t 24:00:00

# load sratoolkit
module load sratoolkit/12Dec2024

# convert to fastq files
for acc in SRR1039508  SRR1039509  SRR1039512  SRR1039513  SRR1039516  SRR1039517  SRR1039520  SRR1039521;

do
    fasterq-dump --split-files --outdir ./fastq_files "$acc"
done
```
When the job has finished (it takes awhile), you'll see a directory called, "fastq_files" that contains eight files, each containing millions of fastq reads. (You might also see that your progress log appears in your error file instead of your output file!)

## Raw read quality control

Use the program you used previously to examine the quality of the reads. You do not need to look at all eight reports. Just pick two or three. Since the samples were all run at the same time, the quality should be similar for all of them.

+ Are there any red flags? Use chat or Claude to help you figure out if you should be concerned.
+ Do you need to trim the sequences?

## Building an index of our genome

We'll use STAR to index our genome. Although star is listed as a module on Explorer, their version doesn't work so we'll use the one I installed in a conda environment in the course "shared" directory. Instead of "source activate \<pathway>" use "conda activate \<pathway>".

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
STAR --runMode genomeGenerate --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.113.gtf --runThreadN 16
```

If it works, you'll see sixteen files added to your directory.

## Map reads to the genome

Let's align the reads to the indexed genome. Once again, **think hard** about how you'll need to customize the command to reflect the location of your files! This is a complicated script and we'll go through the various components in class.

```html
#!/bin/bash
#SBATCH -J alignJH
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o output_%j.txt                   
#SBATCH -e error_%j.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=12:00:00

# load anaconda and activate star conda environment – you know how to do this!


# Directory containing the genome index
GENOME_DIR="GenomeDir"

# Directory containing the fastq files
FASTQ_DIR="SRA_work/fastq_files"

# Create output directory if it doesn't exist
mkdir -p aligned_output

# Loop through all sample prefixes
for sample in $(ls ${FASTQ_DIR}/*_1.fastq | sed 's/_1.fastq$//' | xargs -n1 basename); do
    
    # Define input files
    READ1="${FASTQ_DIR}/${sample}_1.fastq"
    READ2="${FASTQ_DIR}/${sample}_2.fastq"
    
    # Define output prefix
    OUT_PREFIX="aligned_output/${sample}_"
    
    # Run STAR alignment
    STAR --genomeDir ${GENOME_DIR} \
         --readFilesIn ${READ1} ${READ2} \
         --runThreadN 16 \
         --outFileNamePrefix ${OUT_PREFIX} \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 40000000000
         
done

```
Find the new directory "aligned_output" and use `ls -l` to view the contents in long form. Look at the column of sizes to check that all of the files that end in .bam have content. These should be the largest files. We'll use these mapped reads for our upcoming differential gene expression (DGE) assignment.

This assignment is the first in which you have used a full-scale data set, giving you a taste of the time these operations can take and some of the snags that can trip you up.

Keep troubleshooting until it works! For help, use these resources:

+ Copy your error message into chatgpt or Claude
+ Consult with your classmates on the Discussion forum on Canvas

## Submission
+ all files associated with tutorial, organized in user-friendly way
+ all scripts you ran
+ text file with notes on general impressions of quality score report

## Don't forget! 
+ README file
+ no extraneous lines of text (limit output to 10 lines max)
+ permissions set so that I can access everything
+ clean front porch

