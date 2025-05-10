# Working with genome sequences

### Know your sequences

When you think about sequence data, you'll want to ask yourself:
* Which species are the sequences from?
* Is this GENOMIC DNA (gDNA), representing all parts of the genome? 
* Or is this data from mRNA (cDNA), representing parts of the genome that are being expressed under some condition?

Genomic DNA sequences could be a completed assembly representing the reference genome for the species OR reads that have not yet been assembled into a complete genome, OR they could be reads from individuals from a species whose genome has already been assembled, but with the goal of identifying variation in the population ("resequencing").

Always make sure you're clear on what the sequences you are working with represent.

In this exercise you'll download **reference genome sequences** from several species:
* _B burgdorferei_ (Lyme disease bacterium); GCF_000008685.2
* _C elegans_ (lab worm); GCF_000002985.6
* _Parastichopus parvimensis_ (sea cucumber - you already have this one!)
* _Homo sapiens_ (most of us); GCA_000001405.29
* Lambda phage (a famous virus that infects bacteria); NC_001416.1

You'll also download some **other** files containing sequencing reads (fastq format), and a genome annotation file (.gtf) that contains information on genomic features in the human genome.

* reads_1.fq and reads_2.fq (from lambda phage)
* Homo_sapiens.GRCh38.113.gtf.gz

REMEMBER TO FIRST REQUEST A COMPUTING NODE.

### Download your data

**Genome sequences** can be downloaded from various sites. NCBI and Ensembl are two common sources. Download the B. berg., C. eleg., and lambda phage genomes from NCBI. For the first two, use the FTP link near the top of the genome assembly page. This will allow you to see a list of files available. Some of the files have information about the assembly, at least one will be a gtf file (containing genome feature information), and at least one will be the fasta sequence for the genome.

Let's walk through this for the first genome, Borrelia bergdorferei. Google "NCBI" to get to the main page. Next to the search window, from the dropdown menu select "Genome", and enter the accession number GCF_000008685.2. You'll be taken to a page that shows a table with the assembly number on the left. Click the assembly number (starts with ASM) and you'll be taken to a page of useful information about the genome, some of which you'll need for your documentation (see below). 

To download the genome, click the FTP (**f**ile **t**ransfer **p**rotocol) link near the top of the page to go to the list of files available. Scan to look for a fasta nucleotide file, which will end in ".fna" followed by ".gz", since it is compressed. Right-click to Copy Link Address and use the `wget` command to download the sequence to your Explorer directory. Et voilá!

Do the same for _C. elegans_ and Lambda phage, choosing the topmost assembly for each.

For the human genome, go to the Ensembl site (Google "ensembl human genome"). Be sure the correct assembly is shown in the upper left "Genome assembly" pane. Choose "Download DNA sequence (FASTA)". 

Whoa! There are lots of files to choose from. Select the "soft-masked" version (Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz) and download it using `wget`.

Now for our **other** files.

* GTF file for human genome: Go back to the Ensembl human genome page and look at the upper right pane, "Genome annotation". Click on the link to download a gtf file, and download Homo_sapiens.GRCh38.113.gtf.gz
* reads from Lambda phage project: Google "BenLangmead bowtie2 github"; click on the "examples" directory, then "reads". If your click on "reads_1.fq", you'll get a message about not being able to show such a big file. No problem! Right-click on "raw" (on the right) to Copy Link Address, then use `wget` to download the sequence. Do the same for "reads_2.fq".  

Unzip any of the sequences that end in gz using the gunzip command: `gunzip <file.gz>`

### Document your reference genome sequences

Create a directory in your student folder on Explorer for genomic data, with a **Genomes_info** file containing the information below for each genome. Most of the information can be found on the NCBI genomes pages.

*	Size of genome
*	Structure (such as, linear or circular, how many chromosomes, how many plasmids)
*	Estimated number of protein-coding genes
*	Where to access the genome sequence
* At what coverage it was sequenced

### Explore the sequences

To obtain general information about the genome sequences, you'll use a program, emboss, that is available as a conda environment. To use a conda environment you'll need to load the anaconda3/2024.06 module (2024.06 gives info on the version). 

`module load anaconda3/2024.06`

Check that it worked. To see which modules are loaded into your current session:

`module list`

Now you have to activate the conda environment, using the `source activate` command and the path to the files that are in the course "shared" folder:

`source activate /courses/BIOL3411.202540/shared/conda_env/emboss_env`

After entering the command, you'll see that your prompt includes the path of the environment.

Now you can use the emboss command `infoseq` to see what a file contains:

`infoseq <filename>`

* How many sequence files are in the worm genome?
* How many sequence files are in the Borrelia genome?
* Which sequence in the worm genome has the lowest GC%?
* Which sequence in Borrelia has the highest GC%?

**Remember to practice using tab-complete to save time and prevent mistakes!**

## Submission:
+ **Genomes_info** file that includes the five pieces of genome information for all five genomes
+ Well organized directories containing the five genomes, reads, and gtf file; arranged in user-friendly way
+ Saved session including hashed questions and answers to each of the five questions, along with the commands you used and output 

