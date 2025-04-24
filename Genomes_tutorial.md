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

To download the genome, click the FTP link to go to the list of files available. Scan to look for a fasta nucleotide file, which will end in ".fna" followed by ".gz", since it is compressed. Right-click to Copy Link Address and use the `wget` command to download the sequence to your Explorer directory. Eh voilá!

Do the same for C. elegans and Lambda phage.

For the human genome, go to the Ensembl site (Google "ensembl human genome"). Be sure the correct assembly is shown in the upper left "Genome assembly" pane. Choose "Download DNA sequence (FASTA)". 

Whoa! There are lots of files to choose from. Select the "soft-masked" version (Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz) and download it using `wget`.

Now for our **other** files.

* GTF file for human genome: Go back to the Ensembl human genome page and look at the upper right pane, "Genome annotation". Click on the link to download a gtf file, and download Homo_sapiens.GRCh38.113.gtf.gz
* reads from Lambda phage project: Google "BenLangmead bowtie2 github"; click on the "examples" directory, then "reads". If your click on "reads_1.fq", you'll get a message about not being able to show such a big file. No problem! Right-click on "raw" (on the right) to Copy Link Address, then use `wget` to download the sequence. Do the same for "reads_2.fq".  

### Document your reference genome sequences

Create a directory in your student folder on Explorer for genomic data, with a README file containing the information below for each genome. Most of the information can be found on the NCBI genomes pages.

*	Size of genome
*	Structure (such as, linear or circular, how many chromosomes, how many plasmids)
*	Estimated number of protein-coding genes
*	Where to access the genome sequence
* At what coverage it was sequenced

### Explore the sequences

To obtain information about the genome sequences, you'll use a program, emboss, that is available as a conda environment and the other, as a **module** on Explorer: emboss and seqtk.

You need to **load** the seqtk module into your workspace:

`module load seqtk`

And you'll need to load the anaconda
Check that it worked. To see which modules are loaded into your current session:

`module list`

You should see emboss and seqtk among those shown.

You can use the emboss command `infoseq` to see what a file contains:

`infoseq <filename>`

* Does the command work on zipped (compressed) files?
* How many sequence files are in the worm genome?
* How many sequence files are in the Borrelia genome?
* Which sequence in the worm genome has the lowest GC%?
* Which sequence in Borrelia has the highest GC%?

What if you wanted to extract one of the sequences into a new file for further analysis. For example, let's extract the Borrelia sequence with the highest GC%. You can use a command from the program **seqtk** to extract a sequence. First you need to put the sequence identifier (usually the accession number) into a file. You can do this quickly using echo. 

`echo '<identifier>'  >  <make-up-a-filename>`

For example, it might look like this:

`echo 'NC_008524.2' > list.txt`

**Note: I'm using a different accession number as an example**

This will create a file called list.txt that contains just one line: NC_008524.2.

Then use the command that tells seqtk to extract that sequence from the larger file and store it in a new file by itself:

`seqtk subseq <file-with-many-sequences> list.txt > NC_008524.2.txt`

The sequence portion of a file is often presented as all one line, so using the `head` command just outputs the entire sequence. Which command allows you to view the contents in a more controlled way? Use it to note the first three and last three letters of the sequence (write them down).

Now we can generate the **reverse complement** of this sequence, representing the complementary strand. We'll put it into a new file:

`revseq NC_008524.2.txt NC_008524.2.rev`

**Remember to keep using tab-complete, and check after commands that the new file is added to your directory by using ls.**

Look at the contents of the .rev file. Are you convinced this is the reverse complement?

Now perform the same set of operations on the worm sequence that had the lowest GC%!
