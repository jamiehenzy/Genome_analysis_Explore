# Working with genome sequences

### Know your sequences

When you think about sequence data, you'll want to ask yourself:
* Which species are the sequences from?
* Is this GENOMIC DNA (gDNA), representing all parts of the genome? 
* Or is this data from mRNA (cDNA), representing parts of the genome that are being expressed under some condition?

Genomic DNA sequences could be reads that have not yet been assembled into a complete genome, OR they could be reads from individuals from a species whose genome has already been assembled, but with the goal of identifying variation in the population ("resequencing").

Always make sure you're clear on what the sequences you are working with represent.

In this exercise you'll download genome sequences from several species:
* _B burgdorferei_ (Lyme disease bacterium)
* _C elegans_ (lab worm)
* _Parastichopus parvimensis_ (sea cucumber - you already have this one!)

You'll also download a couple of tools from NCBI that allow you to readily download a variety of genomes, genes, and other genomic information.

### Document your sequences

I.  Create a directory in your student folder on Discovery for genomic data, with a README file:

Create a file called README.txt that contains basic information on each genome you add to the directory, including:
*	Size of genome
*	Structure (such as, linear or circular, how many chromosomes, how many plasmids)
*	Estimated number of protein-coding genes
*	Where to access the genome sequence
* At what coverage it was sequenced

Here are the accession numbers for the genome sequences you should use:

B. burgdorferei: GCF_000008685.2

C. elegans: GCF_000002985.6 

### Download the download tools

For this step, you'll work in your temporary "scratch" directory. The ncbi tools you are downloading take up a lot of space and after using them this initial time, there's no reason to keep them in your files.

II.  Download the **ncbi_datasets** tools. These are handy for downloading lots of genome data. Log onto Discovery and cd to your “scratch” account:

`cd /scratch/<username>`

Request a computing node.

Download the NCBI command line tools (datasets and dataformat)

`curl -o datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets`

Then:

`curl -o dataformat https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat`

Make them executable: 

`chmod +x datasets dataformat`

***Use the `ls` command to check that they are there (should see “datasets” and “dataformat” in the directory).***
  
### Use the new tools to download the genomes

III.  Download the genome sequences to your scratch account.

For downloading genomes you'll want to use the accession number for each genome (for example, GCF000002985.6). The NCBI tools make this easy. For each genome accession number, use this command:

`./datasets download genome accession <accession-number>`

After the download, use `ls` to see that an NCBI zip file is now in the directory. Unzip it by entering:

`unzip <filename>.zip`

You’ll see some files being “inflated”. When the command prompt returns, use the `ls` command again and you should see a directory ‘ncbi_dataset’ that contains the unzipped files.

## Stop and think

Look at the name of the file you unzipped into several more files, and at the name of the directory in which those files live. Now think about what would happen if you repeat the same procedure for the second genome accession number. 

Most programs give you options for the output of whatever function you're performing. In this case, and option `--filename` allows you to specify what you want the zipped file to be named. That way, you won't overwrite the first set of files when you download the second set! Just make sure that your output filename ends with .zip. For example, for the Borrelia genome download, you might modify the command above like so:

`./datasets download genome accession GCF_000008685.2 --filename bbur.zip`

Do the same for the _C. elegans_ genome. 

Similarly, when you unzip the files, it creates a directory named ncbi_dataset. To load the files into specific directories so that the one set of files does not replace the other, you can use the `-d` flag with the unzip command. First, make a new directory named something like, "worm_unzip". Then use this command to unzip the files to the newly created directory:

`unzip -d worm_unzip/ worm.zip`

Note the syntax of the command: {do this} {target for output will be a directory} {target directory name} {file to unzip} 

Explore the output files. Which file is the actual genomic sequence? What are the other files?

### First foray into modules!

Before proceeding, move your genome files from your scratch folder to the new directory you set up in your student folder. For each genome, there should be three files that you move.

IV.  Analysis
To analyze the genome sequences, you'll use two programs that are available as **modules** on Discovery: seqtk and emboss.

You need to **load** the modules into your workspace (scratch) to use them.

`module load emboss`

`module load seqtk`

Check that it worked. To see which modules are loaded into your current session:

`module list`

You should see emboss and seqtk among those shown.

When you download your organism's genome sequence, multiple fasta sequences are usually included in one file. You can use the emboss command 'infoseq' to see what a file contains:

`infoseq <filename>`

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
