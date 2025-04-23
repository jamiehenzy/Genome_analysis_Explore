The program ANGSD produces a more reader-friendly variant calling output than the bcftools we used for the lambda phage example. I ran the script below after I had aligned the Ppar reads to the indexed genome using bowtie2. The folder in which I ran the script contained the six SAM files (one for each of the six fastq files aligned). Remember that you'll need to customize the script to reflect the locations of the various files you're using.

```
#!/bin/bash
#SBATCH -J VCangsd                          # Job name
#SBATCH -N 2                                # Number of nodes
#SBATCH -n 16                               # Number of tasks
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file

# Load necessary modules
module load samtools/1.21 bcftools/1.21 anaconda3/2024.06
source activate /courses/BIOL3411.202540/shared/conda_env/angsd_env

# SAM to BAM and then sort
for filename in *.sam
do
	base=$(basename $filename .sam)
	echo ${base}

	samtools view -bhS ${base}.sam | samtools sort -o ${base}.bam
done

# create bam filelist
ls *.bam > bam.filelist

# calculate genotype likelihoods
/courses/BIOL3411.202540/shared/conda_env/angsd_env/bin/angsd -bam bam.filelist -GL 1 -out genotype_likelihoods -doMaf 2 -SNP_pval 1e-2 -doMajorMinor 1
```
