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

> @HD VN:1.0  SO:unsorted
@SQ SN:gi|9626243|ref|NC_001416.1|  LN:48502
@PG ID:bowtie2  PN:bowtie2  VN:2.0.1
r1  0   gi|9626243|ref|NC_001416.1| 18401   42  122M    *   0   0   TGAATGCGAACTCCGGGACGCTCAGTAATGTGACGATAGCTGAAAACTGTACGATAAACNGTACGCTGAGGGCAGAAAAAATCGTCGGGGACATTNTAAAGGCGGCGAGCGCGGCTTTTCCG  +"@6<:27(F&5)9"B):%B+A-%5A?2$HCB0B+0=D<7E/<.03#!.F77@6B==?C"7>;))%;,3-$.A06+<-1/@@?,26">=?*@'0;$:;??G+:#+(A?9+10!8!?()?7C>  AS:i:-5 XN:i:0  XM:i:3  XO:i:0  XG:i:0  NM:i:3  MD:Z:59G13G21G26    YT:Z:UU
r2  0   gi|9626243|ref|NC_001416.1| 8886    42  275M    *   0   0   NTTNTGATGCGGGCTTGTGGAGTTCAGCCGATCTGACTTATGTCATTACCTATGAAATGTGAGGACGCTATGCCTGTACCAAATCCTACAATGCCGGTGAAAGGTGCCGGGATCACCCTGTGGGTTTATAAGGGGATCGGTGACCCCTACGCGAATCCGCTTTCAGACGTTGACTGGTCGCGTCTGGCAAAAGTTAAAGACCTGACGCCCGGCGAACTGACCGCTGAGNCCTATGACGACAGCTATCTCGATGATGAAGATGCAGACTGGACTGC (#!!'+!$""%+(+)'%)%!+!(&++)''"#"#&#"!'!("%'""("+&%$%*%%#$%#%#!)*'(#")(($&$'&%+&#%*)*#*%*')(%+!%%*"$%"#+)$&&+)&)*+!"*)!*!("&&"*#+"&"'(%)*("'!$*!!%$&&&$!!&&"(*"$&"#&!$%'%"#)$#+%*+)!&*)+(""#!)!%*#"*)*')&")($+*%%)!*)!('(%""+%"$##"#+(('!*(($*'!"*('"+)&%#&$+('**$$&+*&!#%)')'(+(!%+ AS:i:-14    XN:i:0  XM:i:8  XO:i:0  XG:i:0  NM:i:8  MD:Z:0A0C0G0A108C23G9T81T46 YT:Z:UU
r3  16  gi|9626243|ref|NC_001416.1| 11599   42  338M    *   0   0   GGGCGCGTTACTGGGATGATCGTGAAAAGGCCCGTCTTGCGCTTGAAGCCGCCCGAAAGAAGGCTGAGCAGCAGACTCAAGAGGAGAAAAATGCGCAGCAGCGGAGCGATACCGAAGCGTCACGGCTGAAATATACCGAAGAGGCGCAGAAGGCTNACGAACGGCTGCAGACGCCGCTGCAGAAATATACCGCCCGTCAGGAAGAACTGANCAAGGCACNGAAAGACGGGAAAATCCTGCAGGCGGATTACAACACGCTGATGGCGGCGGCGAAAAAGGATTATGAAGCGACGCTGTAAAAGCCGAAACAGTCCAGCGTGAAGGTGTCTGCGGGCGAT  7F$%6=$:9B@/F'>=?!D?@0(:A*)7/>9C>6#1<6:C(.CC;#.;>;2'$4D:?&B!>689?(0(G7+0=@37F)GG=>?958.D2E04C<E,*AD%G0.%$+A:'H;?8<72:88?E6((CF)6DF#.)=>B>D-="C'B080E'5BH"77':"@70#4%A5=6.2/1>;9"&-H6)=$/0;5E:<8G!@::1?2DC7C*;@*#.1C0.D>H/20,!"C-#,6@%<+<D(AG-).?&#0.00'@)/F8?B!&"170,)>:?<A7#1(A@0E#&A.*DC.E")AH"+.,5,2>5"2?:G,F"D0B8D-6$65D<D!A/38860.*4;4B<*31?6  AS:i:-22    XN:i:0  XM:i:8  XO:i:0  XG:i:0  NM:i:8  MD:Z:80C4C16A52T23G30A8T76A41   YT:Z:UU
r4  0   gi|9626243|ref|NC_001416.1| 40075   42  184M    *   0   0   GGGCCAATGCGCTTACTGATGCGGAATTACGCCGTAAGGCCGCAGATGAGCTTGTCCATATGACTGCGAGAATTAACNGTGGTGAGGCGATCCCTGAACCAGTAAAACAACTTCCTGTCATGGGCGGTAGACCTCTAAATCGTGCACAGGCTCTGGCGAAGATCGCAGAAATCAAAGCTAAGT(=8B)GD04*G%&4F,1'A>.C&7=F$,+#6!))43C,5/5+)?-/0>/D3=-,2/+.1?@->;)00!'3!7BH$G)HG+ADC'#-9F)7<7"$?&.>0)@5;4,!0-#C!15CF8&HB+B==H>7,/)C5)5*+(F5A%D,EA<(>G9E0>7&/E?4%;#'92)<5+@7:A.(BG@BG86@.G AS:i:-1 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:77C106 YT:Z:UU
r5  0   gi|9626243|ref|NC_001416.1| 48010   42  138M    *   0   0   GTCAGGAAAGTGGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGAGGGAAGAACGCGGGATGTTCATTCTTCATCACTTTTAATTGATGTATATGCTCTCTT  9''%<D)A03E1-*7=),:F/0!6,D9:H,<9D%:0B(%'E,(8EFG$E89B$27G8F*2+4,-!,0D5()&=(FGG:5;3*@/.0F-G#5#3->('FDFEG?)5.!)"AGADB3?6(@H(:B<>6!>;>6>G,."?%  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:138    YT:Z:UU
r6  16  gi|9626243|ref|NC_001416.1| 41607   42  72M2D119M   *   0   0   TCGATTTGCAAATACCGGAACATCTCGGTAACTGCATATTCTGCATTAAAAAATCAACGCAAAAAATCGGACGCCTGCAAAGATGAGGAGGGATTGCAGCGTGTTTTTAATGAGGTCATCACGGGATNCCATGTGCGTGACGGNCATCGGGAAACGCCAAAGGAGATTATGTACCGAGGAAGAATGTCGCT 1H#G;H"$E*E#&"*)2%66?=9/9'=;4)4/>@%+5#@#$4A*!<D=="8#1*A9BA=:(1+#C&.#(3#H=9E)AC*5,AC#E'536*2?)H14?>9'B=7(3H/B:+A:8%1-+#(E%&$$&14"76D?>7(&20H5%*&CF8!G5B+A4F$7(:"'?0$?G+$)B-?2<0<F=D!38BH,%=8&5@+ AS:i:-13    XN:i:0  XM:i:2  XO:i:1  XG:i:2  NM:i:4  MD:Z:72^TT55C15A47  YT:Z:UU
r7  16  gi|9626243|ref|NC_001416.1| 4692    42  143M    *   0   0   TCAGCCGGACGCGGGCGCTGCAGCCGTACTCGGGGATGACCGGTTACAACGGCATTATCGCCCGTCTGCAACAGGCTGCCAGCGATCCGATGGTGGACAGCATTCTGCTCGATATGGACANGCCCGGCGGGATGGTGGCGGGG -"/@*7A0)>2,AAH@&"%B)*5*23B/,)90.B@%=FE,E063C9?,:26$-0:,.,1849'4.;F>FA;76+5&$<C":$!A*,<B,<)@<'85D%C*:)30@85;?.B$05=@95DCDH<53!8G:F:B7/A.E':434> AS:i:-6 XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:98G21C22   YT:Z:UU

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
