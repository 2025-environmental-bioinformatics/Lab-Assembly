# Today in class:

1. [Assembly lecture](https://docs.google.com/presentation/d/19fLb5RsHCdPkPE5QT7y7EJNzvQ2R4d8fH94HHrAu04c/edit#slide=id.p)
3. [Transcriptome assembly lecture](https://docs.google.com/presentation/d/11oK2w1xru0Ih-6tta43E70rnnN49GNv-FydPjlywROk/edit)
4. Transcriptome assembly lab (below)
5. Trimmomatic data cleaning

# Lab-Assembly

First, log onto Poseidon and set up a tmux window so we won't be interrupted:

`tmux new -s txm`

Navigate into your user folder for the class, and pull the "Txm-lab" Git repository:

```
git clone git@github.com:2023-MIT-Environmental-Bioinformatics/Lab-Assembly.git
cd Lab-Assembly
```

Now, let's set up a conda environment to play in:

```
conda create --name trinity -c bioconda trinity
```

In this lab repo, you should see the following fastq file:
`SRR8956770_1_sub_mod.fq`

This is a tiny subset of the forward reads from a 150 bp paired-end transcriptome run on the heart of an adult female red abalone (*Haliotis rufescens*).

Let's make a *de novo* transcriptome with Trinity. First, request some space on the scavenger queue:

`srun -p compute --time=00:30:00 --ntasks-per-node=2 --mem=20gb --pty bash`

Next, we run Trinity. It doesn't need a lot of parameters to run, but you can change many of the defaults as appropriate. Check these out at the Trinity wiki (https://github.com/trinityrnaseq/trinityrnaseq/wiki), or by typing:

`Trinity --help`
(Highlights - most often used settings.)

`Trinity --show_full_usage_info`
(All the gory details.)

Here's a basic Trinity assembly command for single-end reads:

`Trinity --seqType fq --max_memory 20G --single SRR8956770_1_sub_mod.fq --CPU 2`

But, let's make it interesting and customize a bit:

`--min_kmer_cov : default 1` Set higher for faster assembly (suggest 2).

Trinity is chatty ("verbose" in unix terminology) -- it gives you a lot of output and keeps you updated on its progress.

You should have a new directory called trinity_out_dir. Let's see what's in it.

`cd trinity_out_dir`

We have an assembly!

`Trinity.fasta`

How many contigs are in it?

Let's get to know it better, using a handy perl script included with Trinity. Since Trinity includes this script but doesn't have a built-in command to call it with, we need to run it directly with perl. First, we need to find the script....which means figuring out where our conda installation is actually saved on Poseidon.

`echo $CONDA_PREFIX`

Let's take a look:

`ls [YOUR_TRINITY_CONDA_PATH]`

It's not immediately obvious where the Trinity files are actually located, but we need to find the perl file called TrinityStats.pl, so let's look for it in this directory (and all sub-directories):

`find [YOUR_TRINITY_CONDA_PATH] -name "TrinityStats.pl" -print`

Now that we know the full path to this file, we can run it on our shiny new transcriptome:

`perl [YOUR_PATH_TO_TRINITY_FILE]/TrinityStats.pl Trinity.fasta`

Let's take a look. How many genes are there versus transcripts? What's your N50?

Now, we can also see how many reads are represented in our txm:

`bowtie2-build Trinity.fasta Trinity.fasta
bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -U ../SRR8956770_1_sub_mod.fq > Trin_align_stats.txt | samtools view -@10 -Sb -o bowtie2.bam`

Take a look at your new Trin_align_stats.txt file. What proportion of your reads mapped back? If this were real data, would you be happy with it?

We can also quickly map the reads back on the assembly and calculate how well each contig is supported. We'll use salmon (a "wicked-fast" pseuod-aligner) for this. We'll discuss this in more depth in the next class, so for now, just type these commands:

`salmon index -t Trinity.fasta -i [ASSEMBLY_NICKNAME] -k 25`
`salmon quant -i [ASSEMBLY_NICKNAME] -l A -r ../SRR8956770_1_sub_mod.fq --validateMappings -o [ASSEMBLY_NICKNAME]_quant`

You will end up with a new folder with the output name you specified. In it is a file called `quant.sf`. Use command line tools to explore it a little.

EXERCISE (command line fun):
Which contig has the highest TPM? Does it also have the highest read count?
How many contigs have TPM < 2?
