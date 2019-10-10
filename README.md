# Lab-Assembly

First, log onto Poseidon and set up a tmux window so we won't be interrupted:

`tmux new -s txm`

Navigate into your user folder for the class, and pull the "Txm-lab" Git repository:

>SET UP

Now, let's set up a conda environment to play in:

```
module load anaconda
conda create --name trinity
conda activate trinity
conda install -c bioconda trinity
```

You should have pulled a folder named "Txm_assembly", which contains a sub-folder called "clean_reads".

This contains one file:

`SRR8956770_1_sub_mod.fq`

This is a tiny subset of the forward reads from a 150 bp paired-end transcriptome run on the heart of an adult female red abalone (*Haliotis rufescens*).

Let's make a *de novo* transcriptome with Trinity!

`Trinity --seqType fq --max_memory 20G --single SRR8956770_1_sub_mod.fq --CPU 2`

Trinity is chatty ("verbose" in unix terminology) -- it gives you a lot of output and keeps you updated on its progress.

You should have a new directory called trinity_out_dir. Let's see what's in it.

`cd trinity_out_dir`

We have an assembly!

`Trinity.fasta`

How many contigs are in it?

Let's get to know it better, using a handy perl script included with Trinity. Since Trinity includes this script but doesn't have a built-in command to call it with, we need to run it directly with perl. First, we need to find the script....which means figuring out where our conda installation is actually saved on Poseidon.

`echo $CONDA_PREFIX`

Let's go take a look:

`cd /vortexfs1/home/ctepolt/.conda/envs/trinity`

It's not immediately obvious where the Trinity files are actually located, but we need to find the perl file called TrinityStats.pl, so let's look for it in this directory (and all sub-directories):

`find . -name "TrinityStats.pl" -print`

Now that we know the full path to this file, we can go back to our working directory and execute it:

`perl /vortexfs1/home/ctepolt/.conda/envs/trinity/opt/TRINITY_HOME/util/TrinityStats.pl Trinity.fasta`

Let's take a look.

Now, we can also see how many reads are represented in our txm:

```
bowtie2-build Trinity.fasta Trinity.fasta
bowtie2 -p 10 -q --no-unal -k 20 -x Trinity.fasta -U ../SRR8956770_1_sub_mod.fq 2 > Trin_align_stats.txt | samtools view -@10 -Sb -o bowtie2.bam
```

Take a look at your new Trin_align_stats.txt file.
