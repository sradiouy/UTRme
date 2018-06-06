# UTRme: a scoring based tool to annotate UTR regions in trypanosomatid genomes 

*The trypanosomatid linage includes a variety of parasitic protozoans causing a significant worldwide problem burden on human health. Given their peculiar mechanisms of gene expression, these organisms depend on post-transcriptional regulation as the main level of gene expression control. Most signals involved in these regulatory networks are located in the untranslated regions (UTRs) of the mRNAs. To deepen our understanding of gene expression regulation we need to identify these regions with high accuracy. Therefore, we have created UTRme (UTR-mini-exon), a GUI stand-alone application to identify and annotate 5’ and 3’ UTR regions. UTRme implements a multiple scoring system to address the issue of false positive UTR assignment that frequently arise because of the characteristics of the available genomes. The tool offers a way for nonbioinformaticians to precisely determine UTRs from transcriptomic data.*

## How to install UTRme?

**If you do not have [conda/miniconda](https://conda.io/miniconda.html) installed, you must first install it**

 1. wget [https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh](https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh)

 1. bash Miniconda3-latest-Linux-x86_64.sh

**Once conda/miniconda is installed, you must install utrme:**

 1. conda config --add channels defaults

 1. conda config --add channels conda-forge

 1. conda config --add channels bioconda

 1. conda install -c sradiouy utrme

**If you previously have a conda/miniconda installed:**

1. conda install -c sradiouy utrme

**If the previous step does not work, you must first configure the conda channels, and then install utrme, as in the instructions above.**

## How to use UTRme?

The software was written in python (version 3) and depends on cutadapt, bedtools, bowtie2, samtools and several python modules which are automatically configured during installation. 

UTRme needs a reference genome (sequence and annotation) (which can be obtained from http://tritrypdb.org/) and raw reads from a RNA-seq experiment. UTRme can use RNA-seq data either sinlge-end or paired-end. 


### Required arguments:

![UTRME1](https://github.com/sradiouy/UTRme/blob/master/utrme1.png)

* FASTQ files location (1)
  * *Folder where the fastq files (gzipped or not) are located.*
    * *Pair 1 or single-end files.*



* FASTQ files location (2)
  * *Folder where the fastq files (gzipped or not) are located.*
    * *Pair 2 or same folder as FASTQ files location (1) if the experiment is single-end.*
  



* Genome
  * *Reference genome in fasta format.* 
    * *Can be obtained from  [tritrypdb!](http://tritrypdb.org/)*
   


* Annotation
  * *Annotation of the genome in gff/gff3 format.*
    * *Can be obtained from  [tritrypdb!](http://tritrypdb.org/)*
  

* Organism
  * *Define the spliced-ledear sequence used in the program.*
    * *T. cruzi*
    * *T. brucei* 
    * *L. major*
  
* Type of experiment
  * *Single-end*
  * *Paired*
  
 
* Basename
  * *Basename of the output files.*
 

### Optional arguments:



The program starts with the removal of adapter sequences and trimming of low-quality ends from reads using cutadapt. By default, UTRme trims the Illumina TrueSeq adapter, but any sequence can be specified by the user. After that, cutadapt is used to identify the reads containing putative poly(A)  (PA) tails or spliced-leader (SL) sequences (secondary regions), allowing for mismatches. By default, we use an error probability 0.01 for poly(A) sequences (adjustable by the user) and one mismatch for SL sequences. In addition, the user can specify the minimum length of the identified secondary region. In order to correctly identify the trans-splicing sites, the user must also select the organism. Currently, Leishmania major, Trypanosoma brucei and Trypanosoma cruzi are available, however the number of available species can easily grow by including more species specific spliced-leader sequences. 

All reads identified by cutadapt are then aligned to the genome using bowtie2 applying the default very-sensitive local end-to-end alignment mode. The subset of reads aligning to intergenic regions  is selected using bedtools and were previously identified by cutadapt. Each putative splice-acceptor site or putative poly(A) addition site is evaluated to assess its reliability. In both cases a final score is calculated that results from adding an individual and global score. 

The individual score indicates the confidence with which each read predicts a given site; it includes three main components: the primary, secondary and accessory scores.  

The primary score reflects a global assessment of the accuracy of the alignment of the sequence that was left after read clipping (primary region) to the genome (genomic primary region), and its calculated based on the evaluation of their similarity. The secondary score evaluates the difference between the putative poly(A) tail or ME sequence of the read (secondary region) and the corresponding genomic region (genomic secondary region). Finally, the accessory score depends on features such as the confidence that the read was not misplaced during mapping, the presence of undetermined nucleotides (Ns), the presence of unannotated open reading frames (ORFs) in the defined region, the presence of splicing conserved signals (AG acceptor and polypyrimidine tract). 

The global score considers the cumulative evidence of all the reads that support a single processing site. The global score for SL sites is proportional to the number of reads that support the site (“occurrences”). For PA sites, in addition to the previous metric, the sequnces of the putative poly(A) tail of all the reads that support a given site is evaluated. 

 Finally, the reported score is calculated by adding the global score to the value of the third quartile of the individual's scores of all the reads that support that particular site (Fig1). The maximum value for this score is set to 100. The higher the score the more confident is the site prediction. All sites with positive scores are reported as they are supported by a reasonable amount of evidence. By default, if a site has a negative score it is not reported, although this can be modified by the user.  

 In summary, the site score recaps many aspects that influence the certainty that a site can be defined with the provided RNA-seq data. A detailed description of all the scores and their calculation is found in the “Score determination” section of the supplementary material. 

