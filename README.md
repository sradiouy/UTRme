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

![UTRME2](https://github.com/sradiouy/UTRme/blob/master/utrme2.png)

* ID Attribute
  * *GFF attribute to be used as feature ID. i.e. ID=TcCLB.506779.120*
    * *ID=*

* Feature Type


