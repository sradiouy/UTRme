# UTRme: a scoring-based tool to annotate untranslated regions in trypanosomatid genomes

**Motivation**:
*Most  signals  involved  in  post-transcriptional  regulatory  networks  are 
located in the untranslated regions (UTRs) of the mRNAs. Therefore, to deepen our 
understanding of gene expression regulation, delimitation of these regions with high 
accuracy  is  needed.  The  trypanosomatid  lineage  includes  a  variety  of  parasitic 
protozoans  causing  a  significant  worldwide  burden  on  human  health.  Given  their 
peculiar  mechanisms  of  gene  expression,  these  organisms  depend  on  post-
transcriptional regulation as the main level of gene expression control. In this context, 
the definition of the UTR regions becomes of key importance.* 

**Results:**
  *We  have  developed  UTRme  (UTR-mini-exon),  a  GUI  (graphical  user 
interface)  stand-alone  application  to  identify  and  annotate  5’  and  3’  UTR  regions. 
UTRme implements a multiple scoring system tailored to address the issue of false 
positive UTR assignment that frequently arise because of the characteristics of the 
intergenic regions. Even though it was developed for trypanosomatids, the tool can be 
used to predict 3’ sites in any eukaryote and it is easily expanded to predict 5’ UTRs in 
any organism where trans-splicing occurs (such as the model organism 
C.  elegans). UTRme  offers  a  way  for  non-bioinformaticians  to  precisely  determine  UTRs  from 
transcriptomic data.*

## How to install UTRme?

**If you do not have [conda/miniconda](https://conda.io/miniconda.html) installed, you must first install it:**

 1. wget [https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh](https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh)

 1. bash Miniconda3-latest-Linux-x86_64.sh

**Once conda/miniconda is installed, you must create an enviroment and install utrme:**

 1. conda create -n utrme python==3.6
 
 1. source activate utrme (remember to use this every time you want to run utrme)

 1. conda config --add channels defaults

 1. conda config --add channels conda-forge

 1. conda config --add channels bioconda

 1. conda install -c sradiouy utrme

**If you previously have a conda/miniconda installed:**

1. conda install -c sradiouy utrme

**If the previous step does not work, you must first configure the conda channels, and then install utrme, as in the instructions above.**

## How to configure UTRme?

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
    * *Other*
  
* Type of experiment
  * *Single-end*
  * *Paired*
  
* Splice-Leader seq
  * *SL sequence of interest. Remember to put "Other" in Organism*

* Basename
  * *Basename of the output files.*
 

### Optional arguments:

![UTRME2](https://github.com/sradiouy/UTRme/blob/master/utrme2.png)

* ID Attribute
  * *GFF attribute to be used as feature ID. i.e. ID=TcCLB.506779.120*
    * *ID=*

* Feature Type
  * *Feature type (3rd column in GFF/GFF3 file) to be used, all features of other type are ignored*
  * *For example:*
     * CDS
     * gene
     * mRNA

* Min. overlap length
  * *Cutadapt option: shorter secondary regions are ignored.*
    * 3 
    * 4 
    * 5      **default**
    * 6 
    * 7 
    * 8 
    * 9 
    * 10
  
* Max. error rate
  * *Cutadapt option: All searches for secondary regions are error tolerant*
    * 0.05
    * 0.03
    * 0.02
    * 0.01      **default**
    * 0.005

* 5'UTR length
  * *Maximum length of the  5 'UTR region.*
    * 500
    * 1000      **default**
    * 2000
    * 3000
    * 5000
    * 10000
    * no filter

* 3'UTR length
  * *Maximum length of the  5 'UTR region.*
    * 500
    * 1000
    * 2000
    * 3000      **default**
    * 5000
    * 10000
    * no filter  

* Max. ORF length (aa) in UTR
  * *Do not report UTR with ORFs longer than this value.*
    * 30
    * 50
    * 100
    * 200      **default**
    * 300
    * 400
    * no filter

* Adapter
  * *Adapter sequences to filter out. If none leave empty.*
    * AGATCGGAAGAGC **default:  Illumina standard adapters**
    
* Multimapping
  * *Remove multimapping reads.*

* Cores
  * *Number of parallel search cores.*

* Remove temporary directory
  * *Remove the container folder from temporary files created during the execution of UTRme.*

* 5'UTR
  * *Performe 5'UTR detection.*

* 3'UTR
  * *Performe 3'UTR detection.*

* Report UTR's with negative score

* Report UTR's with N's

* Excel
  * *Report results as excel files instead of tabular ones.*
  
## How to run UTRme?

**You only need to click on the start button!**

## How to contact us?

* **Through the issue section on the UTRme github page**
* **To my personal email: sradio91@gmail.com**

**We hope that UTRme will be useful for your research!** :v::v::v::v:

