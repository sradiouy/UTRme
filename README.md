### Now published in Frontiers in Genetics doi: 10.3389/fgene.2018.00671


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

**The first step to install UTRme is to download the UTRme folder of this site. You can use the git command (if you have it installed) or download the zipped folder and unzip it.**

1. git clone https://github.com/sradiouy/UTRme.git 

**or**

1. unzip UTRme-master.zip

**Then, if you do not have [conda/miniconda](https://conda.io/miniconda.html) installed, you must first install it. Otherwise skip this step.**

 1. wget [https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh](https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh)

 1. bash Miniconda3-latest-Linux-x86_64.sh

**Once conda/miniconda is installed, you must create an enviroment and install UTRme. Preferably in the UTRme folder. The file utrme.yml is inside the UTRme folder. So you must enter to it (cd command).**

 1. conda env create -f utrme.yml

**Then you must activate the previously generated environment. Rembeber to activate it every time you want to use UTRme**

1. source activate UTRme


## How to configure UTRme?

The only external files that UTRme need are, a reference genome (sequence and annotation files, which can be obtained from http://tritrypdb.org/) and raw reads from a RNA-seq experiment. UTRme can use RNA-seq data either sinlge-end or paired-end. 

Also, UTRme requires a configuration file to be executed, in which all the necessary parameters for its execution are established. There are two ways to create the configuration file: through a graphical interface or simply by creating a text file. 

Next, we will describe the parameters used by the program showing the graphical interface. The interface can be viewed in any browser.


### Required arguments:

![UTRME1](https://raw.githubusercontent.com/sradiouy/UTRme/master/Images/utrme_Required.png)

* First pair folder (1)
  * *Folder where the fastq files (gzipped or not) are located.*
    * *Pair 1 or single-end files.*



* Second pair folder (2)
  * *Folder where the fastq files (gzipped or not) are located.*
    * *Pair 2 or leave empty if the experiment is single-end.*
  



* Full path to genome file
  * *Reference genome in fasta format.* 
    * *Can be obtained from  [tritrypdb!](http://tritrypdb.org/)*
   


* Full path to annotation file
  * *Annotation of the genome in gff/gff3 format.*
    * *Can be obtained from  [tritrypdb!](http://tritrypdb.org/)*
  
* Type of experiment
  * *Paired-end*
  * *Single-end*

* Organism
  * *Define the spliced-ledear sequence used in the program.*
    * *T. cruzi*
    * *T. brucei* 
    * *L. major*
    * *Other*
  
  
* Spliced-leader sequence
  * *SL sequence of interest. Remember to put "Other" in Organism*

* Basename
  * *Basename of the output files.*
 

### Optional arguments:

![UTRME2](https://raw.githubusercontent.com/sradiouy/UTRme/master/Images/utrme_Optional.png)

* Feature type
  * *Feature type (3rd column in GFF/GFF3 file) to be used, all features of other type are ignored*
  * *For example:*
     * CDS
     * gene
     * mRNA
     * Polypeptide


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

* Error probability
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


* ID Attribute
  * *GFF attribute to be used as feature ID. i.e. ID=TcCLB.506779.120*
    * *ID=*


* Adapter
  * *Adapter sequences to filter out. If none leave empty.*
    * AGATCGGAAGAGC **default:  Illumina standard adapters**

* Threads
  * *Number of parallel search cores.*


* MultiMapping
  * *YES*  **default**
  * *NO*

* Perform analysis in 3'UTR
  * *YES*  **default**
  * *NO*

* Perform analysis in 5'UTR
  * *YES*  **default**
  * *NO*

* Report UTR's with negative score
  * *YES*
  * *NO*  **default**

* Report UTR's with N's
  * *YES*
  * *NO*  **default**

* Excel
  * *Report results as excel files instead of tabular ones.*
        * *YES*
        * *NO*  **default**
  
* Remove temporary directory
  * *Remove the container folder from temporary files created during the execution of UTRme.*
        * *YES*
        * *NO*  **default**

## How to run the configuration of UTRme?

**python utrme_configuration.py**

*And then enter in a web borwser and to this direction http://127.0.0.1:8050/ .Create the configuration file and exit*

*Or create the configuration file without using the GUI. There is no difference. Please check the example file (Configuration_Files/Example_configuration_file.txt) to know how to create it!*

"""
temporary_value:YES

excel_value:NO

N_value:NO

score_value:NO

utr3_value:YES

utr5_value:YES

mmap_value:YES

core_value:4

adapter_value:AGATCGGAAGAGC

id_value:ID=
orf_value:200
len3_value:3000
len5_value:1000
error_value:0.05
overlap_value:5
feature_value:CDS
basename_value:UTRme-Run
sl_value:
organism_value:T. cruzi
experiment_value:Paired-end
annotation_value:/home/Maria/Reference/annotation.gff
genome_value:/home/Maria/Reference/genome.fasta
second_value:/home/Maria/Experiment/FASTQ2
first_value:/home/Maria/Experiment/FASTQ1"""

## How to contact us?

* **Through the issue section on the UTRme github page**
* **To my personal email: sradio91@gmail.com**

**We hope that UTRme will be useful for your research!** :v::v::v::v:

