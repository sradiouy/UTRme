#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
from collections import Counter, defaultdict
import pandas as pd
import xlsxwriter
import re
import subprocess
import resource 
import numpy as np
import itertools
import time
import sys
import pysam
import regex as re
from time import localtime, strftime
import logging
import matplotlib
matplotlib.use('Agg')
from fuzzywuzzy import fuzz
import glob
import os
import json
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from argparse import ArgumentParser
from gooey import Gooey, GooeyParser
import warnings
from Bio import BiopythonWarning
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import scipy
from matplotlib import rcParams
import gc



@Gooey(optional_cols=2, program_name="UTRme",default_size=(810, 630))

def parse_args():
    """ Use GooeyParser to build up the arguments we will use in UTRme
    Save the arguments in a default json file so that we can retrieve them
    every time we run UTRme.
    """
    stored_args = {}
    # get the script name without the extension & use it to build up
    # the json filename
    script_name = os.path.splitext(os.path.basename(__file__))[0]
    args_file = "{}-args.json".format(script_name)
    if os.path.isfile(args_file):
        with open(args_file) as data_file:
            stored_args = json.load(data_file)
    parser = GooeyParser(description='UTRme: a scoring based tool to annotate UTR regions in trypanosomatid genomes\n\n')
    parser.add_argument('-FASTQ1',
                        required=True,
                        metavar='FASTQ files location',
                        action='store',
                        widget='DirChooser',
                        default=stored_args.get('FASTQ1'),
                        help="Pair 1 (or Single-End)")
    parser.add_argument('-FASTQ2',
                        required=True,
                        metavar='FASTQ files location',
                        action='store',
                        widget='DirChooser',
                        default=stored_args.get('FASTQ2'),
                        help=" Pair 2 (or repeat Single-End location)")
    parser.add_argument('-Genome',
                        required=True,
                        metavar='Genome',
                        action='store',
                        widget='FileChooser',
                        default=stored_args.get('Genome'),
                        help='Fasta format.')
    parser.add_argument('-Annotation',
                        required=True,
                        metavar="Annotation",
                        action='store',
                        widget='FileChooser',
                        default=stored_args.get('Annotation'),
                        help='GFF format.')
    parser.add_argument('-Basename',
                        required=True,
                        metavar = "Basename",
                        action='store',
                        default='Example',
                        help='For output files.')
    parser.add_argument('-i',"--identificator",
                        metavar="id attribute",
                        action='store',
                        help='GFF attribute to be used as feature ID. i.e. ID=TcCLB.506779.120',
                        default="ID=")
    parser.add_argument('-Species',
                        required=True,
                        metavar="Organisim",
                        help="For SL search",
                        action='store',
                        choices=['T. cruzi','T. brucei','L. major'],
                        default = "T. cruzi",
                        widget="Dropdown")
    parser.add_argument("-r", "--rmtemp",
                        metavar = "Remove temp folder",
                        default=False,
                        action="store_true")
    parser.add_argument('-a','--adapter',
                        metavar = "Adapter",
                        action='store',
                        default='AGATCGGAAGAGC',
                        help='Adapter sequences to filter out. If none leave empty.')
    parser.add_argument("-x", "--excel",
                        metavar = "Excel",
                        action="store_true",
                        help="Write output as Excel file")
    parser.add_argument("-f", "--feature",
                        metavar="feature type",
                        action='store',
                        help='Feature type (3rd column in GFF file) to be used, all features of other type are ignored',
                        choices=['gene', 'CDS','polypeptide'],
                        default = "CDS",
                        widget="Dropdown")
    parser.add_argument("-experiment",
                        required=True,
                        metavar="Type of experiment",
                        choices=["Single","Paired"],
                        widget="Dropdown",
                        help="Single or Paired-End",
                        action='store',
                        default = "Paired")
    parser.add_argument("-3UTR", "--polya",
                        metavar = "3'UTR",
                        action="store_false",
                        default=True,
                        help="Performe 3'UTR detection.")
    parser.add_argument("-5UTR", "--splicedleader",
                        metavar = "5'UTR",
                        action="store_false",
                        default=True,
                        help="Performe 5'UTR detection.")
    parser.add_argument("-s", "--score",
                        metavar="Report UTR's with negative score",
                        action="store_true",
                        default=False)    
    parser.add_argument("-uf","--fiveutrlen",
                        metavar="5'UTR length",
                        help="Max. 5'UTRs LENGTH",
                        choices=['500','1000','2000','3000','5000','10000','no filter'],
                        default = "1000",
                        widget="Dropdown")
    parser.add_argument("-ut","--threeutrlen",
                        metavar="3'UTR length",
                        help="Max 3'UTRs LENGTH",
                        choices=['500','1000','2000','3000','5000','10000','no filter'],
                        default = "3000",
                        widget="Dropdown")
    parser.add_argument("-P","--orf",
                        metavar="Max. ORF length (aa) in UTR",
                        help="Do not report UTR with ORFs longer than this value.",
                        choices=['30','50','100','200','300','400','no filter'],
                        default = "200",
                        widget="Dropdown")
    parser.add_argument("-n", "--nbases",
                        metavar="Report UTR's with N's",
                        default=False,
                        action="store_true")
#    parser.add_argument("-m", "--minimum_length",
#                        choices=['30','50','70','90','100','120','150'],
#                        default = "50",
#                        widget="Dropdown",
#                        metavar="Length",
#                        help="Cutadapt option: Discard processed reads that are shorter than this value.")
    parser.add_argument('-o','--overlap',
                        metavar = "Min. overlap length",
                        action='store',
                        help='Cutadapt option: shorter secondary regions are ignored.',
                        choices=['3','4','5','6','7','8','9','10'],
                        default = "5",
                        widget="Dropdown")
    parser.add_argument('-e', '--error_probability',
                        metavar="Max. error rate",
                        action='store',
                        help='Cutadapt option: All searches for secondary regions are error tolerant',
                        choices=['0.05','0.03','0.02','0.01','0.005'],
                        default = "0.01",
                        widget="Dropdown")
    parser.add_argument('-p','--threads',
                        metavar = "Cores",
                        help='Bowtie 2 option: Number of parallel search cores.',
                        action="store",
                        default=4)
    args = parser.parse_args()
    with open(args_file, 'w') as data_file:
        # Using vars(args) returns the data as a dictionary
        json.dump(vars(args), data_file)
    return args



class UTR:
    def __init__(self):
        self.id = ""
        self.chr = ""
        self.start = 0
        self.end = 0
        self.strand = ""
        self.gene = ""
        self.utr_len = 0
        self.intergenic_id = ""
        self.num_errors = 0
        self.srread = ""
        self.srreference = ""
        self.name_adapter = ""
        self.read_start = 0
        self.read_end = 0
        self.read_sequence = ""
        self.reference_sequence = ""
        self.is_multimapped=False
        self.ssite = 0
            
    def hamming_distance(self,s1, s2):
        assert len(s1) == len(s2)
        return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
           
    def score_by_multimapping(self):
        if self.is_multimapped:
            score = -5
        else:
            score = 5
        return score
        
    def score_by_numerrors(self):
        if self.num_errors == 0:
            score = 0
        else:
            score = - self.num_errors
        return score
                    
    def score_AG(self):
        dinucl = self.srreference[-2:]
        AG = -3
        if dinucl == "AG":
            AG = 3
        return dinucl,AG
                
    def primary_score(self):
        remove = False
        fzscore = 0
        fzratio = fuzz.ratio(self.read_sequence,self.reference_sequence)
        if fzratio == 100:
            fzscore = 30
        elif fzratio >= 98:
            fzscore = 25
        elif fzratio >= 95:
            fzscore = 15
        elif fzratio >= 90:
            fzscore = 5
        elif fzratio >= 85:
            fzscore = -15
        elif fzratio >= 75:
            fzscore = -30
        else:
            remove = True
            fzscore = -100
        return fzscore,remove
        
    def secondary_score_sl(self):
        lenScore = 0
        score = 0
        remove = False
        adapter_len = len(self.srread)
        fzratio = fuzz.ratio(self.srread,self.srreference)
        if fzratio >= 95 or "N" in self.srreference:
            remove = True
            return score,remove
        else:
            if adapter_len >= 15:
                lenScore = 5
            score = 15
        hd = float(self.hamming_distance(self.srread,self.srreference))
        hd = (hd/adapter_len)*15 # ~1/6 of the score at maximum is based in the hamming distance
        score += hd + lenScore
        return int(score),remove
        
    def remove_short(self):
        remove = True
        if len(self.srread) >= 4 and not "N" in self.srreference:            
            remove = False
        return remove
        
    def secondary_score_pa(self):
        score = 0
        lenScore = 0
        adapter_len = len(self.srread)
        remove = self.remove_short()
        if remove:
            return score,remove
        else:
            fzratio = fuzz.ratio(self.srread,self.srreference)
            if fzratio >= 95:
                remove = True
                return score,remove
            else:
                if adapter_len >= 15:
                    lenScore = 5
                score = 15 # ~1/3 of the score at maximum is based in the Levenshtein distance
            hd = float(self.hamming_distance(self.srread,self.srreference))
            hd = (hd/adapter_len)*15 # ~1/6 of the score at maximum is based in the Hamming distance
            percentA = (float(self.srreference.count("A"))/adapter_len*10) # percentage of A based correction.
            score += hd + lenScore - percentA
        return int(score),False

def timer(start,end):
    """
    This function calcualtes the elapsed time between the start of the function until the end of this.
    :param: float: starting time.
    :param: float: ending time.
    :return: str: time elapsed.
    """
    hours, rem = divmod(end-start, 3600) #calculate the hours
    minutes, seconds = divmod(rem, 60) #calculate the minutes and seconds
    return ("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

def loginFunction(filename,programname):
    # create logger with 'application'
    logger = logging.getLogger(programname)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(filename)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def remove_folder(path):
    # check if folder exists
    if os.path.exists(path):
         # remove if exists
         shutil.rmtree(path)


def runcml(cml,program,error,logger,bool):
    """Runs a command and logs the output to a specified log handle"""
    startTime = time.time()
    info = "> Running " + program + "!\n"
    logger.info(info)
    logger.info(cml)
    print(info)
    process = subprocess.Popen(cml.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout, stderr = process.communicate()
    stderr = stderr.decode("utf-8")
    stdout = stdout.decode("utf-8")
    if stderr != "":
        mssg = error + ":   " + stderr
        print(cml)
        print(mssg)
        raise Exception(mssg)
    else:
        endTime = time.time()
        elapsedtime =  program + " had a execution time of: " + timer(startTime,endTime)
        logger.info(elapsedtime)
        if bool:
            return stdout

def splice_leader(species):
    sl = {"T. cruzi":"AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG","T. brucei":"AACTAACGCTATTATTAGAACAGTTTCTGTACTATATTG","L. major":"AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG"}
    return sl[species]


def merge_fastq(Fastq1Dir,Fastq2Dir,basename,out_dir,experiment):
    out1 = "path.fastq"
    out2 = "path.fastq"
    path, dirs, files1 = next(os.walk(Fastq1Dir))
    file1_count = len(files1)
    path, dirs, files2 = next(os.walk(Fastq2Dir))
    file2_count = len(files2)
    if file1_count != file2_count:
        error ="The number of fastq files is diferent in", Fastq1Dir, "and", Fastq2Dir
        raise Exception(error)
    elif file1_count == 0 or file2_count == 0:
        error ="One Folder is empty!!!!"
        raise Exception(error)
    else:
        files1 = list((os.path.join(Fastq1Dir,x) for x in files1))
        files2 = list((os.path.join(Fastq2Dir,x) for x in files2))
        files1.sort()
        files2.sort()
        if file1_count == 1:
            return files1[0], files2[0]
        else:
            files1 = list((os.path.join(Fastq1Dir,x) for x in files1))
            files2 = list((os.path.join(Fastq2Dir,x) for x in files2))
            files1.sort()
            files2.sort()
            extension =os.path.splitext(files1[0])[1]
            if extension == ".gz":
                name = basename + "-merged_1.fastq.gz"
            else:
                name = basename + "-merged_1.fastq"
            out1 = os.path.join(out_dir,name)
            with open(out1,'wb') as wfd:
                for fn in files1:
                    with open(fn,'rb') as fd:
                        shutil.copyfileobj(fd, wfd, 1024*1024*10)
            if experiment != "Single":
                if extension != os.path.splitext(files2[0])[1]:
                    raise Exception(error)
                else:
                    if extension == ".gz":
                        name = basename + "-merged_2.fastq.gz"

                    else:
                        name = basename + "-merged_2.fastq"
                    out2 = os.path.join(out_dir,name)
                    with open(out2,'wb') as wfd:
                        for fn in files2:
                            with open(fn,'rb') as fd:
                                shutil.copyfileobj(fd, wfd, 1024*1024*10)
    return out1,out2

def checkmerge(dirfastq1,dirfastq2,experiment):
    check = True
    if experiment == "Single":
        mssg = "You select single-end mode, but your directories are not the same: " + dirfastq1 + " - " + dirfastq2
        check = (dirfastq1 == dirfastq2)
    else:
        check = (dirfastq1 != dirfastq2) 
        mssg = "You select paired-end mode, but your directories are the same: " + dirfastq1 + " - " + dirfastq2
    if not check:        
        raise NameError(mssg)
    return check

def checkadapter(adapter):
    check = True
    if adapter != "":          
        check = bool(re.match('^[ACTG]+$',adapter.upper()))
    if not check:
        raise NameError('Check your Adapter sequence! Must have only ACTG chars.')
    return check

def adapterquality(fastq1,fastq2,adapter,overlap,e,threads,logger,dir,experiment):
    tra = os.path.join(dir,"qr1.fastq.gz")
    if experiment == "Single":
        if adapter != "":
            commandLine = "cutadapt -q 3,3 --quiet --cores " + str(threads) + " -e " + str(e) + " --overlap "+ str(overlap) + " -a " + adapter + "$" + " --quality-base 33 -o " + tra + " " + fastq1
        else:
            commandLine = "cutadapt -q 3,3 --quiet --cores " + str(threads) + " -e " + str(e) + " --overlap "+ str(overlap) + " --quality-base 33 -o " + tra + " " + fastq1
    else:
        trb = os.path.join(dir,"qr2.fastq.gz")
        if adapter != "":
            commandLine = "cutadapt --quiet -e " + str(e) + " --cores " + str(threads) + " --overlap "+ str(overlap) + " -a " + adapter + "$" + " -A " + adapter + "$" + " --quality-base 33 -o " + tra + " -p " + trb + " " + fastq1 + " " + fastq2
        else:
            commandLine = "cutadapt --quiet -e " + str(e) + " --cores " + str(threads) + " --overlap "+ str(overlap) + " --quality-base 33 -o " + tra + " -p " + trb + " " + fastq1 + " " + fastq2
    runcml(commandLine,"cutadapt quality and adapter filtering","cutadapt quality and adapter filtering error",logger,False)
    if experiment == "Single":
        return tra,tra
    else:
        return tra,trb

def keepGeneGFF(annotation,dirout,features,identificador,logger):
    """
    This function return a gff file with only genes features. If the file has sequence inside the program will
    terminate.
    :param: annotation: file : annotation file.
    :param: logger: logger : the logger.
    """
    output_file = os.path.join(dirout,"temp.gff")
    gffgenes = open(output_file,"w")
    with open(annotation, "r") as gff:
        for line in gff:
            if line.startswith(">"):
                sys.exit("Your annotation files contain fasta or perhaps have some lines starting with the character '>'")
            elif line.startswith("#"):
                            pass
            else:
                seqName,source,feature,start,end,score,strand,frame,attribute = line.split("\t")
                if feature == features:
                    gene =  attribute.split(";")[0].replace(identificador,"")
                    try:
                        if ":" in gene:
                            gene =  gene.split(":")[0]
                        elif "-" in gene:
                            gene =  gene.split("-")[0]
                    except:
                        pass
                    gffgenes.write(seqName + "\t" + source + "\t" + feature + "\t" + str(start) + "\t" + str(end) + "\t" + str(score) + "\t"+ strand + "\t" + frame + "\t" + gene + "\n")
                else:
                    pass
    gffgenes.close()
    out = sortGenes(output_file,logger,dirout)
    os.remove(output_file)
    return out

def sortGenes(gff_file,logger,dirout):
    time.sleep(1)
    commandLine = "bedtools sort -i " + gff_file
    stdout = runcml(commandLine,"bedtools sorting","bedtools sort error",logger,True)
    output_file = os.path.join(dirout,"Annot.gff")
    with open(output_file, 'w') as f:
        f.write(stdout)
    return output_file

def chrom_size(genome,dirout,logger):
    commandLine = "faidx " + genome + " -i chromsizes"
    stdout = runcml(commandLine,"faidx","faidx error",logger,True)
    output_file = os.path.join(dirout,"sizes.genome")
    with open(output_file, 'w') as f:
        f.write(stdout)
    return output_file

def sort_chrom_size(genome,gff_file,chrom_size,dirout,logger):
    dfg = pd.read_table(gff_file,sep="\t",header=None)
    dfg.columns = ["seqName","source","feature","start","end","score","strand","frame","gene"]
    dfs = pd.read_table(chrom_size,sep="\t",header=None)
    dfc = pd.DataFrame(list(dfg.seqName.unique()))
    dfc.columns = ["Chrom"]
    dfs.columns = ["Chrom","size"]
    df = pd.merge(dfc,dfs,on="Chrom",how="left")
    output_file = os.path.join(dirout,"sizes.genome")
    df.to_csv(output_file, sep='\t',header=False, index=False)
    return output_file

def get_intergenic(genome,gff_file,chrom_size,dirout,logger):
    commandLine = "bedtools complement -i " +gff_file + " -g " + chrom_size
    stdout = runcml(commandLine,"bedtools complement","bedtools complement error",logger,True)
    
    output_file = os.path.join(dirout,"intergenic.fasta")
    tempfile = os.path.join(dirout,"temp.txt")
    with open(output_file, 'w') as f:
        f.write(stdout)
    clean_intergenic(output_file)
    return output_file

def clean_intergenic(output_file):
    df = pd.read_table(output_file,sep="\t",header=None,names=["id","start","end"])
    df = df[(df["end"] >= 0) & (df["start"] >= 0)]
    df.to_csv(output_file,header=False,index=False,sep="\t")   
    return output_file


def bedtools(genome,annotation,dirout,features,identificador,logger):
    gff_file = keepGeneGFF(annotation,dirout,features,identificador,logger)
    chrom_sizes = chrom_size(genome,dirout,logger)
    sort_chrom_sizes = sort_chrom_size(genome,gff_file,chrom_sizes,dirout,logger)
    intergenic = get_intergenic(genome,gff_file,sort_chrom_sizes,dirout,logger)
    commandLine = "bedtools closest -D ref -a " + intergenic  + " -b " + gff_file
    stdout = runcml(commandLine,"bedtools closest","bedtools closest error",logger,True)
    output_file = os.path.join(dirout,"closest.tab")
    with open(output_file, 'w') as f:
        f.write(stdout)
    os.remove(gff_file)
    os.remove(chrom_sizes)
    os.remove(intergenic)
    return output_file

#START: Detection of Secondary Regions

def detect_secondary_regions(tra,trb,overlap, e, threads,logger,dir,species,side,experiment):
    report = os.path.join(dir,"report")
    if tra.endswith('.gz'):
        tra2 = os.path.join(dir,"sr1.fastq.gz")
        sr = os.path.join(dir,"sr.fastq.gz")
        srfinal = os.path.join(dir,"srfinal.fastq.gz")
    else:
        tra2 = os.path.join(dir,"sr1.fastq")
        sr = os.path.join(dir,"sr.fastq")
        srfinal = os.path.join(dir,"srfinal.fastq.gz")
    if side == 3:
        A = 'A{100}'
        T = 'T{100}'
        e = str(e)
    else:
        T = splice_leader(species)
        A = Seq(T).reverse_complement()._data
        e = "0.03"
    commandLine = "cutadapt --quiet --trimmed-only --minimum-length 30 --no-trim -x P1 --cores " + str(threads) +  " -e " + e + " --overlap "+ str(overlap) + " -a " + A + " -g " + T + " --quality-base 33 -o " + tra2 + " " + tra
    runcml(commandLine,"cutadapt secondary region detection","cutadapt detection error",logger,False)
    if experiment != "Single":
        if tra.endswith('.gz'):
            trb2 = os.path.join(dir,"sr2.fastq.gz")
        else:
            trb2 = os.path.join(dir,"sr2.fastq")
        commandLine = "cutadapt --quiet --trimmed-only --minimum-length 30 --no-trim -x P2 --cores " + str(threads) +  " -e " + e + " --overlap "+ str(overlap) + " -a " + A + " -g " + T + " --quality-base 33 -o " + trb2 + " " + trb
        runcml(commandLine,"cutadapt secondary region detection of second pair","cutadapt detection error",logger,False)
        filenames = [tra2, trb2]                
        with open(sr,'wb') as wfd:
            for fn in filenames:
                with open(fn,'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)
    else:
        sr = tra2
    commandLine = "cutadapt --quiet --trimmed-only --minimum-length 30 --info-file " + report + " -e " + e + " --overlap "+ str(overlap) + " -a " + A + " -g " + T + " --quality-base 33 -o " + srfinal + " " + sr
    runcml(commandLine,"cutadapt reporting","cutadapt reporting error",logger,False)
    return srfinal,report

#END: Detection of Secondary Regions

#START: READ ALINGMENTS

def check_bowtie2_indices(genome):
    """
    check for bowtie2 indices and
    create if needed
    """
    # check index for target species
    genome = os.path.expanduser(genome)
    basename = 0
    genomeName = os.path.splitext(genome)[basename]
    # if index exists, stop here
    if os.path.exists('%s.1.bt2' % genomeName):
        Flag = True
    else:
        Flag = False
    return Flag

def bowtie2(genome,fastq,dirout,basename,threads,side,logger):
    """
    The goal of this function is using bowtie2 for the alignment of the reads.
    First it will check the existence of a bowtie index. Then will align the
    reads to the genome.
    :param genome: file: genome fasta file.
    :param: fastq: file: reads.
    :param basename: file: basename of the alignment file.
    """
    base = basename + "-"+str(side)+".sam"
    output_file = os.path.join(dirout,base)
    genomaBaseName = genome.split(".fasta")[0]
    chk_index = check_bowtie2_indices(genome)
    if not chk_index:
        commandLine = "bowtie2-build --quiet -f " + genome + " " + genomaBaseName
        runcml(commandLine,"bowtie2 build","bowtie2 build error",logger,False)
    commandLine= "bowtie2 --quiet -q --very-sensitive -x " + genomaBaseName + " -S " + output_file + " -U " + fastq + " -p " + str(threads)
    runcml(commandLine,"bowtie2","bowtie2 error",logger,False)
    return output_file


def filter_sort_index(genome,fastq,dirout,basename,threads,side,logger):
    samfile = bowtie2(genome,fastq,dirout,basename,threads,side,logger)
    basename = os.path.splitext(samfile)[0]
    bamfile = basename + ".bam"
    commandLine = "samtools view -F 4 -bh -@ " + str(threads)+  " -o " + bamfile  + " " + samfile
    runcml(commandLine,"samtools view","samtools view error",logger,False)
    os.remove(samfile)
    sortedbam = basename + "_sort" + ".bam"
    out = os.path.join(dirout,"UTRmeSort")
    commandLine = "samtools sort -O bam -m 8G -@ " + str(threads)+  " -T "+ out + " -o  " + sortedbam + " " + bamfile
    runcml(commandLine,"samtools sort","samtools sort error",logger,False)
    commandLine = "samtools index -@ " + str(threads) +" "  + sortedbam
    runcml(commandLine,"samtools index","samtools index error",logger,False)
    os.remove(bamfile)
    return sortedbam

#END: READ ALINGMENTS

#ACCESORY FUNCTION

def ObtainSeq(genome):
    """
    Transforme a genome in fasta format into a dictionary..
    :param: file: genome: genome file in fasta format
    :return: genome sequence dictionary
    """
    handle = open(genome, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict


#START: FILTERING SECONDARY REGIONS

def load_Closest_Table(file,identificator,side):
    my_cols = ["Reference_Name","Start","End","chr","source","feature","start","end","score","strand","frame","gene","distance"]
    use_cols = ["Reference_Name","Start","End","gene","strand","start","end"]
    df = pd.read_table(file,sep="\t",header=None,names=my_cols,usecols=use_cols,dtype={"Reference_Name": 'category', "Start": np.uint32,"End":np.uint32,"gene":'category',"strand":'category',"start":np.uint32,"end":np.uint32})
    try:
        df = df[~(df.gene == ".")]
    except:
        pass
    if side == 3:
        df = df[((df.strand == "+") & (df.end == df.Start)) | ((df.strand == "-") & (df.start - df.End == 1))] #Me quedo con las regiones 3 UTR
    else:
        df = df[((df.strand == "+") & (df.start - df.End == 1)) | ((df.strand == "-") & (df.end == df.Start))] #Me quedo con las regiones 5 UTR
    return df[["Reference_Name","Start","End","gene","strand"]]

def load_Adapter_Report(infofile):
    my_cols =["QueryName","NumErrors","Zero_BasedStartAdapterMatch","Zero_BasedEndAdapterMatch","LeftToAdapter","Adapter","RightToAdapter","NameOfAdapter","QualLeft","QualAdapter","QualRight"]
    use_cols = ["QueryName","NumErrors","Adapter","NameOfAdapter"]
    df1 = pd.read_table(infofile,sep="\t",header=None,names=my_cols,usecols=use_cols,dtype={"QueryName":'object',"NumErrors":np.uint32,"Adapter":'category',"NameOfAdapter":np.uint32})
    df1.QueryName = df1.QueryName.str.split(" ").str.get(0)
    return df1

def mem_usage(pandas_obj):
    if isinstance(pandas_obj,pd.DataFrame):
        usage_b = pandas_obj.memory_usage(deep=True).sum()
    else: # we assume if not a df it's a series
        usage_b = pandas_obj.memory_usage(deep=True)
    usage_mb = usage_b / 1024 ** 2 # convert bytes to megabytes
    return "{:03.2f} MB".format(usage_mb)


def opt_int_df(df):
    df_int = df.select_dtypes(include=['int'])
    df_int = df_int.apply(pd.to_numeric,downcast='unsigned')
    df[df_int.columns] = df_int
    return df


def opt_object_df(df):
    gl_obj = df.select_dtypes(include=['object']).copy()
    converted_obj = pd.DataFrame()
    for col in gl_obj.columns:
        num_unique_values = len(gl_obj[col].unique())
        num_total_values = len(gl_obj[col])
        if num_unique_values / num_total_values < 0.5:
            converted_obj.loc[:,col] = gl_obj[col].astype('category')
        else:
            converted_obj.loc[:,col] = gl_obj[col]
    df[converted_obj.columns] = converted_obj
    return df


def get_reads(bam,refname,start,end,gene,strand):
    list_variable = [refname,start,end,gene,strand] + [x for x in bam.fetch(refname,start,end)]
    return list_variable


def set_polya_site(strand,adapter,reference,start,end):
    seq = adapter
    ref_seq = reference 
    if strand == "+":
        ssite = end
    else:
        ssite = start
    try:
        difference = [i for i in range(len(seq)) if seq[i] != ref_seq[i]]
    except:
        difference = []
    if len(difference) > 0:
        ref_seq = ref_seq[difference[0]:]
        seq = seq[difference[0]:] 
        if strand == "-":
            ssite -= difference[0]
        else:
            ssite += difference[0]
    return ssite,seq,ref_seq

def get_adapter(genome,strand,chr,start,end,nameadapter,adapter,side):
    len_adapter = len(adapter)
    if side ==5:
        if strand == "+":
            ref_seq = genome[chr][(start-len_adapter):start].seq._data
            ssite = start
        else:
            ref_seq = genome[chr][end:(end+len_adapter)].reverse_complement().seq._data
            ssite = end
        if nameadapter == 1:
            seq = Seq(adapter).reverse_complement()._data
        else:
            seq = adapter
    else:
        if nameadapter == 2:
            seq = Seq(adapter).reverse_complement()._data
        else:
            seq = adapter
        if strand == "+":
            if seq[0] != "A":
                end += 1
                seq = seq[1:]  
                len_adapter=len(seq)
            elif seq[-1] != "A":
                seq = seq[:-1]
                len_adapter=len(seq)          
            ref_seq = genome[chr][end:(end + len_adapter)].seq._data       
            ssite,seq,ref_seq = set_polya_site(strand,seq,ref_seq,start,end)
        else:
            if seq[0] != "A":
                start -= 1
                seq = seq[1:]
                len_adapter= len(seq)
            elif seq[-1] != "A":
                seq = seq[:-1]
                len_adapter=len(seq)
            ref_seq = genome[chr][(start-len_adapter):start].reverse_complement().seq._data     
            ssite,seq,ref_seq = set_polya_site(strand,seq,ref_seq,start,end)      
    return seq,ref_seq,ssite


def create_secondary_region_table(closest,report,genomefile,bamfile,identificator,side):
    listr=[]
    dfclosest = load_Closest_Table(closest,identificator,side)
    dfinfo = load_Adapter_Report(report)
    genome = ObtainSeq(genomefile)
    bam = pysam.AlignmentFile(bamfile, "rb")    
    vfunc = np.vectorize(get_reads)
    region = vfunc(bam,dfclosest["Reference_Name"].values, dfclosest["Start"].values, dfclosest["End"].values,dfclosest["gene"].values,dfclosest["strand"].values)
    for r in region:
        chr,start,end,gene,strand = r[:5]
        for read in r[5:]:
            query_seq = read.query_sequence
            real_start = read.reference_start #Start of the alingment in the reference
            real_end = read.reference_end #End of the alingment in the reference
            try:
                if read.get_tag("AS") == read.get_tag("XS"):
                    multimap = True
                else:
                    multimap = False
            except:
                multimap = False
            if strand == "+":
                ref_seq = genome[chr][real_start:real_end].seq._data
            else:
                ref_seq = genome[chr][real_start:real_end].reverse_complement().seq._data
                query_seq = Seq(query_seq).reverse_complement()._data
            lista=[chr,start,end,strand,gene,read.query_name,real_start,real_end,query_seq,ref_seq,read.is_reverse,multimap]
            listr.append(lista)
    columns = ["chr","start","end","strand","gene","QueryName","read_start","read_end","read_sequence","reference_sequence","is_reversed","is_multimapped"]
    del(dfclosest)
    dfread = pd.DataFrame(listr, columns=columns)
    del(listr)
    df = pd.merge(dfread,dfinfo,on=["QueryName"],how="inner")
    del(dfinfo)
    df = opt_int_df(df)
    df.drop("QueryName",axis=1,inplace=True)        
    df = opt_object_df(df)
    if side ==5:
        df = df[((df.is_reversed ==True) & (df.strand == "+") & (df.NameOfAdapter == 1)) | ((df.is_reversed ==False) & (df.strand == "+") & (df.NameOfAdapter == 2)) | ((df.is_reversed ==False) & (df.strand == "-") & (df.NameOfAdapter == 1)) | ((df.is_reversed ==True) & (df.strand == "-") & (df.NameOfAdapter == 2))]
    else:
        df = df[((df.is_reversed ==False) & (df.strand == "+") & (df.NameOfAdapter == 1)) | ((df.is_reversed ==True) & (df.strand == "+") & (df.NameOfAdapter == 2)) | ((df.is_reversed ==True) & (df.strand == "-") & (df.NameOfAdapter == 1)) | ((df.is_reversed ==False) & (df.strand == "-") & (df.NameOfAdapter == 2))]
    vfunc = np.vectorize(get_adapter)
    df["srread"],df["srreference"],df["ssite"] = vfunc(genome,df["strand"],df["chr"],df["read_start"],df["read_end"],df["NameOfAdapter"],df["Adapter"],side)           
    df = opt_int_df(df)
    df = opt_object_df(df)
    return df


#END: FILTERING SECONDARY REGIONS


#START: CREATE SCORE TABLE 

def addscoreinclass(df,genomefile,bamfile,side):
    list_row  = df.values.tolist()
    lists = []
    genome = ObtainSeq(genomefile)
    for item in list_row:
        r = UTR()
        r.chr,r.start,r.end,r.strand,r.gene,r.read_start,r.read_end,r.read_sequence,r.reference_sequence,_,r.is_multimapped,r.num_errors,_,r.name_adapter,r.srread,r.srreference,r.ssite= item
        r.intergenic_id = r.chr  + "_" + str(r.start) + "_" + str(r.end)
        adapter_len = len(r.srread)
        smmap = r.score_by_multimapping()
        snerrors = r.score_by_numerrors()
        saccesory = smmap + snerrors
        sprimary,fzremove = r.primary_score()
        if side == 5:
            AG,sag = r.score_AG()
            saccesory += sag
            ssecondary,remove = r.secondary_score_sl()
            if r.strand == "+":
                r.utr_len = r.end - r.ssite
            else:
                r.utr_len = r.ssite - r.start
            si = sag + smmap + ssecondary + snerrors
            if r.ssite > r.start and r.ssite < r.end:
                pass
            else:
                remove = True
            if not remove and not fzremove:
                si += sprimary
                l = [r.chr,r.start,r.end,r.strand,r.gene,r.utr_len,r.intergenic_id,r.read_start,r.read_end,r.read_sequence,r.reference_sequence,r.is_multimapped,r.num_errors,r.srread,r.srreference,r.ssite,AG,sag,smmap,snerrors,sprimary,ssecondary,saccesory,si]
                lists.append(l)
        else:
            ssecondary,remove = r.secondary_score_pa()            
            if r.strand == "+":
                r.utr_len = r.ssite - r.start
            else:
                r.utr_len = r.end - r.ssite
            si = smmap + ssecondary + snerrors            
            if r.ssite > r.start and r.ssite < r.end:
                pass
            else:
                remove = True
            if not remove and not fzremove:
                si += sprimary
                l = [r.chr,r.start,r.end,r.strand,r.gene,r.utr_len,r.intergenic_id,r.read_start,r.read_end,r.read_sequence,r.reference_sequence,r.is_multimapped,r.num_errors,r.srread,r.srreference,r.ssite,smmap,snerrors,sprimary,ssecondary,saccesory,si]
                lists.append(l)                
    #print('LIST: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if side == 5:
        try:
            columns =["chr","start","end","strand","gene","utr_len","intergenic_id","read_start","read_end","read_sequence","reference_sequence","is_multimapped","num_errors","read_ss_seq","reference_ss_seq","ssite","acceptor","sacceptor","smmap","snerrors","sprimary","ssecondary","saccesory","sread"]
            df = pd.DataFrame(lists,columns=columns)
        except:
            error ="There were no evidence of trans-splicing in your data"
            raise Exception(error)
    else:
        try:
            columns = ["chr","start","end","strand","gene","utr_len","intergenic_id","read_start","read_end","read_sequence","reference_sequence","is_multimapped","num_errors","read_ss_seq","reference_ss_seq","ssite","smmap","snerrors","sprimary","ssecondary","saccesory","sread"]
            df = pd.DataFrame(lists,columns=columns)
        except:
            error ="There were no evidence of polyA in your data"
            raise Exception(error)
    del(l)
    del(lists)
    df = opt_int_df(df)
    df = opt_object_df(df)
    #print('OPT DF: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    #print(mem_usage(df))
    return df


#ADD SCORE 

def get_ORF(seq):
    seq = Seq(seq)
    table = 1
    min_pro_len = 30
    ORF = []
    for nuc in [seq]:
        for frame in range(3):
            for pro in nuc[frame:].translate(table).split("*"):
                if len(pro) >= min_pro_len:
                    if "X" in pro._data:
                        pass
                    else:
                        ORF.append(pro._data)
    try:
        s = max(ORF)
    except:
        s = "-"
    return s, len(s)


def set_occurrence_score(ocurrence,individual_score):
   score = 0
   if individual_score <= 0:       
       return score
   if ocurrence >= 300:
      score = 35
   if ocurrence >= 100:
      score = 30
   if ocurrence >= 50:
      score = 25
   elif ocurrence >= 25:
      score = 20
   elif ocurrence >= 10:
      score = 10
   elif ocurrence >= 2:
      score = 5
   return score

## PA

def getMaxOccuringChar(str):
    # Create array to keep the count of individual characters
    # Initialize the count array to zero
    ASCII_SIZE = 256
    count = [0] * ASCII_SIZE
    # Utility variables
    max = -1
    c = ''
    # Traversing through the string and maintaining the count of
    # each character
    for i in str:
        count[ord(i)]+=1
    ints = [i for i, j in enumerate(count) if j == np.max(count)]
    characters = [chr(n) for n in ints]
    return characters

def counterrors(lists):
    totalreads = len(lists)
    totalerrors = "".join(map(str,lists)).count("T") + "".join(map(str,lists)).count("C") + "".join(map(str,lists)).count("G") +  "".join(map(str,lists)).count("N")
    if totalerrors == 0:
        score = 0
    else:
        score = round(totalerrors/totalreads*10)
    return score

def pileup_score(pileup_pa,error_probability):
    totalreads = pileup_pa.count(",") + 1
    bases = pileup_pa.replace(",","")
    basescount = len(bases)
    terror =  bases.count("T")
    cerror = bases.count("C")
    gerror = bases.count("G")
    nerror = bases.count("N")
    totalerrors = terror + cerror + gerror + nerror
    if totalerrors > 0:
        if totalerrors == terror or  totalerrors == cerror or totalerrors == gerror or totalerrors == nerror:
            if totalreads == totalerrors:
                if totalreads >= 3:
                    score = -500
                    print(pileup_pa,score)
                else:
                    score = -10
            else:
                score = -10
        else:
            score = -round(totalerrors/basescount)*100
            print(pileup_pa,score)
    else:
        if totalreads < 10:
            score = totalreads
        else:
            score = 10
    return score


def pileup_seq(df):
    df["gene"] = df["gene"].astype(object)
    #gdf = df.groupby(["gene","ssite"],as_index=False)["read_ss_seq"].apply(lambda x: x.tolist())
    gdf = df.groupby(["gene","ssite"],as_index=False)["read_ss_seq"].apply(lambda x: ",".join(x))
    gdf = gdf.to_frame(name="pileup_pa")
    gdf.reset_index(level=["gene","ssite"], inplace=True)
    df = pd.merge(df, gdf, on=['gene','ssite'], how='inner')
    #gdf = df.groupby(["gene","ssite"],as_index=False)["reference_ss_seq"].apply(lambda x: max(x.tolist(),key=len))
    #gdf = gdf.to_frame(name="pileup_genome")
    #gdf.reset_index(level=["gene","ssite"], inplace=True)
    #print("2!")
    #print(mem_usage(gdf))
    #df = pd.merge(df, gdf, on=['gene','ssite'], how='inner')
    #print("2merge!")
    df = opt_object_df(df)
    df = opt_int_df(df)
    del(gdf) 
    #print('pup: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return df

#def pileup_score(pileup_pa,pileup_genome):
#    ref_seq = pileup_genome
#    pileup_seq = pileup_pa.split(",")
#    max_len = len(ref_seq)
#    remove = False
#    score = 0
#    mylist=[""]*max_len
#    for item in pileup_seq:
#        for i,y in enumerate(item):
#            mylist[i] = mylist[i] + y #me quedo con todas las letras para cada posicion        
#    for x in range(0,max_len):
#            if len(mylist[x]) > 2: #dos o mas bases para comparar
#                moc = getMaxOccuringChar(mylist[x]) #me quedo con las bases mas abundantes en cada indice.
#                if len(moc) == 1: #Si hay una sola base maxima.
#                    if moc[0] != "A": #Si no es A.                        
#                        if mylist[x].count(moc[0]) == len(mylist[x]):
#                            remove = True #All the bases are equal, so the probability that is due to a sequencing error is low.
#                            continue
#                        elif moc[0] == ref_seq[x]: #Comparation with the references
#                            if mylist[x].count("A") == 0:                          
#                                score -= 2*(mylist[x].count(moc[0])/1)
#                            else:
#                                score -= 2*(mylist[x].count(moc[0])/mylist[x].count("A"))
#                    else:
#                        score += 0 # If the base is in all the xth position in the reads an A.
#                else:
#                    if not "A" in moc:
#                        remove = True #Hay más de una base mayoritaria y ninguna es "A".
#                        continue
#                    else:
#                        score -=  1*(mylist[x].count(moc[0])/mylist[x].count("A"))
#    if remove:
#        score = -500
#        return score
#    if score == 0:
#        score = max_len - counterrors(pileup_pa)
#    return round(score)

def add_pileup_score(df,error_probability):
    df = pileup_seq(df)
    vfunc = np.vectorize(pileup_score)
    #df["pileup_score"] = vfunc(df["pileup_pa"].values,df["pileup_genome"].values)
    df["pileup_score"] = vfunc(df["pileup_pa"].values,error_probability)
    df.drop("pileup_pa",1,inplace = True)
    #df.drop("pileup_genome",1,inplace = True)
    df = df[(df['pileup_score'] !=  -500)]
    df = opt_int_df(df)
    df = opt_object_df(df)   
    return df

def global_score_pa(df,error_probability):
    print("> Adding scores in 3' UTRs!\n")    
    print("> Adding pileup scores in 3' UTRs!\n")
    df = add_pileup_score(df,error_probability)
    print("> Adding pileup scores in 3' UTRs:Done!\n")
    print("> Adding ocurrence and individual scores in 3' UTRs!\n")
    df["occurrence"] = df.groupby(["gene","ssite"])["end"].transform("count")
    df["sprimary"] = df.groupby(["gene","ssite"])["sprimary"].transform(lambda x: x.drop_duplicates().quantile(q=0.75))
    df["ssecondary"] = df.groupby(["gene","ssite"])["ssecondary"].transform(lambda x: x.drop_duplicates().quantile(q=0.75))
    df["saccesory"] = df.groupby(["gene","ssite"])["saccesory"].transform(lambda x: x.drop_duplicates().quantile(q=0.75))
    df["individual_score"] = round(df["sprimary"] + df["ssecondary"] +  df["saccesory"])
    vfunc = np.vectorize(set_occurrence_score)
    df["oscore"] = vfunc(df["occurrence"].values,df["individual_score"].values)    
    df["global_score"] = round(df["oscore"] + df['pileup_score'])
    df["score"] = round(df["global_score"] + df['individual_score'])
    df.drop("oscore",1,inplace=True)
    df.drop("pileup_score",1,inplace=True)
    df.drop("sread",1,inplace=True)
    df = opt_int_df(df)
    df = opt_object_df(df)
    print("> Adding ocurrence and individual scores in 3' UTRs: Done!\n")
    return df




## SL

def global_score_sl(df):
    print("> Adding ocurrence and individual scores in 5' UTRs!\n")
    df["occurrence"] = df.groupby(["gene","ssite"])["end"].transform("count")
    df["sprimary"] = df.groupby(["gene","ssite"])["sprimary"].transform(lambda x: x.drop_duplicates().quantile(q=0.75))
    df["ssecondary"] = df.groupby(["gene","ssite"])["ssecondary"].transform(lambda x: x.drop_duplicates().quantile(q=0.75))
    df["saccesory"] = df.groupby(["gene","ssite"])["saccesory"].transform(lambda x: x.drop_duplicates().quantile(q=0.75))
    df["individual_score"] = round(df["sprimary"] + df["ssecondary"] +  df["saccesory"])
    #print('IS: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)    
    vfunc = np.vectorize(set_occurrence_score)
    df["oscore"] = vfunc(df["occurrence"].values,df["individual_score"].values)    
    df["global_score"] = round(df["oscore"])
    df["score"] = round(df["global_score"] + df['individual_score'])
    df.drop("oscore",axis = 1,inplace=True)
    df.drop("sread",axis = 1,inplace=True)
    df = opt_int_df(df)
    df = opt_object_df(df)
    return df

def groupby_ssite(df,side):
    print("> Grouping by splicing site!\n")    
    if side == 5:        
        df = df.filter(['chr', 'start', 'end', 'strand', 'gene', 'utr_len',   
       'ssite', 'acceptor','sprimary', 'ssecondary',
       'saccesory', 'occurrence', 'individual_score', 'global_score','score'])
    else:
        df = df.filter(['chr', 'start', 'end', 'strand', 'gene', 'utr_len',   
       'ssite','sprimary', 'ssecondary',
       'saccesory', 'occurrence', 'individual_score', 'global_score','score'])
    df["id"] = df['gene'].astype(str) + "-" + df['utr_len'].astype(str)
    df = df.sort_values(['gene','score'],ascending=[False,True]).groupby('id',as_index=False).first()
    df.drop('id', axis = 1, inplace=True)
#    print(mem_usage(df))
    print("> Grouping by splicing site: DONE!\n")
    return df


#polypyrimidine tract


def get_ppt_seq(genome,ssite,strand,chr,start,end):
    if strand == "+":
        seq = genome[chr][start:ssite].seq._data
    else:
        seq = genome[chr][ssite:end].reverse_complement().seq._data
    return seq

def ppt_pattern(genome,ssite,strand,chr,start,end,gscore):
    sequence = get_ppt_seq(genome,ssite,strand,chr,start,end)
    score=0
    sequence = sequence.upper()     
    seq = sequence.replace("T","1")
    seq = seq.replace("C","1")
    seq = seq.replace("A","0")
    seq = seq.replace("G","0")
    seq = seq.replace("N","0")
    if not "11111" in seq:
        gscore -= 10
        return -1,-1,-1,gscore        
    pieces = seq.split("0")
    patt = max('1'.join('0'.join(pieces[i:i+2]) for i in range(0, len(pieces), 2)).split("0"))
    r = '1'.join('0'.join(pieces[i:i+2]) for i in range(0, len(pieces), 2))
    M = re.search(patt,r)
    ppt = sequence[M.start():M.end()]
    pptlen=len(ppt)
    if strand == "+":
        pptdist = len(sequence) - M.end()
    if strand == "-":
        pptdist = len(sequence) - M.start() - pptlen
    if pptlen >= 11:
        score += 1
    else:
        score -= 5
    pptstrech = len(max(re.compile("(T+T)*").findall(ppt.replace("C","T"))))
    if pptstrech >= 15:
        score += 4
    elif pptstrech >= 10:
        score += 2
    elif pptstrech >= 6:
        score += 1
    else:
        score += -5
    if pptdist >= 10 and pptdist <= 40:
        score += 5
    else:
        score -= 5
    sdict={'A':-2,'C':1,'T':2,'G':-2,'N':-5}
    for x in ppt:
        try:
            score += sdict[x]
        except:
            gscore -= 10
            return -1,-1,-1,gscore
    if score > 100:
        score = 100
    score = round(score/10)
    gscore += score
    return ppt,pptdist,pptlen,int(gscore)



def add_ptt_info(genome,df):
    print("> Detection and scoring of PPY tracts!\n")
#    vfunc = np.vectorize(get_ppt_seq)
#    df["pptseq"] = vfunc(genome,df["ssite"],df["strand"],df["chr"],df["start"],df["end"])  
#    print("> Scoring of PPY tracts!\n")
    vfunc = np.vectorize(ppt_pattern)
    df["ppt"],df["ppt_distance"],df["ppt_len"],df["saccesory"] = vfunc(genome,df["ssite"],df["strand"],df["chr"],df["start"],df["end"],df["saccesory"])
    df["individual_score"] = round(df["sprimary"] + df["ssecondary"] +  df["saccesory"])    
    df["score"] = round(df["individual_score"] + df['global_score'])
    df = opt_int_df(df)
    df = opt_object_df(df)
    print("> PPY tracts: DONE!\n")
    return df 

#ORFs and Ns

#UTRs sequence

def get_utr_seq(genome,site,strand,chr,start,end,side):
    if side == 5:
        if strand == "+":
            seq = genome[chr][site:end].seq._data       
        else:
            seq = genome[chr][start:site].reverse_complement().seq._data
    else:
        if strand == "+":
            seq = genome[chr][start:site].seq._data
        else:
            seq = genome[chr][site:end].reverse_complement().seq._data
    return seq


def calc_norf_score(genome,site,strand,chr,start,end,side,score):
    seq = get_utr_seq(genome,site,strand,chr,start,end,side)
    if "N" in seq:
        n = 1
        score += -50
    else:
        n = 0
        score += 0
    s,l = get_ORF(seq)
    if l >= 30:
        orf = s
    else:
        orf = "-"
    if l >= 100:
        score += -50
    elif l >= 50:
        score += -10
    elif l >= 30:
        score += -5    
    return score,n,orf

#def add_seq(genome,df,side):      
#    vfunc = np.vectorize(get_utr_seq)
#    print("> extracting UTRs sequences!\n")
#    print('Pre UTR: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
#    print(mem_usage(df))
#    print(df.info())
#    df["seq"] = vfunc(genome,df["ssite"].values,df["strand"].values,df["chr"].values,df["start"].values,df["end"].values,side)
#    print('Pre UTR: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
#    print(mem_usage(df))
#    print('UTR: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
#    print(mem_usage(df))
#    print(df.info())
#    print("> extracting UTRs sequences: DONE!\n")
#    return df



def add_norf_score(genome,df,side):
    vfunc = np.vectorize(calc_norf_score)
    print("> Computing scores for N and ORFs!\n")
    df["saccesory"],df["n"],df["orf"] = vfunc(genome,df["ssite"].values,df["strand"].values,df["chr"].values,df["start"].values,df["end"].values,side,df["saccesory"].values)    
    df["individual_score"] = round(df["sprimary"] + df["ssecondary"] +  df["saccesory"])    
    df["score"] = round(df["individual_score"] + df['global_score'])
    df = opt_int_df(df)
    df = opt_object_df(df)
    #print(mem_usage(df))
    #print('NORF PPY: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print("> Computing scores for N and ORFs: DONE!\n")
    return df


def final_score(genomefile,report,closest,bamfile,identificator,side,error_probability):
    genome = ObtainSeq(genomefile)             
    if side == 5:
        print("> Filtering out 5' secondary regions!\n")  
        df = create_secondary_region_table(closest,report,genomefile,bamfile,identificator,side)      
        print("> Adding scores in 5'UTRs!\n")
        df = addscoreinclass(df,genomefile,bamfile,side)
        df = global_score_sl(df)
        df = groupby_ssite(df,side)
        df = add_norf_score(genome,df,side)
        df = add_ptt_info(genome,df)
    else:
        print("> Filtering out  3' secondary regions!\n")   
        df = create_secondary_region_table(closest,report,genomefile,bamfile,identificator,side)      
        df = addscoreinclass(df,genomefile,bamfile,side)
        df =global_score_pa(df,error_probability)
        df = groupby_ssite(df,side)
        #df = add_seq(genome,df,side)
        df = add_norf_score(genome,df,side)
    return df


########## GENERATE REPORTS

def convert(x):
    try:
        return x.astype(int)
    except:
        return x    

def gettop(df):
    top ="Yes"
    notop = "-"
    df = df.sort_values(['gene','occurrence','score'],ascending=False)
    listgene = []
    listtop = []
    for row in df.itertuples():
        if row.gene in listgene:
            listtop.append(notop) 
        else:
            listgene.append(row.gene)
            listtop.append(top)     
    df["max_occurrence"] = listtop
    return df

def getbest(df):
    top ="Yes"
    notop = "-"
    df = df.sort_values(['gene','score','occurrence'],ascending=False)
    listgene = []
    listtop = []
    for row in df.itertuples():
        if row.gene in listgene:
            listtop.append(notop) 
        else:
            listgene.append(row.gene)
            listtop.append(top)     
    df["max_score"] = listtop
    return df    

def Generate_Indivuals_Reprots(df,fiveutrlen,threeutrlen,excel,output_directory,basename,score,orf,nbases,side):
    df["gene"] = df["gene"].astype(object)
    if len(df) == 0:
         return False,False,False     
    if not score:
        df = df[df.score.astype(int) > -1]
    if orf != "no filter":
        df = df[df.orf.str.len() <= int(orf)]
    if not nbases:
        df = df[df["n"] == 0]
    if side ==5:
        if fiveutrlen != "no filter":
             df = df[df.utr_len.astype(int) < int(fiveutrlen)]
    else:
        if threeutrlen != "no filter":
             df = df[df.utr_len.astype(int) < int(threeutrlen)]
    if len(df) == 0:
         return False,False,False
    mask = df.score >= 100
    df.loc[mask, "score"] = 100
    df['sites'] = df.groupby("gene")["ssite"].transform("nunique") #todos los sitios.
    df = gettop(df)
    df = getbest(df)
    df = opt_int_df(df)
    df = opt_object_df(df)
    if side==5:
        name = basename + "-5-full-report"
    else:
        name = basename + "-3-full-report"
    output = os.path.join(output_directory,name)
    if side == 5:
        save_df(df[["chr","start","end","strand","gene","utr_len","acceptor","ppt","ppt_distance","ppt_len","sprimary","ssecondary","saccesory","individual_score","global_score","score","occurrence","sites","n","orf","max_occurrence","max_score"]].apply(convert),excel,output)
    else:
        save_df(df[["chr","start","end","strand","gene","utr_len","sprimary","ssecondary","saccesory","individual_score","global_score","score","occurrence","sites","n","orf","max_occurrence","max_score"]].apply(convert),excel,output)    
    df2 = df.sort_values(['gene','score','occurrence',"utr_len"],ascending=[False,False,False,True]).groupby('gene', as_index=False).first()
    if side==5:
        name = basename + "-5-best-summary"
    else:
        name = basename + "-3-best-summary"
    output = os.path.join(output_directory,name)
    if side ==5:        
        save_df(df2[["gene","utr_len","acceptor","score","occurrence","sites"]].apply(convert),excel,output)
    else:
        save_df(df2[["gene","utr_len","score","occurrence","sites"]].apply(convert),excel,output)
    return True,df,df2


def GenerateFasta(df,basename,output_directory,genomefile,side,sites):
    genome = ObtainSeq(genomefile)
    if sites == "Best":
        output = basename + "-" + str(side) + "-best-score.fasta"
    else:
        print("> Creating Fasta File!\n")
        output = basename + "-" + str(side) + "-all.fasta"
    file_location = os.path.join(output_directory,output)
    handle = open(file_location,"w")
    for row in df.itertuples():
        try:
            seq = get_utr_seq(genome,row.ssite,row.strand,row.chr,row.start,row.end,side)
        except:
            print(row)
            continue
        if side==5:
            if row.strand == "+":
                description = " | 5' UTR | Score: " + str(row.score) + " | " + "Ocurrence: " + str(row.occurrence) + " | " + "Length: " + str(len(seq)) + " | Strand: " + row.strand
                id=">"+ row.gene + description
                if len(seq) > 5:
                    handle.write(id + "\n" + seq + "\n")
            else:
                description = " | 5' UTR | Score: " + str(row.score) + " | " + "Ocurrence: " + str(row.occurrence) + " | " + "Length: " + str(len(seq)) + " | Strand: " + row.strand
                id=">"+ row.gene + description
                if len(seq) > 5:
                    handle.write(id + "\n" + seq + "\n")
        else:
            if row.strand == "+":
                description = " | 3' UTR | Score: " + str(row.score) + " | " + "Ocurrence: " + str(row.occurrence) + " | " + "Length: " + str(len(seq)) + " | Strand: " + row.strand
                id=">"+ row.gene + description
                if len(seq) > 5:
                    handle.write(id + "\n" + seq + "\n")
            else:
                description = " | 3' UTR | Score: " + str(row.score) + " | " + "Ocurrence: " + str(row.occurrence) + " | " + "Length: " + str(len(seq)) + " | Strand: " + row.strand
                id=">"+ row.gene + description
                if len(seq) > 5:
                    handle.write(id + "\n" + seq + "\n")
    handle.close()
    return 0


def GenerateGFF(annotation,df,dirout,basename,excel,feature,identificador,logger,side,sites):
    if sites == "Best":
        output = basename + "-" + str(side) + "-best-score.gff"
    else:
        print("> Creating Annotation File!\n")
        output = basename + "-" + str(side) + "-all.gff"
    source = "UTRme"
    frame = "."
    lists = []
    my_cols = ["seqid","source","feature","start","end","score","strand","frame","attributes"]
    dg = pd.read_table(annotation,sep="\t", names=my_cols,comment='#')
    for row in dg.itertuples():
        l = list(row)[1:]
        lists.append(l)
    for row in df.itertuples():
        attributes = "ID=utr_"+row.gene+":mRNA_1;Parent="+row.gene+":mRNA"
        if side == 5:
            if row.strand == "+":
                start = (int(row.ssite) +1)            
                l = [row.chr,source,"five_prime_UTR",start,row.end,row.score,row.strand,frame,attributes]
                lists.append(l)
            else:
                start = (int(row.start) +1)
                l = [row.chr,source,"five_prime_UTR",start,row.ssite,row.score,row.strand,frame,attributes]
                lists.append(l)
        else:
            if row.strand == "+":
                l = [row.chr,source,"three_prime_UTR",row.start,row.ssite,row.score,row.strand,frame,attributes]
                lists.append(l)
            else:
                l = [row.chr,source,"three_prime_UTR",row.ssite,row.end,row.score,row.strand,frame,attributes]
                lists.append(l)
    li = lists
    df = pd.DataFrame(li)
    df.columns = my_cols
    df = df.sort_values(["seqid","start","end"])
    output = os.path.join(dirout,output)
    df.to_csv(output, sep='\t',header=False, index=False)
    return 0

#Graphics
#For each column, first it computes the Z-score of each value in the column, relative to the column mean and standard #deviation. Then is takes the absolute of Z-score because the direction does not matter, only if it is below the #threshold. .all(axis=1) ensures that for each row, all column satisfy the constraint. Finally, result of this condition is #used to index the dataframe.


def utrmekdeplot(df,p,zdesv,label,legend,dirout,output):
    try:
        plt.clf()
    except:
        pass
    rcParams['axes.labelsize'] = 9
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    rcParams['legend.fontsize'] = 9
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Computer Modern Roman']
    rcParams['text.usetex'] = True
    rcParams['figure.figsize'] = 7.3, 4.2
    sns.set(color_codes=True)
    df = df.select_dtypes(['number']) # I keep only numerical columns
    df = df[(np.abs(stats.zscore(df[p])) < zdesv)]
    l = df[p]
    fig = sns.kdeplot(l,  shade=True, color="r");
    x,y = fig.get_lines()[0].get_data()
    cdf = scipy.integrate.cumtrapz(y, x, initial=0)
    nearest_05 = np.abs(cdf-0.5).argmin()
    x_median = x[nearest_05]
    y_median = y[nearest_05]
    mv = " median: " + str(int(x_median))
    plt.vlines(x_median, 0, y_median, colors='gray',linestyles='dashed')
    legend += " - median: " + str(int(df[p].median()))
    if p == "utr_len":
        legend += " nt"
    plt.legend([legend])
    plt.xlabel(label)
    o = output + ".pdf"
    output = os.path.join(dirout,o)
    plt.savefig(output)


def utrmejointplot(df,zdesv,dirout,output):
    try:
        plt.clf()
    except:
        pass
    rcParams['axes.labelsize'] = 9
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    rcParams['legend.fontsize'] = 9
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Computer Modern Roman']
    rcParams['text.usetex'] = True
    rcParams['figure.figsize'] = 7.3, 4.2
    sns.set(color_codes=True)
    df = df.select_dtypes(['number']) # I keep only numerical columns
    df = df[(np.abs(stats.zscore(df["occurrence"])) < zdesv) & (np.abs(stats.zscore(df["score"])) < zdesv)]
    x = df["occurrence"]
    y = df["score"]
    fig = sns.jointplot(x="occurrence", y="score", data=df, kind="hex", space=0, color="r",stat_func=None);
    o = output + ".pdf"
    output = os.path.join(dirout,o)
    plt.savefig(output)


def vsplot(df,df2,p,label,legtxt,zdesv,dirout,output):
    try:
        plt.clf()
    except:
        pass
    rcParams['axes.labelsize'] = 9
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    rcParams['legend.fontsize'] = 9
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Computer Modern Roman']
    rcParams['text.usetex'] = True
    rcParams['figure.figsize'] = 7.3, 4.2
    sns.set(color_codes=True)
    df = df.select_dtypes(['number']) # I keep only numerical columns
    df = df[(np.abs(stats.zscore(df[p])) < zdesv)]
    x = df[p]
    fig = sns.kdeplot(x,  shade=True, color="r");
    mv = str(int(x.median()))
    if p == "utr_len":
        legend1 = "5 " + legtxt + " - median: " + mv + " nt"
    else:
        legend1 = "5 " + legtxt + " - median: " + mv    
    df2 = df2.select_dtypes(['number']) # I keep only numerical columns
    df2 = df2[(np.abs(stats.zscore(df2[p])) < zdesv)]
    y = df2[p]
    fig = sns.kdeplot(y,  shade=True, color="b");
    mv = str(int(y.median()))
    if p == "utr_len":
        legend2 = "3 " + legtxt + " - median: " + mv + " nt"
    else:
        legend2 = "3 " + legtxt + " - median: " + mv        
    plt.xlabel(label)
    plt.legend([legend1,legend2])
    o = output + ".pdf"
    output = os.path.join(dirout,o)
    plt.savefig(output)


###Save DF

def save_df(df,excel,output):
    """ Perform a summary of the data and save the data as an excel file
    """
    if excel:
        output_file = output + ".xlsx"
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='Sheet1',header=True, index=False)
        writer.save()
    else:
        output_file = output + ".tab"
        df.to_csv(output_file, sep='\t',header=True, index=False)


########


def main():
    start_time_program = time.time() # initial time counter
    warnings.simplefilter('ignore', BiopythonWarning)
    startTime = strftime("%Y-%m-%d;%H:%M:%S", localtime()) # for the log file name
    conf = parse_args()
    filename = conf.Basename + "-" + strftime("%Y-%m-%d;%H:%M:%S", localtime()) + ".log"
    programname = "UTRme"
    pd.options.mode.chained_assignment = None
    output_directory = os.getcwd() + "/UTRme_" + strftime("%Y-%m-%d_%H-%M-%S", localtime())
    if not os.path.exists(output_directory):
       os.makedirs(output_directory)
    else:
       error = output_directory + ": Already exist. We need to create one, but we dont want to remove your own folder so, copy to another location or change the name of your old directory"
       raise Exception(error)
    dir = os.path.join(output_directory,"temp")
    os.makedirs(dir)
    annot_dir = os.path.join(output_directory,"GFF")
    os.makedirs(annot_dir)
    report_dir = os.path.join(output_directory,"Reports")
    os.makedirs(report_dir)
    fasta_dir = os.path.join(output_directory,"Fasta")
    os.makedirs(fasta_dir)
    fig_dir = os.path.join(output_directory,"Figures")
    os.makedirs(fig_dir)
    mssg ="\n > WELCOME TO UTRme: UTRme is a scoring based tool to annotate UTR regions in trypanosomatid genomes\n" 
    filename = os.path.join(output_directory,filename)
    logger = loginFunction(filename,programname)
    print(mssg)
    logger.info(mssg)    
    time.sleep(2)
    sp = "> In this case " + conf.Species + "\n"
    print(sp)
    time.sleep(2)
    print("> UTRme is starting! \n")
    time.sleep(2)
    print("> Reads files configuration!\n")
    checkmerge(conf.FASTQ1,conf.FASTQ2,conf.experiment)
    fastq1,fastq2 = merge_fastq(conf.FASTQ1,conf.FASTQ2,conf.Basename,dir,conf.experiment)
    print("> Checking adapter sequence!\n")
    checkadapter(conf.adapter)
    tra,trb = adapterquality(fastq1,fastq2,conf.adapter,conf.overlap,conf.error_probability,conf.threads,logger,dir,conf.experiment)
    closest = bedtools(conf.Genome,conf.Annotation,dir,conf.feature,conf.identificator,logger)
    ##################################################
    #################################################
    #5' UTRs
    if not conf.splicedleader:
        print("> Starting with 5'UTRs...\n")
        time.sleep(2)
        fastq,report = detect_secondary_regions(tra,trb,conf.overlap, conf.error_probability,conf.threads,logger,dir,conf.Species,5,conf.experiment)
        bamfile = filter_sort_index(conf.Genome,fastq,dir,conf.Basename,conf.threads,5,logger)
        #print('Pre Score: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        df = final_score(conf.Genome,report,closest,bamfile,conf.identificator,5,conf.error_probability)
        #print('Score: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        booleansl,dfa,dfbsl = Generate_Indivuals_Reprots(df,conf.fiveutrlen,conf.threeutrlen,conf.excel,report_dir,conf.Basename,conf.score,conf.orf,conf.nbases,5) #Both best score, and all sites.
        del(df)
        #print('Reports: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        mssg = "> Reporting 5'UTRs..."
        logger.info(mssg)
        if not booleansl:
            print("> Can not find 5' UTRs!\n")
        else:
            GenerateFasta(dfa,conf.Basename,fasta_dir,conf.Genome,5,"All")
            GenerateFasta(dfbsl,conf.Basename,fasta_dir,conf.Genome,5,"Best")
            GenerateGFF(conf.Annotation,dfa,annot_dir,conf.Basename,conf.excel,conf.feature,conf.identificator,logger,5,"All")
            del(dfa)
            GenerateGFF(conf.Annotation,dfbsl,annot_dir,conf.Basename,conf.excel,conf.feature,conf.identificator,logger,5,"Best")
            #print('GFF: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            utrmekdeplot(dfbsl,"score",3,"score","score",fig_dir,"5utr-score")
            utrmekdeplot(dfbsl,"sites",3,"sites","sites",fig_dir,"5utr-sites")
            utrmekdeplot(dfbsl,"utr_len",3,"nucleotide","length",fig_dir,"5utr-length")
            utrmejointplot(dfbsl,3,fig_dir,"5utr-score_vs_ocurrence")
            #print('Graphics: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            end_time = time.time()
            fiveInfo = "The analysis of the 5' UTRs has a execution time of: " + timer(start_time_program,end_time) + "."
            logger.info(fiveInfo)    
            print("> The analysis of the 5' UTRs, is finished!\n")
    else:
        booleansl = False
        print("> The analysis of the 5' UTRs, was not selected!\n")
    if not conf.polya:
        print("> Starting with 3'UTRs...\n")
        time.sleep(2)
        fastq,report = detect_secondary_regions(tra,trb,conf.overlap, conf.error_probability,conf.threads,logger,dir,conf.Species,3,conf.experiment)
        bamfile = filter_sort_index(conf.Genome,fastq,dir,conf.Basename,conf.threads,3,logger)
        mssg = "> Scoring 3'UTRs..."
        logger.info(mssg)
        df = final_score(conf.Genome,report,closest,bamfile,conf.identificator,3,conf.error_probability) 
        booleanpa,dfa,dfbpa = Generate_Indivuals_Reprots(df,conf.fiveutrlen,conf.threeutrlen,conf.excel,report_dir,conf.Basename,conf.score,conf.orf,conf.nbases,3) #Both best score, and all sites.
        del(df)
        mssg = "> Reporting 3'UTRs..."
        logger.info(mssg)
        if not booleanpa:
            print("> Can not find 3' UTRs!\n")
        else:
            GenerateFasta(dfa,conf.Basename,fasta_dir,conf.Genome,3,"All")
            GenerateFasta(dfbpa,conf.Basename,fasta_dir,conf.Genome,3,"Best")
            GenerateGFF(conf.Annotation,dfa,annot_dir,conf.Basename,conf.excel,conf.feature,conf.identificator,logger,3,"All")
            del(dfa)
            GenerateGFF(conf.Annotation,dfbpa,annot_dir,conf.Basename,conf.excel,conf.feature,conf.identificator,logger,3,"Best")
            utrmekdeplot(dfbpa,"score",3,"score","score",fig_dir,"3utr-score")
            utrmekdeplot(dfbpa,"utr_len",3,"nucleotide","length",fig_dir,"3utr-length")
            utrmekdeplot(dfbpa,"sites",3,"sites","sites",fig_dir,"3utr-sites")
            utrmejointplot(dfbpa,3,fig_dir,"3utr-score_vs_ocurrence")
            end_time = time.time()
            threeInfo = "> The analysis of the 3' UTRs has a execution time of: " + timer(start_time_program,end_time) + "."
            logger.info(threeInfo) 
            print("> The analysis of the 3' UTRs, is finished!\n")
    else:
        booleanpa = False
        print("> The analysis of the 3' UTRs, was not selected!\n")
    if booleanpa and booleansl:
        vsplot(dfbsl,dfbpa,"score","Score","score",3,fig_dir,"vsplot-score")
        vsplot(dfbsl,dfbpa,"utr_len","UTR Length","length",3,fig_dir,"vsplot-length")
        vsplot(dfbsl,dfbpa,"sites","Gene sites","sites",3,fig_dir,"vsplot-sites")
        del(dfbsl)
        del(dfbpa)    
    #Finishing    
    os.remove(closest)
    if conf.rmtemp:
        if "temp" in dir:
            remove_folder(dir)
    end_time_program = time.time()
    exitProgram= "> The whole program has a execution time of: " + timer(start_time_program,end_time_program) + ".\n"
    logger.info(exitProgram)
    print(exitProgram)
    print('Program: Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print("> Thanks you so much for using UTRme....Good Bye!")
