# !/usr/bin/env python2
# -*- coding: utf-8 -*-

import logging
import sys
from argparse import ArgumentParser
from collections import namedtuple
import vcf
import pysam
import re
import pandas as pd
import numpy as np
from multiprocessing import Pool
import subprocess
program = 'Cis_Trans_Check'
version = '1.0.0.2'
update_info = """
-> v1.0.0.2 :  
"""
parser = ArgumentParser(prog=program)
parser.add_argument('-vcf', dest='vcf_file', required=False, action='store', type=str, help='path to input vcf file')
parser.add_argument('-var', dest='filt_var_file', required=False, action='store', type=str, help='path to input vcf file')
parser.add_argument('-bam', dest='bam', required=True, action='store', type=str, help='path to bam')
parser.add_argument('-ref', dest='ref', required=True, action='store', type=str, help='path to reference')
parser.add_argument('-pre', dest='prefix', required=True, action='store', type=str, help='path to output file')
parser.add_argument('-rc', dest='read_cut', required=False, action='store', type=int, default=5, help='Threshold of the reads number support Trans/Cis  ')
parser.add_argument('-fc', dest='freq_cut', required=False, action='store', type=float, default=0.01, help='Threshold of the Trans/Cis Frequence')
parser.add_argument('-cfb', dest='cover_flank_base', required=False, action='store', type=int, default=5, help='flank base num, only the reads cover the \{cover_flank_base\}  Upstream and downstream of the mutation(ins/del) site will be considerd in the Cis/Trans Check and frequence calculation')
parser.add_argument('-filter_flag', dest='filter_flag', required=False, action='store', default="all", type=str,
                    help='path to input vcf file')
parser.add_argument('-nd', dest='MaxnNucleicDistance', required=False, action='store', default=3, type=int,
                    help='Max base distance for each mutation')
parser.add_argument('-ad', dest='MaxproteinDistance', required=False, action='store', default=0, type=int,
                    help='Max base distance for each mutation')
parser.add_argument('-version', action='version', version='%(prog)s ' + version)
args = parser.parse_args()

assert args.filter_flag in ['all', 'samtools', 'nofilter']
logger = logging.getLogger(program)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)

logger.info(program + ' ' + version)

Read = namedtuple('Read', ['read_name', 'pair', 'strand'])
CisTransPair = namedtuple('CisTransPair', ["Achr", "Asite", "Aref", "Aalt", "Bchr", "Bsite", "Bref", "Balt"])
CisTransSingle = namedtuple('CisTransSingle', ["chr", "site", "ref", "alt"])


class VcfMut():
    def __init__(self):
        self.CHROM = 'chr0'
        self.POS = '0'
        self.REF = 'N'
        self.ALT = 'N'
        self.Mut = ""
        self.TYPE = ""
        self.DP = 0
        self.pHGVS = ""
        self.pHGVSsite = ""
        ''' # 测试用例 SNV
        VcfMut = VcfMut()
        VcfMut.POS = 29704715
        VcfMut.CHROM = "chr2"
        VcfMut.REF = "G"
        VcfMut.ALT = "A"
        VcfMut.Mut = ""
        VcfMut.TYPE = "SNV"
        # 测试用例 INS
        VcfMut = VcfMut()
        VcfMut.POS = 16461518
        VcfMut.CHROM = "chr1"
        VcfMut.REF = "T"
        VcfMut.ALT = "TC"
        VcfMut.Mut = ""
        VcfMut.TYPE = "INS"
        # 测试用例 DEL
        VcfMut = VcfMut()
        VcfMut.POS = 86672819
        VcfMut.CHROM = "chr5"
        VcfMut.REF = "TGTT"
        VcfMut.ALT = "T"
        VcfMut.Mut = ""
        VcfMut.TYPE = "DEL"
        '''




def get_snv_support_reads(VcfMut, bam_file, ref, OutBam, mapq=20, baseq=20, overlaps=True, stepper="all", orphans=True):
    support_reads = []
    cover_reads = []
    start_reads = {}
    EndSite = VcfMut.POS + len(VcfMut.REF)
    for pileup_column in bam_file.pileup(region=str(VcfMut.CHROM) + ':' + str(VcfMut.POS) + '-' + str(VcfMut.POS),mapq=mapq , baseq = baseq,
                                         stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        if pileup_column.nsegments > 0:
            for pileup_read in pileup_column.pileups:
                aln = pileup_read.alignment
                read_name = aln.query_name
                pair = 'pe1' if aln.is_read1 else 'pe2'
                strand = '-' if aln.is_reverse else '+'
                read = Read(read_name, pair, strand)
                if pileup_read.is_del or pileup_read.is_refskip or (aln.flag > 1024) or (aln.mapping_quality < mapq) or \
                        aln.query_qualities[pileup_read.query_position] < baseq:
                    continue
                start_reads[read] = [pileup_read.query_position, aln]
    for pileup_column in bam_file.pileup(region=str(VcfMut.CHROM) + ':' + str(EndSite) + '-' + str(EndSite),
                                         stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        if pileup_column.nsegments > 0:
            for pileup_read in pileup_column.pileups:
                aln = pileup_read.alignment
                read_name = aln.query_name
                pair = 'pe1' if aln.is_read1 else 'pe2'
                strand = '-' if aln.is_reverse else '+'
                read = Read(read_name, pair, strand)
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                if read in start_reads:
                    start_query_position, start_aln = start_reads[read]
                    seq = start_aln.query_sequence[start_query_position:pileup_read.query_position]
                    cover_reads.append(aln)
                    if seq.upper() == VcfMut.ALT.upper():
                        support_reads.append(aln)
                        OutBam.write(aln)
    readID_list = []
    readID2Read = {}
    cover_readID = []
    for aln in cover_reads:
        cover_readID.append(aln.query_name)
    for aln in support_reads:
        readID_list.append(aln.query_name)
        readID2Read[aln.query_name] = aln
    return [readID_list,readID2Read,cover_readID]


def get_ins_support_reads(VcfMut, bam_file, ref, OutBam, mapq=20, baseq=20, overlaps=True, stepper="all", orphans=True):
    support_reads = []
    cover_reads = []
    bam = {}
    EndSite = VcfMut.POS + len(VcfMut.REF)
    CoverStart = VcfMut.POS-coverReadFlank
    CoverEnd = EndSite + coverReadFlank
    insLength=len(VcfMut.ALT)-len(VcfMut.REF)
    for pileup_column in bam_file.pileup(region=str(VcfMut.CHROM) + ':' + str(VcfMut.POS) + '-' + str(VcfMut.POS), mapq=mapq, baseq=baseq, stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        if pileup_column.nsegments > 0:
            for pileup_read in pileup_column.pileups:
                aln = pileup_read.alignment
                bam[aln.query_name] = pileup_read
                if (CoverStart in aln.positions) and (CoverEnd in aln.positions):
                    cover_reads.append(aln)
                    if pileup_read.query_position and aln.cigarstring.find("I") > 0:
                        start = pileup_read.query_position-1
                        altstop = pileup_read.query_position - 1 +len(VcfMut.ALT)
                        refstop = pileup_read.query_position-1 + len(VcfMut.REF)
#                    print(aln.query_name+"\t"+VcfMut.CHROM+":"+str(VcfMut.POS)+"\tref:"+str(start)+"-"+str(refstop)+aln.get_reference_sequence()[start:refstop].upper()+">"+aln.query_sequence[start:altstop].upper()+"\talt:"+str(+start)+"-"+str(altstop)+VcfMut.REF.upper()+">"+VcfMut.ALT.upper())
                        if aln.query_sequence[start:altstop].upper() == VcfMut.ALT.upper() and \
                                aln.get_reference_sequence()[start:refstop].upper() == VcfMut.REF.upper():
                            OutBam.write(aln)
                            support_reads.append(aln)
                        elif aln.query_sequence[pileup_read.query_position-insLength:pileup_read.query_position -insLength+ len(VcfMut.ALT)].upper() == VcfMut.ALT.upper() and \
                            aln.get_reference_sequence()[pileup_read.query_position-insLength:pileup_read.query_position - insLength + len(VcfMut.REF)].upper() == VcfMut.REF.upper():
                            OutBam.write(aln)
                            support_reads.append(aln)
                        elif aln.query_sequence[pileup_read.query_position:pileup_read.query_position + len(VcfMut.ALT)].upper() == VcfMut.ALT.upper() and \
                            aln.get_reference_sequence()[pileup_read.query_position:pileup_read.query_position + len(VcfMut.REF)].upper() == VcfMut.REF.upper():
                            OutBam.write(aln)
                            support_reads.append(aln)
                    #else:
                        #print("Deal with Mutation : " + VcfMut.CHROM + ":" + str(VcfMut.POS) + "." + VcfMut.REF + "->" + VcfMut.ALT)
                        """print(aln.query_name + "\t" + aln.query_sequence[
                                                      pileup_read.query_position:pileup_read.query_position + len(
                                                          VcfMut.ALT)].upper() +
                              "\t" + VcfMut.ALT.upper() + "\t##\t" +
                              aln.get_reference_sequence()[
                              pileup_read.query_position:pileup_read.query_position + len(VcfMut.REF)].upper() + "\t" + VcfMut.REF.upper())
                        print(aln.query_name + "\t" + aln.query_sequence[pileup_read.query_position - len(
                            VcfMut.ALT) + 1:pileup_read.query_position + 1].upper() + "\t" + VcfMut.ALT.upper() + "\t##\t" + aln.get_reference_sequence()[
                                                                                                                             pileup_read.query_position - 1:pileup_read.query_position -1 + len(
                                                                                                                                 VcfMut.REF) ].upper() + "\t" + VcfMut.REF.upper())
                                                                                                                                 """
    readID_list = []
    readID2Read = {}
    cover_readID = []
    for aln in cover_reads:
        cover_readID.append(aln.query_name)
    for aln in support_reads:
        readID_list.append(aln.query_name)
        readID2Read[aln.query_name] = aln
    return [readID_list,readID2Read,cover_readID]


def get_del_support_reads(VcfMut, bam_file, ref, OutBam, mapq=20, baseq=20, overlaps=True, stepper="all", orphans=True):
    support_reads = []
    cover_reads = []
    bam = {}
    EndSite = VcfMut.POS + len(VcfMut.REF)
    CoverStart = VcfMut.POS-coverReadFlank
    CoverEnd = EndSite + coverReadFlank
    for pileup_column in bam_file.pileup(region=str(VcfMut.CHROM) + ':' + str(VcfMut.POS) + '-' + str(EndSite), mapq=mapq , baseq = baseq,
                                         stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        if pileup_column.nsegments > 0:
            for pileup_read in pileup_column.pileups:
                aln = pileup_read.alignment
                bam[aln.query_name]=pileup_read
                if (CoverStart in aln.positions) and (CoverEnd in aln.positions):
                    cover_reads.append(aln)
                    if pileup_read.query_position_or_next and aln.cigarstring.find("D") > 0:
                        start = pileup_read.query_position_or_next - 1
                        refstop = pileup_read.query_position_or_next + len(VcfMut.REF) - 1
                        altstop = pileup_read.query_position_or_next +len(VcfMut.ALT) -1
                        #print(aln.query_sequence[start:altstop].upper() + "\t" + VcfMut.ALT.upper() + "\t##\t"+aln.get_reference_sequence()[start:refstop].upper() + "\t" + VcfMut.REF.upper())
                        if aln.get_reference_sequence()[start:refstop].upper() == VcfMut.REF.upper() and aln.query_sequence[start:altstop].upper() == VcfMut.ALT.upper():
                            OutBam.write(aln)
                            support_reads.append(aln)
    readID_list = []
    readID2Read = {}
    cover_readID = []
    for aln in cover_reads:
        cover_readID.append(aln.query_name)
    for aln in support_reads:
        readID_list.append(aln.query_name)
        readID2Read[aln.query_name] = aln
    return [readID_list,readID2Read,cover_readID]


def cis_check(MutA, MutB, bam_pysam,ref, OutBam, CisMutBam):
    MutA_readID_list = []
    MutB_readID_list = []
    coverA_readID = []
    coverB_readID = []
    MutA_readID2Read = {}
    MutB_readID2Read = {}
    Cis_Reads = []
    Cis_reverse_num = 0
    Cis_forward_num = 0
    TransA_Reads = []
    TransA_reverse_num = 0
    TransA_forward_num = 0
    TransB_Reads = []
    TransB_reverse_num = 0
    TransB_forward_num = 0
    if (MutA.TYPE == "SNV"):
        [MutA_readID_list, MutA_readID2Read, coverA_readID] = get_snv_support_reads(MutA, bam_pysam, ref, OutBam)
    elif (MutA.TYPE == "DEL"):
        #print(str(MutA.POS)+":"+MutA.REF+">"+MutA.ALT)
        [MutA_readID_list, MutA_readID2Read, coverA_readID] = get_del_support_reads(MutA, bam_pysam, ref, OutBam)
    elif (MutA.TYPE == "INS"):
        [MutA_readID_list, MutA_readID2Read, coverA_readID] = get_ins_support_reads(MutA, bam_pysam, ref, OutBam)
    if (MutB.TYPE == "SNV"):
        [MutB_readID_list, MutB_readID2Read, coverB_readID] = get_snv_support_reads(MutB, bam_pysam, ref, OutBam)
    elif (MutB.TYPE == "DEL"):
        [MutB_readID_list, MutB_readID2Read, coverB_readID] = get_del_support_reads(MutB, bam_pysam, ref, OutBam)
    elif (MutB.TYPE == "INS"):
        [MutB_readID_list, MutB_readID2Read, coverB_readID] = get_ins_support_reads(MutB, bam_pysam, ref, OutBam)
    #print(str(MutA.POS))
    #print(MutA_readID_list)
    #print(str(MutB.POS))
    #print(MutB_readID_list)
    coverboth_readID = set(coverA_readID) & set(coverB_readID) 
    for Cis_ReadID in (set(MutA_readID_list) & set(MutB_readID_list) & coverboth_readID):
        Cis_Reads.append(Cis_ReadID)
        CisMutBam.write(MutA_readID2Read[Cis_ReadID])
        if MutA_readID2Read[Cis_ReadID].is_reverse:
            Cis_reverse_num += 1
        else:
            Cis_forward_num += 1
    for TransA_ReadID in ((set(MutA_readID_list) - set(MutB_readID_list))&coverboth_readID):
        TransA_Reads.append(TransA_ReadID)
        #CisMutBam.write(MutA_readID2Read[TransA_ReadID])
        if MutA_readID2Read[TransA_ReadID].is_reverse:
            TransA_reverse_num += 1
        else:
            TransA_forward_num += 1
    for TransB_ReadID in ((set(MutB_readID_list) - set(MutA_readID_list))&coverboth_readID):
        TransB_Reads.append(TransB_ReadID)
        #CisMutBam.write(MutA_readID2Read[TransB_ReadID])
        if MutB_readID2Read[TransB_ReadID].is_reverse:
            TransB_reverse_num += 1
        else:
            TransB_forward_num += 1
    return (Cis_Reads, Cis_forward_num, Cis_reverse_num, TransA_Reads , TransA_forward_num, TransA_reverse_num, TransB_Reads, TransB_forward_num, TransB_reverse_num)


def cis_mut_combie(MutA, MutB, reference, bam_pysam, ref, mapq=20, baseq=20, overlaps=True, stepper="all", orphans=True):
    MutA_reads = {}
    forward_readnum = 0
    reverse_readnum = 0
    for pileup_column in bam_pysam.pileup(region=str(MutA.CHROM) + ':' + str(MutA.POS) + '-' + str(MutA.POS),
                                         mapq=mapq, baseq=baseq,
                                         stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        for pileup_read in pileup_column.pileups:
            aln = pileup_read.alignment
            read_name = aln.query_name
            pair = 'pe1' if aln.is_read1 else 'pe2'
            strand = '-' if aln.is_reverse else '+'
            read = Read(read_name, pair, strand)
            MutA_reads[read] = [pileup_read.query_position, aln]
    for pileup_column in bam_pysam.pileup(region=str(MutB.CHROM) + ':' + str(MutB.POS) + '-' + str(MutB.POS),
                                         stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        for pileup_read in pileup_column.pileups:
            aln = pileup_read.alignment
            read_name = aln.query_name
            pair = 'pe1' if aln.is_read1 else 'pe2'
            strand = '-' if aln.is_reverse else '+'
            read = Read(read_name, pair, strand)
            if read in MutA_reads:
                if(read[2] == "+"):
                    forward_readnum += 1
                else:
                    reverse_readnum += 1
    if(MutA.POS > MutB.POS):
        MutA, MutB = MutB, MutA
    cis_mut = VcfMut()
    Mut_Inside_start = MutA.POS+len(MutA.REF)-1
    Mut_Inside_end = MutB.POS - 1
    Mut_Inside_seq = reference.fetch(MutA.CHROM, Mut_Inside_start, Mut_Inside_end).upper()
    cis_mut.CHROM = MutA.CHROM
    cis_mut.POS = MutA.POS
    cis_mut.REF = MutA.REF+Mut_Inside_seq+MutB.REF
    cis_mut.ALT = MutA.ALT + Mut_Inside_seq + MutB.ALT
    cis_mut.TYPE = "DelIns"
    cis_mut.Mut =MutA.Mut
    return(cis_mut, forward_readnum, reverse_readnum)

def single_mut_combie(MutA, bam_pysam, ref, mapq=20, baseq=20, overlaps=True, stepper="all", orphans=True):
    forward_readnum = 0
    reverse_readnum = 0
    for pileup_column in bam_pysam.pileup(region=str(MutA.CHROM) + ':' + str(MutA.POS) + '-' + str(MutA.POS),
                                            mapq=mapq, baseq=baseq, stepper=stepper, fastaFile=ref, max_depth=200000, **{"truncate": True}):
        for pileup_read in pileup_column.pileups:
            aln = pileup_read.alignment
            read_name = aln.query_name
            pair = 'pe1' if aln.is_read1 else 'pe2'
            strand = '-' if aln.is_reverse else '+'
            read = Read(read_name, pair, strand)
            if(read[2] == "+"):
                forward_readnum += 1
            else:
                reverse_readnum += 1

    cis_mut = VcfMut()
    cis_mut.CHROM = MutA.CHROM
    cis_mut.POS = MutA.POS
    cis_mut.REF = MutA.REF
    cis_mut.ALT = MutA.ALT
    cis_mut.TYPE = "DelIns"
    cis_mut.Mut =MutA.Mut
    return(cis_mut, forward_readnum, reverse_readnum)


def PaserVCF(vcfFile, bam_file, ref, outBamfile, CisSupportBamfile, CisSupportresult):
    bam_pysam = pysam.AlignmentFile(bam_file)
    vcf_reader = vcf.Reader(open(vcfFile, "r"))
    PreMut = VcfMut()
    CisMutBam = pysam.AlignmentFile(CisSupportBamfile, "wb", header=bam_pysam.header)
    OutBam = pysam.AlignmentFile(outBamfile, "wb", header=bam_pysam.header)
    CisMutVcf = open(CisSupportresult, "w")
    vcfHead = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=StandardFilter,Description="Depth>=10 && Min_var_freq=0.01 && Min_var_reads=3">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ<25 are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chrM,length=16571,assembly=hg19>
##reference=file:///jdfstj1/B2C_COM_P1/pipeline/pancancer/db/alignment/hg19/hg19.fa
##bcftools_normVersion=1.2+htslib-1.2.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcancer
chr1\t2488039\t.\tT\tA\t.\t.\t.\tGT:AD:AF:DP\t0/1:1538,6:0.004:1546
"""
    CisMutVcf.write(vcfHead)
    hg19 = pysam.FastaFile(ref)
    for vcf_line in vcf_reader:
        # 读取vcf中每条变异结果的信息。
        Mut = VcfMut()
        Mut.POS = vcf_line.POS
        Mut.CHROM = vcf_line.CHROM
        Mut.REF = vcf_line.REF
        Mut.Mut = vcf_line
        #print(PreMut.CHROM +"\t"+Mut.CHROM +"\t"+ str(Mut.POS) +"\t"+ str(PreMut.POS) + "\t"+ str(args.MaxCisRegion))
        if (PreMut.CHROM == Mut.CHROM and abs(Mut.POS-PreMut.POS-len(PreMut.REF)) < args.MaxnNucleicDistance and Mut.POS-PreMut.POS-len(PreMut.REF) >= 0 and len(vcf_line.ALT) <2):
            for ALT in vcf_line.ALT:
                Mut.ALT = ALT.sequence
                #print(Mut.CHROM + ":" + str(Mut.POS) + Mut.REF + ">" + Mut.ALT)
                if (Mut.REF == "." or (len(Mut.ALT)>len(Mut.REF))):
                    Mut.TYPE = "INS"
                elif (Mut.ALT == "." or (len(Mut.REF)>len(Mut.ALT))):
                    Mut.TYPE = "DEL"
                elif (len(Mut.REF) == 1 and len(Mut.ALT) == 1):
                    Mut.TYPE = "SNV"
                #print(str(PreMut.POS) +"\t"+str(Mut.POS)+"\t"+PreMut.TYPE+"\t"+Mut.TYPE+"\n")
                Cis_Support_read, alt_forward_num, alt_reverse_num, TransA_Reads, TransA_forward_num, TransA_reverse_num , TransB_Reads,TransB_forward_num, TransB_reverse_num = cis_check(PreMut, Mut, bam_pysam, ref, OutBam,
                                                                             CisMutBam)
                #print(Cis_Support_read)
                if (len(Cis_Support_read) > Cutoff_4Read):
                    mut_after_combine, read_forwardnum, read_reversenum = cis_mut_combie(PreMut, Mut, hg19, bam_pysam, ref)
                    # mut_after_combine.POS = mut_after_combine.POS -1
                    if mut_after_combine.REF == "":
                        mut_after_combine.REF = "."
                    if mut_after_combine.ALT == "":
                        mut_after_combine.ALT = "."
                    ref_forwardnum = read_forwardnum - alt_forward_num
                    ref_reversenum = read_reversenum - alt_reverse_num
                    Totalref_dnum = ref_forwardnum + ref_reversenum
                    ToralDepth = read_forwardnum + read_reversenum
                    GT_type="0/1"
                    rate = len(Cis_Support_read)/ToralDepth
                    if(rate>0.75):
                         GT_type="1/1"
                    CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                        str(mut_after_combine.POS) + 
                        "\t.\t" + mut_after_combine.REF + "\t" + mut_after_combine.ALT +
                        "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                        str(alt_forward_num) + "," + str(alt_reverse_num) + 
                        ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                        Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                        str(Totalref_dnum) + "," + str(len(Cis_Support_read)) + ":" + str(ToralDepth) + "\n")

                    if (len(TransA_Reads) > Cutoff_4Read):
                        GT_type="0/1"
                        rate = len(TransA_Reads)/ToralDepth
                        if(rate>0.75):
                             GT_type="1/1"
                        CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                            str(PreMut.POS) + 
                            "\t.\t" + PreMut.REF + "\t" + PreMut.ALT +
                            "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                            str(TransA_forward_num) + "," + str(TransA_reverse_num) + 
                            ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                            Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                            str(Totalref_dnum) + "," + str(len(TransA_Reads)) + ":" + str(ToralDepth) + "\n")
                    if (len(TransB_Reads) > Cutoff_4Read):
                        GT_type="0/1"
                        rate = len(TransB_Reads)/ToralDepth
                        if(rate>0.75):
                             GT_type="1/1"
                        CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                            str(Mut.POS) + 
                            "\t.\t" + Mut.REF + "\t" + Mut.ALT +
                            "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                            str(TransB_forward_num) + "," + str(TransB_reverse_num) + 
                            ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                            Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                            str(Totalref_dnum) + "," + str(len(TransB_Reads)) + ":" + str(ToralDepth) + "\n")
                    #print(PreMut.CHROM +"."+ str(PreMut.POS)+":"+ PreMut.REF+">"+PreMut.ALT+"\t"+Mut.CHROM +"."+ str(Mut.POS) + ":"+ Mut.REF+">"+Mut.ALT+ "\t"+ str(args.MaxnNucleicDistance))
                    #print(Cis_Support_read)
                else:
                     pass 
#                    CisMutVcf.write(PreMut.Mut)
            PreMut = Mut
        else:
#            CisMutVcf.write(PreMut.Mut)
            for ALT in vcf_line.ALT:
                Mut.ALT = ALT.sequence
            if (Mut.REF == "." or (len(Mut.ALT)>len(Mut.REF))):
                Mut.TYPE = "INS"
            elif (Mut.ALT == "." or (len(Mut.REF)>len(Mut.ALT))):
                Mut.TYPE = "DEL"
            elif (len(Mut.REF) == 1 and len(Mut.ALT) == 1):
                Mut.TYPE = "SNV"
            PreMut = Mut


def filt_var_paser(filt_var_file, bam_file, ref, outBamfile, CisSupportBamfile, CisSupportresult):
    bam_pysam = pysam.AlignmentFile(bam_file)
    filt_var = open(filt_var_file, "r")
    filt_var_mut_list = []
    CisMutBam = pysam.AlignmentFile(CisSupportBamfile, "wb", header=bam_pysam.header)
    OutBam = pysam.AlignmentFile(outBamfile, "wb", header=bam_pysam.header)
    CisMutVcf = open(CisSupportresult, "w")
    vcfHead = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=StandardFilter,Description="Depth>=10 && Min_var_freq=0.01 && Min_var_reads=3">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ<25 are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=249250621,assembly=hg19>
##contig=<ID=chr2,length=243199373,assembly=hg19>
##contig=<ID=chr3,length=198022430,assembly=hg19>
##contig=<ID=chr4,length=191154276,assembly=hg19>
##contig=<ID=chr5,length=180915260,assembly=hg19>
##contig=<ID=chr6,length=171115067,assembly=hg19>
##contig=<ID=chr7,length=159138663,assembly=hg19>
##contig=<ID=chr8,length=146364022,assembly=hg19>
##contig=<ID=chr9,length=141213431,assembly=hg19>
##contig=<ID=chr10,length=135534747,assembly=hg19>
##contig=<ID=chr11,length=135006516,assembly=hg19>
##contig=<ID=chr12,length=133851895,assembly=hg19>
##contig=<ID=chr13,length=115169878,assembly=hg19>
##contig=<ID=chr14,length=107349540,assembly=hg19>
##contig=<ID=chr15,length=102531392,assembly=hg19>
##contig=<ID=chr16,length=90354753,assembly=hg19>
##contig=<ID=chr17,length=81195210,assembly=hg19>
##contig=<ID=chr18,length=78077248,assembly=hg19>
##contig=<ID=chr19,length=59128983,assembly=hg19>
##contig=<ID=chr20,length=63025520,assembly=hg19>
##contig=<ID=chr21,length=48129895,assembly=hg19>
##contig=<ID=chr22,length=51304566,assembly=hg19>
##contig=<ID=chrX,length=155270560,assembly=hg19>
##contig=<ID=chrY,length=59373566,assembly=hg19>
##contig=<ID=chrM,length=16571,assembly=hg19>
##reference=file:///jdfstj1/B2C_COM_P1/pipeline/pancancer/db/alignment/hg19/hg19.fa
##bcftools_normVersion=1.2+htslib-1.2.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcancer
chr1\t2488039\t.\tT\tA\t.\t.\t.\tGT:AD:AF:DP\t0/1:1538,6:0.004:1546
"""
    CisMutVcf.write(vcfHead)
    hg19 = pysam.FastaFile(ref)
    for mutation in filt_var:
        try:
            mut = mutation.split("\t")
            # 读取vcf中每条变异结果的信息。
            Mut = VcfMut()
            Mut.POS = int(mut[26])+1
            Mut.CHROM = mut[25]
            Mut.REF = mut[28]
            Mut.ALT = mut[29]
            Mut.pHGVS = mut[3]
            if(Mut.pHGVS == "."):
                Mut.pHGVSsite = 0
            else:
                Mut.pHGVSsite = int(re.findall("(\d+)",mut[3])[0])
            Mut.Mut = mut
            if (Mut.REF == "."):
                Mut.REF = ""
                Mut.TYPE = "INS"
            elif (len(Mut.ALT) > len(Mut.REF)):
                Mut.TYPE = "INS"
            elif (Mut.ALT == "."):
                Mut.ALT = ""
                Mut.TYPE = "DEL"
            elif (len(Mut.REF) > len(Mut.ALT)):
                Mut.TYPE = "DEL"
            elif (len(Mut.REF) == 1 and len(Mut.ALT) == 1):
                Mut.TYPE = "SNV"
            filt_var_mut_list.append(Mut)
        except:
            print("Miss Mut:" + mutation)
    for PreMut in filt_var_mut_list:
        for Mut in filt_var_mut_list:
            if(Mut.CHROM != PreMut.CHROM or (Mut.CHROM == PreMut.CHROM and PreMut.POS==Mut.POS and PreMut.REF==Mut.REF and PreMut.ALT==Mut.ALT)):
                pass
            elif((Mut.POS >= PreMut.POS and (Mut.POS+len(Mut.REF))<=(PreMut.POS+len(PreMut.REF)) )):
                    Cis_Support_read, alt_forward_num, alt_reverse_num, TransA_Reads, TransA_forward_num, TransA_reverse_num , TransB_Reads, TransB_forward_num, TransB_reverse_num = \
                    cis_check(PreMut, Mut, bam_pysam, ref, OutBam, CisMutBam)
                    mut_after_combine, read_forwardnum, read_reversenum = single_mut_combie(PreMut, bam_pysam, ref)

                    ref_forwardnum = read_forwardnum - alt_forward_num
                    ref_reversenum = read_reversenum - alt_reverse_num
                    Totalref_dnum = ref_forwardnum + ref_reversenum
                    ToralDepth = read_forwardnum + read_reversenum
                    #PreMut contain Mut, export PreMut 
                    if(Mut.REF==""):
                        Mut.REF = "."
                    if(Mut.ALT ==""):
                        Mut.ALT = "."
                    if (len(Cis_Support_read) > Cutoff_4Read):
                        GT_type="0/1"
                        rate = len(TransA_Reads)/ToralDepth
                        if(rate>0.75):
                             GT_type="1/1"
                        CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                            str(PreMut.POS) + 
                            "\t.\t" + PreMut.REF + "\t" + PreMut.ALT +
                            "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                            str(alt_forward_num) + "," + str(alt_reverse_num) + 
                            ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                            Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                            str(Totalref_dnum) + "," + str(len(Cis_Support_read)) + ":" + str(ToralDepth) + "\n")
                
                        if (len(TransB_Reads) > Cutoff_4Read):
                            GT_type="0/1"
                            rate = len(TransB_Reads)/ToralDepth
                            if(rate>0.75):
                                    GT_type="1/1"
                            CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                                str(Mut.POS) + 
                                "\t.\t" + Mut.REF + "\t" + Mut.ALT +
                                "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                                str(TransB_forward_num) + "," + str(TransB_reverse_num) + 
                                ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                                Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                                str(Totalref_dnum) + "," + str(len(TransB_Reads)) + ":" + str(ToralDepth) + "\n")

            #print(PreMut.CHROM +"\t"+Mut.CHROM +"\t"+ str(Mut.POS) +"\t"+ str(PreMut.POS) + "\t"+ str(PreMut.pHGVSsite)+"\t"+str(Mut.pHGVSsite))
            elif (PreMut.CHROM == Mut.CHROM and (abs(Mut.POS - PreMut.POS - len(PreMut.REF)) < args.MaxnNucleicDistance and
                    Mut.POS > PreMut.POS and Mut.POS - PreMut.POS - len(
                    PreMut.REF) >= 0 and (abs(PreMut.pHGVSsite - Mut.pHGVSsite) <= args.MaxproteinDistance))) :
                    print("Cis Check : " + PreMut.CHROM +"\t"+ str(PreMut.POS)+ "\t"+PreMut.REF+ "\t"+PreMut.ALT + "\t"+Mut.CHROM +"\t"+ str(Mut.POS) + "\t"+Mut.REF+ "\t"+Mut.ALT)
                    # print(Mut.CHROM + ":" + str(Mut.POS) + Mut.REF + ">" + Mut.ALT)
                    Cis_Support_read, alt_forward_num, alt_reverse_num, TransA_Reads, TransA_forward_num, TransA_reverse_num , TransB_Reads, TransB_forward_num, TransB_reverse_num = \
                        cis_check(PreMut, Mut, bam_pysam, ref, OutBam, CisMutBam)
                    if(Mut.REF==""):
                        Mut.REF = "."
                    if(Mut.ALT ==""):
                        Mut.ALT = "."
                    if (len(Cis_Support_read) > Cutoff_4Read):
                        mut_after_combine, read_forwardnum, read_reversenum = cis_mut_combie(PreMut, Mut, hg19, bam_pysam, ref)
                        #mut_after_combine.POS = mut_after_combine.POS -1
                        if mut_after_combine.REF == "":
                            mut_after_combine.REF = "."
                        if mut_after_combine.ALT == "":
                            mut_after_combine.ALT = "."
                        DP = int(mut_after_combine.Mut[17])+int(mut_after_combine.Mut[18])
                        ref_forwardnum = read_forwardnum - alt_forward_num
                        ref_reversenum = read_reversenum - alt_reverse_num
                        Totalref_dnum = ref_forwardnum + ref_reversenum
                        ToralDepth = read_forwardnum + read_reversenum

                        CisMutVcf.write(mut_after_combine.CHROM + "\t" + str(
                            mut_after_combine.POS) + "\t.\t" + mut_after_combine.REF + "\t" + mut_after_combine.ALT +
                                        "\t.\t.\tDP=" + str(ToralDepth) + ";DP4=" + str(ref_forwardnum) + "," + str(
                            ref_reversenum) + "," + str(
                            alt_forward_num) + "," + str(alt_reverse_num) + ";MergeMut=" + PreMut.CHROM + "-" + str(
                            PreMut.POS) + "-" + PreMut.REF + "-" + PreMut.ALT + ":" + Mut.CHROM + "-" + str(
                            Mut.POS) + "-" + Mut.REF + "-" + Mut.ALT + ";\tGT:AD:DP\t1/1:" + str(
                            Totalref_dnum) + "," + str(len(Cis_Support_read)) + ":" + str(ToralDepth) + "\n")
                        if (len(TransA_Reads) > Cutoff_4Read):
                            GT_type="0/1"
                            rate = len(TransA_Reads)/ToralDepth
                            if(rate>0.75):
                                 GT_type="1/1"
                            CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                                str(PreMut.POS) + 
                                "\t.\t" + PreMut.REF + "\t" + PreMut.ALT +
                                "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                                str(TransA_forward_num) + "," + str(TransA_reverse_num) + 
                                ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                                Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                                str(Totalref_dnum) + "," + str(len(TransA_Reads)) + ":" + str(ToralDepth) + "\n")
                        if (len(TransB_Reads) > Cutoff_4Read):
                            GT_type="0/1"
                            rate = len(TransB_Reads)/ToralDepth
                            if(rate>0.75):
                                 GT_type="1/1"
                            CisMutVcf.write(mut_after_combine.CHROM + "\t" + 
                                str(Mut.POS) + 
                                "\t.\t" + Mut.REF + "\t" + Mut.ALT +
                                "\t.\t.\tDP="+str(ToralDepth)+";DP4=" + str(ref_forwardnum) + "," + str(ref_reversenum) + "," + 
                                str(TransB_forward_num) + "," + str(TransB_reverse_num) + 
                                ";MergeMut=" + PreMut.CHROM +"-"+ str(PreMut.POS)+"-"+ PreMut.REF+"-"+PreMut.ALT+":"+ 
                                Mut.CHROM +"-"+ str(Mut.POS) + "-"+ Mut.REF+"-"+Mut.ALT + ";\tGT:AD:DP\t"+ GT_type+":" + 
                                str(Totalref_dnum) + "," + str(len(TransB_Reads)) + ":" + str(ToralDepth) + "\n")



def main():
    AllMutSave = {}
    # Filt_Variants=pd.read_csv(args.variant);
    global Cutoff_4Read  
    global Cutoff_4Freq
    global coverReadFlank
    coverReadFlank = args.cover_flank_base 
    Cutoff_4Read = args.read_cut 
    Cutoff_4Freq = args.freq_cut
    CisSupportBamfile = args.prefix + ".MutSupport.Cis.Reads.bam"
    CisSupportresult = args.prefix + ".vcf"
    outBamfile = args.prefix + ".MutSupport.Raw.Reads.bam"
    if(args.vcf_file ):
        PaserVCF(args.vcf_file, args.bam, args.ref, outBamfile, CisSupportBamfile, CisSupportresult)
    elif(args.filt_var_file):
        filt_var_paser(args.filt_var_file, args.bam, args.ref, outBamfile, CisSupportBamfile, CisSupportresult)

if __name__ == '__main__':
    main()



