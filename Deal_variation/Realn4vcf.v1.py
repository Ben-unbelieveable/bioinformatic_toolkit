#!/jdfstj1/B2C_COM_P1/PipeAdmin/02.software/Conda/bin/python
# -*- coding: utf-8 -*- 
# @Time : 2020/6/2
# @Contact: liubo4@genomics.cn
# @Description: 3' aln for Variations 
# @Version: v1.0.0
# @Update:

from Bio import SeqIO
import argparse
import pandas as pd
import numpy as np
import re
import os
import sys

# args
parser = argparse.ArgumentParser(description="Ref Genome like Hg19")
parser.add_argument("-R", "--refgenome", type=str, default="/jdfstj1/B2C_COM_P1/pipeline/oseq/db/alignment/hg19/hg19.fa",help="path of Ref Genome")
parser.add_argument("-I", "--inputFile", type=str, help="InPutFIle")
parser.add_argument("-B", "--BedFile", type=str, help="BedFIle")
parser.add_argument("-O", "--output_path", type=str, help="the path of output, default: current directory")
args = parser.parse_args()

class Record(object):
	'''
	One line information in vcf file
	'''
	def __init__(self, line):
		info = line.split("\t")
		self.CHROM = info[0] 
		self.Start = int(info[1])
		self.End = int(info[1])+1
		self.ID = info[2]
		self.REF = info[3]
		self.ALT = info[4]
		self.QUAL = info[5]
		self.FILTER = info[6]
		self.INFO = info[7]
		self.Other = info[8:]
		if self.REF[0] == self.ALT[0] and len(self.REF)==1:
			self.REF = "."
			self.ALT = self.ALT[1:]
			self.Type = "ins"
			self.Start = int(info[1])+1
		elif self.REF[0] == self.ALT[0] and len(self.ALT)==1:
			self.REF = self.REF[1:]
			self.ALT = "."
			self.Type = "del"
			self.Start = int(info[1])+1
			self.End = self.Start+len(self.REF)
		else:
			self.Type = "delins"
class VCF(object):
	'''
	VCF class, read VCF, write VCF, get VCF information
	'''
	def __init__(self, uncompress_vcf):
		self.header = []
		self.reader = open(uncompress_vcf, 'r')
		self.line = self.reader.readline().strip()
		while self.line.startswith('#'):
			self.header.append(self.line)
			self.line = self.reader.readline().strip()
		if self.line.startswith("chr"):
			self.record = Record(self.line) 
		else:
			self.reader.close()
	def __iter__(self): 
		return self
	def __next__(self): 
		self.line = self.reader.readline().strip()
		if self.line != "":
			self.record = Record(self.line) 
			return self.record
		else:
			self.reader.close()
			#raise StopIteration()
	def reader_close(self):
		self.reader.close()
def HGVSAlnChrSite(index,mut,refGDict):
	if (str(mut.loc[index, "Strand"]) == "-"):
		#print(mut)
		if (str(mut.loc[index, "BI_MutType"]) == "del"):
			Filt_mut = leftAlnDel(index,mut,refGDict)
			#print("leftAlnDel(index,mut,refGDict)")
		else:
			Filt_mut = leftAlnIns(index,mut, refGDict)
			#print("leftAlnIns(index,mut, refGDict)")
	else:
		if (str(mut.loc[index, "BI_MutType"]) == "del"):
			Filt_mut = rightAlnDel(index,mut,refGDict)
			#print("rightAlnDel(index,mut,refGDict)")
		else:
			Filt_mut = rightAlnIns(index,mut, refGDict)
			#print("rightAlnIns(index,mut, refGDict)")
	return(Filt_mut)

def leftAlnDel(index,mut,refGDict):
	Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "Start"]-2:mut.loc[index, "Start"]-1].seq)
	while(Foward_Base.upper() == mut.loc[index,"Ref"][-1].upper() ):
#		print(mut)
		mut.loc[index, "Start"] = mut.loc[index, "Start"]  - 1
		mut.loc[index, "End"] = mut.loc[index, "End"] - 1
		mut.loc[index,"Ref"] = mut.loc[index,"Ref"][-1]+mut.loc[index,"Ref"][:-1]
		Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "Start"]-2:mut.loc[index, "Start"]-1].seq)
	return(mut)

def leftAlnIns(index,mut,refGDict):
	Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "Start"]-2:mut.loc[index, "Start"]-1].seq)
	while(Foward_Base.upper() == mut.loc[index,"Alt"][-1].upper() ):
#		print(mut)
		mut.loc[index, "Start"] = mut.loc[index, "Start"]  - 1
		mut.loc[index, "End"] = mut.loc[index, "End"] - 1
		mut.loc[index,"Alt"] = mut.loc[index,"Alt"][-1]+mut.loc[index,"Alt"][:-1]
		Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "Start"]-2:mut.loc[index, "Start"]-1].seq)
	return(mut)

def rightAlnDel(index,mut,refGDict):
	Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "End"]-1:mut.loc[index, "End"]].seq)
	while(Foward_Base.upper() == mut.loc[index,"Ref"][0].upper() ):
#		print(mut)
		mut.loc[index, "Start"] = mut.loc[index, "Start"]  + 1
		mut.loc[index, "End"] = mut.loc[index, "End"] + 1
		mut.loc[index,"Ref"] = mut.loc[index,"Ref"][1:]+mut.loc[index,"Ref"][0]
		Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "End"]-1:mut.loc[index, "End"]].seq)
	return(mut)

def rightAlnIns(index,mut,refGDict):
	Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "End"]-1:mut.loc[index, "End"]].seq)
	while (Foward_Base.upper() == mut.loc[index, "Alt"][0].upper()):
		#print(mut)
		mut.loc[index, "Start"] = mut.loc[index, "Start"] + 1
		mut.loc[index, "End"] = mut.loc[index, "End"] + 1
		mut.loc[index, "Alt"] = mut.loc[index,"Alt"][1:]+mut.loc[index,"Alt"][0]
		Foward_Base = str(refGDict[mut.loc[index, "Chr"]][mut.loc[index, "End"]-1:mut.loc[index, "End"]].seq)
	return (mut)

def main():
	BedData = pd.read_csv(args.BedFile, header=0, sep="\t")
	vcf = VCF(args.inputFile)
	O = open(args.output_path,"w")
	O.write("\n".join(vcf.header))
	O.write("\n")
	DataInput = pd.DataFrame(columns = ['Chr','Start','End','Strand','Ref','Alt','BI_MutType'])
	refGDict = SeqIO.to_dict(SeqIO.parse(args.refgenome, "fasta"))
	#DataInput=pd.read_table("/jdfstj1/B2C_COM_P1/Clinical_product/pancancer/20200528/hgvs_aln.test", header=0, sep="\t")
	index=0
	while(vcf.line != ""):
		Chr   = vcf.record.CHROM
		Start = vcf.record.Start
		End   = vcf.record.End
		Ref   = vcf.record.REF
		Alt   = vcf.record.ALT
		QUAL  = vcf.record.QUAL
		FILTER  = vcf.record.FILTER
		MINFO  = vcf.record.INFO
		Type  = vcf.record.Type
		Otherinfo  = "\t".join(vcf.record.Other)
		Strand = "".join(BedData[(BedData["Chr"] == Chr)&(BedData["Start"] <= Start)&(BedData["End"] >= Start)]["Strand"].values)
#		MINFO ="Strand="+Strand+";"+MINFO 
		TMPdata = pd.Series({ "Chr":Chr,"Start":Start,"End":End,"Strand":Strand,"Ref":Ref,"Alt":Alt,"BI_MutType":Type})
#		print(TMPdata)
		DataInput=DataInput.append(TMPdata,1)
#	while(index<DataInput.shape[0]):
#		strand = BedInfo[BedInfo["Chr"] == DataInput.loc[[index],["Chr"]]]
		#print(index)
		if( str(DataInput.loc[index,"BI_MutType"]) == "del"):
			mut = DataInput.loc[[index], ["Chr","Start","End","Strand","Ref","Alt","BI_MutType"]]
			Filt_mut = HGVSAlnChrSite(index,mut,refGDict)
			Start = Filt_mut.loc[index,"Start"]-1
			Alt = str(refGDict[Chr][Start-1:Start].seq).upper()
			Ref = Alt+Filt_mut.loc[index,"Ref"]
			OutLine=Chr+"\t"+str(Start)+"\t.\t"+str(Ref)+"\t"+str(Alt)+"\t"+QUAL+"\t"+FILTER+"\t"+MINFO+"\t"+Otherinfo+"\n"
			O.write(OutLine)
		elif( str(DataInput.loc[index,"BI_MutType"]) == "ins"):
			mut = DataInput.loc[[index], ["Chr","Start","End","Strand","Ref","Alt","BI_MutType"]]
			Filt_mut = HGVSAlnChrSite(index,mut,refGDict)
			Start = Filt_mut.loc[index,"Start"]-1
			Ref = str(refGDict[Chr][Start-1:Start].seq).upper()
			Alt = Ref+Filt_mut.loc[index,"Alt"]
			OutLine=Chr+"\t"+str(Start)+"\t.\t"+str(Ref)+"\t"+str(Alt)+"\t"+QUAL+"\t"+FILTER+"\t"+MINFO+"\t"+Otherinfo+"\n"
			O.write(OutLine)
#			DataInput.loc[[index],["Start"]] = Filt_mut.loc[[index],["Start"]]
#			DataInput.loc[[index],["End"]] = Filt_mut.loc[[index],["End"]]
#			DataInput.loc[[index],["Ref"]] = Filt_mut.loc[[index],["Ref"]]
#			DataInput.loc[[index],["Alt"]] = Filt_mut.loc[[index],["Alt"]]
		else:
			OutLine=Chr+"\t"+str(Start)+"\t.\t"+str(Ref)+"\t"+str(Alt)+"\t"+QUAL+"\t"+FILTER+"\t"+MINFO+"\t"+Otherinfo+"\n"
			O.write(OutLine)
		index+=1
		vcf.__next__()
#	DataInput.to_csv(args.output_path,sep="\t",index=False)
	O.close()
if __name__ == '__main__':
	   main()
