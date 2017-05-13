import random
import shutil, os
from os import listdir, system
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from util import managefiledir as man

class SimInsertSeq():
    def __init__(self,seqData,len_seqData,amountIS):
        self.seqData = seqData
        self.len_seqData = len_seqData
        self.amountIS = amountIS
        self.totalISinfo = None
        self.genome = None
        self.ISname = None
        self.ISfamily = None
        self.ISgroup = None
        self.len_DR = None
        self.ISseq = None
        self.randomLen_DR = None
        self.randomPos = None
        self.IR = None

    def __randomLen_DR(self):
        if int(len(self.len_DR)) < 1:
            min_DR = 0
            max_DR = 0
        elif "-" in self.len_DR:
            min_DR = int(self.len_DR.split("-")[0])
            max_DR = int(self.len_DR.split("-")[1])
        else:
            min_DR = int(self.len_DR)
            max_DR = int(self.len_DR)
        return random.randint(min_DR,max_DR)

    def __randomPos(self):
        posList = []
        cur_pos_ok = None
        distance = int(self.len_seqData)/1000

        while len(posList) < self.amountIS:
            cur_pos = random.randint(distance,self.len_seqData-distance)
            if len(posList) == 0:
                posList.append(cur_pos)
            elif all(cur_pos >= abs(pos+distance) or cur_pos <= abs(pos-distance) for pos in posList):
                posList.append(cur_pos)
        self.randomPos = sorted(posList)

    def __getISinfo(self,ISinfo_file,i):
        openFile_info = open(ISinfo_file, "r")
        readFile_info = openFile_info.readlines()[1:]
        openFile_info.close()

        self.totalISinfo =  len(readFile_info)
        info = readFile_info[i]

        if len(info) > 0:
            self.genome = info.split("\t")[0]
            self.ISname = info.split("\t")[1]
            self.ISfamily = info.split("\t")[2]
            self.ISgroup = info.split("\t")[3]
            self.len_DR = info.split("\t")[4]
            self.ISseq = info.split("\t")[5]
            self.IR = info.split("\t")[8]
            ##[randomly assign values]
            self.randomLen_DR = self.__randomLen_DR()

    def newSeqRecord(self,ISinfo_file,outISinfo_file):
        newSeq = ''
        nextstart = None
        ##[generate random position once at first time]
        self.__getISinfo(ISinfo_file,0)
        self.__randomPos()

        for i in range(self.amountIS):
            ##[randomly select IS from ISinfo readfile]
            randIS = random.randint(0,self.totalISinfo-1)

            self.__getISinfo(ISinfo_file,randIS)

            if len(newSeq) == 0:##[case 1st time of insertion]
                newSeq = self.seqData.seq[:self.randomPos[i]+self.randomLen_DR]+self.ISseq
                nextstart = self.randomPos[i]
            else:
                newSeq += self.seqData.seq[nextstart:self.randomPos[i]+self.randomLen_DR]+self.ISseq
                nextstart = self.randomPos[i]

            self.__saveISinfo(i,outISinfo_file)

        ##[add termination after inserion]
        newSeq += self.seqData.seq[nextstart:]

        ##[prepare data as record]
        out_record_seq = SeqRecord(newSeq, id = self.genome , description = '')
        # print("Genome: "+self.genome+" Total: "+str(len(newSeq))+" original: "+str(len(self.seqData.seq)))
        return out_record_seq

    def __saveISinfo(self,i,outISinfo_file):
        ISfile = open(outISinfo_file,"a")#+self.genome.split(".")[0],"a")
        header = ["Genome","IS Name","IS Family","IS Group","DR-Length","Random DR-Length",\
        "DR-Seq","IS-Length","IS Position","IR"]
        if i == 0:
            ISfile.write("\t".join(header)+"\n")

        out_info = ([self.genome,self.ISname,self.ISfamily,self.ISgroup,str(self.len_DR),\
            str(self.randomLen_DR),str(self.seqData.seq[self.randomPos[i]:self.randomPos[i]+self.randomLen_DR]),\
        str(len(self.ISseq)),str(self.randomPos[i]+1),self.IR])
        ISfile.write("\t".join(out_info))
        ISfile.close()


def readSeqData(fastaFile,amountIS):
    for read_seqData in SeqIO.parse(fastaFile,"fasta"):
        readData = SimInsertSeq(read_seqData,len(read_seqData.seq),amountIS)
        return readData

def writeSeqRecord(outsimseq_dir, out_record_seq):
    fileName = outsimseq_dir+"simIS.fasta"
    write_file = open(fileName,"w")
    SeqIO.write(out_record_seq, write_file, "fasta")
    write_file.close()
    return fileName


def getSimIS(ISinfo_file,fastaFile,outsimseq_dir,outISinfo_file,amountIS):
    outSeq = readSeqData(fastaFile,amountIS)
    out_record_seq = outSeq.newSeqRecord(ISinfo_file,outISinfo_file)
    out_seq_file = writeSeqRecord(outsimseq_dir,out_record_seq)
    return out_seq_file