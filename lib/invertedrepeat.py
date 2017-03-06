#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import subprocess
import tempfile
import traceback

class InvertRepeat():
    """Class containing inverted repeat at the 5' and the 3'"""
    def __init__(self, prime5_pos, prime5_seq, prime3_pos, prime3_seq):
        self.prime5_pos = prime5_pos ##position beginning ir from clip end consensus (1)
        self.prime5_seq = prime5_seq ##Sequence ir from clip end consensus
        self.prime3_pos = prime3_pos ##position - len(seq) -1 beginning ir from clip start consensus (-1)
        self.prime3_seq = prime3_seq ##sequence ir from clip start consensus (reverse)

    def __str__(self):
        if self.prime5_seq is None:
            return "No IR"
        else :
            return " ".join(self.prime5_pos, self.prime5_seq, self.prime3_pos, self.prime3_seq)

    def __repr__(self):
        if self.prime5_seq is None:
            return "No IR"
        else :
            return str(self.prime5_pos)+" "+self.prime5_seq+" "+str(self.prime3_pos)+" "+self.prime3_seq



def createtempfile():
    newtempfile = tempfile.NamedTemporaryFile()
    return newtempfile

        
def closetempfile(closefile):
    closefile.close()


def writetempfile(writedata):
    temp = tempfile.NamedTemporaryFile()
    try:
        temp.write(writedata)
        temp.seek(0)    
    except Exception:
        exe = traceback.format_exc()
        print exe
    return temp


def readtempfile(tempread):
    readdata = []
    try:
        tempread.seek(0)
        readdata = tempread.readlines()
    except Exception:
        exe = traceback.format_exc()
        print exe
    return readdata



def searchir(prime5seq, prime3seq):
    """Searched ir from 5 to 3 prie consensus and return InvertRepeat class"""
    ##TODO by super Panisa
    ##Only store IR that is close to the 5 and 3 prime (<10bp)

    rangeIRpos = 10
    invertRepeat = []
    prime5_min = rangeIRpos
    prime3_max = -rangeIRpos
    prime5_seq_result = None
    prime3_seq_result = None
    len_seq = len(prime5seq+prime3seq)

    ##Prepare input for einverted
    intempfile = writetempfile('\n'.join(['>',prime5seq+prime3seq]))

    outtempfile = createtempfile()
    seqtempfile = createtempfile()

    ##Call einverted from emboss to find IR
    commandir = ['einverted -sequence '+ intempfile.name +' -gap 12 -threshold 15 -match 3 -mismatch -4 -outfile ' + outtempfile.name + ' -outseq '+ seqtempfile.name + ' -auto Y -warning N']
    processir = subprocess.Popen(commandir, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    processir.communicate()

    ##Read data from outseq of einverted
    outseq_data = readtempfile(seqtempfile)

    ##Extract IR information
    if len(outseq_data) > 3:
        for block in range(len(outseq_data)/4):
            line = 4*block

            prime5_pos = int(outseq_data[0+line].split('_')[-2])
            prime5_seq = outseq_data[1+line].splitlines()[0]
            prime3_pos = int(outseq_data[2+line].split('_')[-1])-len_seq-1
            prime3_seq = outseq_data[3+line].splitlines()[0]

            ##Check IR in range
            if prime5_pos <= rangeIRpos and prime3_pos >= -rangeIRpos :
                
                ##Find the best IR
                if (prime5_pos <= prime5_min) and (prime3_pos >= prime3_max):
                    prime5_min = prime5_pos
                    prime3_max = prime3_pos
                    prime5_seq_result = prime5_seq
                    prime3_seq_result = prime3_seq
            else:
                pass
    else:
        pass

    ##Close all temp files
    closetempfile(intempfile)
    closetempfile(outtempfile)
    closetempfile(seqtempfile)


    invertRepeat.append(InvertRepeat(prime5_min, prime5_seq_result, prime3_max, prime3_seq_result)) 

    return invertRepeat
