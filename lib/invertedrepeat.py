#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import subprocess
from . import temporalfile as temp

class InvertRepeat():
    """Class containing inverted repeat at the 5' and the 3'"""
    def __init__(self, prime5_pos, prime5_seq, prime3_pos, prime3_seq):
        ##position beginning ir from clip end consensus (1)
        self.prime5_pos = prime5_pos
        ##Sequence ir from clip end consensus
        self.prime5_seq = prime5_seq
        ##position - len(seq) -1 beginning ir from clip start consensus (-1)
        self.prime3_pos = prime3_pos
        ##sequence ir from clip start consensus (reverse)
        self.prime3_seq = prime3_seq

    def __str__(self):
        return " ".join([str(self.prime5_pos), self.prime5_seq.upper(), "--", \
                         self.prime3_seq.upper(), str(self.prime3_pos)])

    def isinrange(self, range_ir):
        if self.prime5_pos <= range_ir and self.prime3_pos >= -range_ir :
            return True

    def isbestir(self, other_ir):
        diffa = (self.prime5_pos - self.prime3_pos)
        diffb = (other_ir.prime5_pos - other_ir.prime3_pos)
        if diffa < diffb:
            return False
        elif diffa > diffb:
            return True
        else:
            raise Exception("Implement differentiation with similar IR pos")
        

def searchir(prime5seq, prime3seq):
    """Searched ir from 5 to 3 prime consensus and return InvertRepeat class or None"""
    range_ir_pos = 15

    ##test if consensus sequence is not empty
    if prime5seq == "" or prime3seq:
        return None
    
    outseq_data = __performedirsearchwitheinverted(prime5seq, prime3seq)
    len_cons = __lenofconsensus(prime5seq, prime3seq)

    for invert in __getnextirfromeinverted(outseq_data,len_cons):
        if invert is None:
            return None
        invertRepeat = None

        ##Check IR in range and choice best one
        if invert.isinrange(range_ir_pos):
            if invertRepeat is None or invert.isbestir(invertRepeat):
                invertRepeat = str(invert)
                
    return invertRepeat


def __performedirsearchwitheinverted(prime5seq, prime3seq):
    """Take 5 and 3 prime consensus to search invert repeat with einverted"""
    intempfile = temp.writefile("\n".join([">",prime5seq+prime3seq]))

    outtempfile = temp.createfile()
    seqtempfile = temp.createfile()

    ##Call einverted from emboss to find IR
    commandir = ['einverted -sequence '+ intempfile.name +' -gap 12 -threshold 15 -match 3 -mismatch -4 -outfile ' \
    + outtempfile.name + ' -outseq '+ seqtempfile.name + ' -auto Y -warning N']
    processir = subprocess.Popen(commandir, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    processir.communicate()

    ##Read data from outseq of einverted
    outseq_data = temp.readfile(seqtempfile,">")


    ##Close all temp files
    temp.closefile(intempfile)
    temp.closefile(outtempfile)
    temp.closefile(seqtempfile)    

    return outseq_data


def __getnextirfromeinverted(outseq_data,len_cons):
    """read the fasta export of einverted to return invertrepeat object"""
    if len(outseq_data) < 1: ##No IR
        yield None
    
    for block in range(int(len(outseq_data)/2)):#a block refers a result set that contains 2 lists
        line = 2*block
        ##create IR
        prime5_pos_t = int(outseq_data[0+line].split("\n")[0].split('_')[-2])
        prime5_seq_t = "".join(outseq_data[0+line].split("\n")[1:])
        prime3_pos_t = int(outseq_data[1+line].split("\n")[0].split('_')[-1])-len_cons-1
        prime3_seq_t = "".join(outseq_data[1+line].split("\n")[1:])
        invert = InvertRepeat(prime5_pos_t, prime5_seq_t, prime3_pos_t, prime3_seq_t)
        yield invert

def __lenofconsensus(prime5seq, prime3seq):
    return len(prime5seq+prime3seq)

if __name__=='__main__':
    import doctest
    doctest.testmod()
