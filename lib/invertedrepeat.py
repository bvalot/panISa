#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import subprocess
import temporalfile as temp

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
        return " ".join([str(self.prime5_pos), self.prime5_seq, "--", self.prime3_seq, str(self.prime3_pos)])

    def isinrange(self, range_ir):
        if self.prime5_pos <= range_ir and self.prime3_pos >= -range_ir :
            return True

    def isbestis(self, other_ir):
        diffa = (self.prime5_pos - self.prime3_pos)
        diffb = (self.prime5_pos - self.prime3_pos)
        if diffa < diffb:
            return False
        elif diffa > diffb:
            return True
        else:
            raise Exception("Implement differentiation with similar IR pos")
        

    def extractIRinfo(self, len_consensus, irdata):
        ##Extract IR information

        rangeIRpos = 15
        len_bestir = rangeIRpos*2

        if len(irdata) > 3:
            for block in range(len(irdata)/4):#a block refers a result set that contains 4 lines
                line = 4*block

                prime5_pos_temp = int(irdata[0+line].split('_')[-2])
                prime5_seq_temp = irdata[1+line].splitlines()[0]
                prime3_pos_temp = int(irdata[2+line].split('_')[-1])-len_consensus-1
                prime3_seq_temp = irdata[3+line].splitlines()[0]

                ##Check IR in range
                if prime5_pos_temp <= rangeIRpos and prime3_pos_temp >= -rangeIRpos :
                    
                    ##Find the best IR
                    if prime5_pos_temp-prime3_pos_temp < len_bestir:
                        len_bestir = prime5_pos_temp-prime3_pos_temp
                        self.prime5_pos = prime5_pos_temp
                        self.prime3_pos = prime3_pos_temp
                        self.prime5_seq = prime5_seq_temp
                        self.prime3_seq = prime3_seq_temp



def searchir(prime5seq, prime3seq):
    """Searched ir from 5 to 3 prime consensus and return InvertRepeat class or None"""
    range_ir_pos = 15

    ##Prepare input for einverted
    intempfile = temp.writefile("\n".join([">",prime5seq+prime3seq]))

    outtempfile = temp.createfile()
    seqtempfile = temp.createfile()

    ##Call einverted from emboss to find IR
    commandir = ['einverted -sequence '+ intempfile.name +' -gap 12 -threshold 15 -match 3 -mismatch -4 -outfile ' \
    + outtempfile.name + ' -outseq '+ seqtempfile.name + ' -auto Y -warning N']
    processir = subprocess.Popen(commandir, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    processir.communicate()

    ##Read data from outseq of einverted
    outseq_data = temp.readfile(seqtempfile)

    ##Close all temp files
    temp.closefile(intempfile)
    temp.closefile(outtempfile)
    temp.closefile(seqtempfile)

    if len(outseq_data) < 3: ##No IR
        return None
    
    invertRepeat = None
    for block in range(len(outseq_data)/4):#a block refers a result set that contains 4 lines
        line = 4*block
        ##create IR
        prime5_pos_t = int(outseq_data[0+line].split('_')[-2])
        prime5_seq_t = outseq_data[1+line].splitlines()[0]
        prime3_pos_t = int(outseq_data[2+line].split('_')[-1])-len(prime5seq+prime3seq)-1
        prime3_seq_t = outseq_data[3+line].splitlines()[0]
        invert = InvertRepeat(prime5_pos_t, prime5_seq_t, prime3_pos_t, prime3_seq_t)

        ##Check IR in range and choice best one
        if invert.isinrange(range_ir_pos):
            if invertRepeat is None or invert.isbestis(invertRepeat):
                invertRepeat = invert
                
    return invertRepeat

if __name__=='__main__':
    import doctest
    doctest.testmod()
