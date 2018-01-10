#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL


class Position():
    """Class that store clip read at a genome position"""
    def __init__(self, chrom, pos):
        self.chrom = chrom
        self.pos = pos
        self.clipstart = []
        self.clipend = []

    def addclipread(self, clipread):
        if clipread.isstartclip():
            self.clipstart.append(clipread)
        else:
            self.clipend.append(clipread)

    def removeclipside(self, minread):
        if len(self.clipstart) < minread:
            self.clipstart = []
        if len(self.clipend) < minread:
            self.clipend = []
            
    def getconsensus(self, percent, start):
        if start:
            seqs = [clip.getclipseq() for clip in self.clipstart]
        else:
            seqs = [clip.getclipseq() for clip in self.clipend]

        while None in seqs:
            seqs.remove(None)
        if len(seqs) == 0:
            return ""
        cons = []
        for i in range(0, max(map(len, seqs))):
            allele = set()
            values = []
            for seq in seqs:
                if len(seq) <= i:
                    continue
                if start:
                    base = seq[-i-1]
                else:
                    base = seq[i]
                allele.add(base)
                values.append(base)
            toadd = "N"
            for base in allele:
                if len(values) >= 2 and \
                   float(values.count(base))/len(values) >= percent:
                    toadd = base
            cons.append(toadd)
        if start:
            return "".join(cons[::-1])
        else:
            return "".join(cons)
    
    def __len__(self):
        return max(len(self.clipstart), len(self.clipend))
        
    def __repr__(self):
    	return self.chrom + " : " + str(self.pos)

    def __str__(self):
        return "\t".join([self.chrom, str(self.pos), str(len(self.clipstart)), \
                          str(len(self.clipend))])
    def __eq__(self, other):
        if self.chrom == other.chrom and self.pos == other.pos:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.pos < other.pos
        else:
            return self.chrom < other.chrom
        

class Positions():
    """Container of all positions with a clip read"""
    def __init__(self):
        self.posdict = {}

    def getposition(self, chrom, pos):
        return self.posdict.setdefault(chrom+"_"+str(pos), Position(chrom, pos))

    def filterposition(self, minread):
        """Filter the positions to remove position insufficiant minread"""
        tmpdict = {}
        for item in self.posdict.items():
            if len(item[1]) >= minread:
                item[1].removeclipside(minread)
                tmpdict[item[0]] = item[1]
        self.posdict = tmpdict

    def nextposition(self):
        """Create iterator that return position in the order of chromosome"""
        posits = list(self.posdict.values())
        posits.sort()
        for pos in posits:
            yield pos
        
                
if __name__=='__main__':
    import doctest
    doctest.testmod()
