#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL


class Position():
    """Object that store clip read at a genome position"""
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

    def getconsensusstart(self):
        ##TODO
        return ""

    def getconsensusend(self):
        ##TODO
        return ""
    
    def __len__(self):
        return max(len(self.clipstart), len(self.clipend))
        
    def __repr__(self):
    	return self.chrom + " : " + str(self.pos)
    

class Positions():
    """Container of all positions with a clip read"""
    def __init__(self):
        self.posdict = {}

    def getposition(self, chrom, pos):
        return self.posdict.setdefault(chom+"_"+str(pos), Position(chom, pos))

    def filterposition(self, minread):
        """Filter the positions to remove position insufficiant minread"""
        tmpdict = {}
        for item in self.posdict.items():
            if len(item[1]) >= minread:
                tmpdict[item[O]] = item[1]
        self.posdict = tmpdict

    def nextposition(self):
        """Create iterator that return position in the order of chromosome"""
        ##TODO 
        yield None
                
if __name__=='__main__':
    import doctest
    doctest.testmod()
