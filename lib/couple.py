#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

class Couple():
    """Class containing a start and end clip position with potential IS"""
    def __init__(self, posstart, posend):
        self.posstart = posstart
        self.posend = posend
        self.consstart = None
        self.consend = None
    
    def createconsensus(self):
        """Create consensus of clip read at both side"""
        self.consstart = self.posstart.getconsensusstart()
        self.consend = self.posend.getconsensusend()

    def searchir(self):
        if self.consstart is None or self.consend is None:
            self.createconsensus()
        ##TODO search ir on consensus

    def __str__(self):
        return str(self.posstart) + "\t" + str(self.posend)
        
    def __len__(self):
        return self.posend.pos - self.posstart.pos

    def __repr__(self):
        return self.posstart + " - " + self.posend


class Couples():
    """Class containing all the couple"""
    def __init__(self):
        self.couples = []

    def groupeposition(self, positions, size):
        """Fonction to parse the positions and regroup start and end clip positions to a couple"""
        posstart = []
        posend = []
        chrom = None
        for pos in positions.nextposition():
            if chrom is None: ##first position
                chrom = pos.chrom
            if chrom == pos.chrom: ##same chromosome, add position to list
                if len(pos.clipstart) != 0:
                    posstart.append(pos)
                if len(pos.clipend) != 0:
                    posend.append(pos)
            else: ##change chromosome, create couple and reset
                self.__createcouples(posstart, posend, size)
                posstart = []
                posend = []
                chrom = pos.chrom
        ##add last chromomose couple
        self.__createcouples(posstart, posend, size)

    
    def __createcouples(self, posstart, posend, size):
        """Create all couple from same chromomosome"""
        for s in posstart:
            for e in posend:
                diff = (e.pos - s.pos)
                if diff >= 0 and diff <= size:
                    self.couples.append(Couple(s,e))

    def filteroverlapcouple(self):
        """Function to filter overlap couples and return best couple"""
        couplestmp = []
        for couple in self.couples:
            couplestmp.append(couple)
        self.couples = couplestmp

    def __iter__(self):
        for couple in self.couples:
            yield couple

if __name__=='__main__':
    import doctest
    doctest.testmod()
