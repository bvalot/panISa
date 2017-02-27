#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

class Couple():
    """Class containing a start and end clip position with potential IS"""
    def __init__(self, chrom, posstart, posend):
        self.chrom = chrom
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

    def rangeposition(self):
        return range(self.posstart.pos, self.posend.pos+1)

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
        """Fonction to parse the positions 
        and regroup start and end clip positions to a couple"""
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
                self.__createcouples(chrom,posstart, posend, size)
                posstart = []
                posend = []
                chrom = pos.chrom
        ##add last chromomose couple
        self.__createcouples(chrom, posstart, posend, size)

    
    def __createcouples(self, chrom, posstart, posend, size):
        """Create all couple from same chromomosome"""
        for s in posstart:
            if s.chrom != chrom:
                raise Exception("Only position with same chromosome can be coupled")
            for e in posend:
                if e.chrom != chrom:
                    raise Exception("Only position with same chromosome can be coupled")
                diff = (e.pos - s.pos)
                if diff >= 0 and diff <= size:
                    self.couples.append(Couple(chrom,s,e))

    def filteroverlapcouple(self):
        """Function to filter overlap couples and return best couple"""
        couplestmp = []
        cur_couple = None
        for couple in self.couples:
            if cur_couple is None:
                cur_couple = couple
            else:
                if self.__overlap(couple, cur_couple):
                    cur_couple = self.__best(couple, cur_couple)
                else:
                    couplestmp.append(cur_couple)
                    cur_couple = couple
        if cur_couple is not None:
            couplestmp.append(cur_couple)
        self.couples = couplestmp

    def __overlap(self, couple1, couple2):
        if couple1.chrom != couple2.chrom:
            return False
        else:
            a = set(couple1.rangeposition())
            b = set(couple2.rangeposition())
            if len(a.intersection(b)) == 0:
                return False
            else:
                return True

    def __best(self, couple1, couple2):
        diff = len(couple1.posstart.clipstart) - len(couple2.posstart.clipstart) + \
               len(couple1.posend.clipend) - len(couple2.posend.clipend)
        if diff < 0:
            return couple2
        elif diff > 0:
            return couple1
        else:
            if len(couple1) > len(couple2):
                return couple1
            else:
                return couple2
    
    def __iter__(self):
        for couple in self.couples:
            yield couple

if __name__=='__main__':
    import doctest
    doctest.testmod()
