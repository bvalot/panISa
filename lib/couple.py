#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import invertedrepeat as ir

class Couple():
    """Class containing a start and end clip position with potential IS"""
    def __init__(self, chrom, posstart, posend):
        self.chrom = chrom
        self.posstart = posstart
        self.posend = posend
        self.cons5prime = None
        self.cons3prime = None
        self.dr = None
        self.ir = None
        
    def createconsensus(self, percent):
        """Create consensus of clip read at both side"""
        self.cons3prime = self.posstart.getconsensus(percent, True)
        self.cons5prime = self.posend.getconsensus(percent, False)
        self.dr = self.__createdr(percent)
        
    def searchir(self):
        if self.cons5prime is None or self.cons3prime is None:
            raise Exception("You must create consensus before search ir")
        self.ir = ir.searchir(self.cons5prime, self.cons3prime)

    def rangeposition(self):
        return range(self.posstart.pos, self.posend.pos+1)

    def __createdr(self, percent):
        size = self.posend.pos - self.posstart.pos
        if size == 0:
            return ""
        seqs = []
        for clip in self.posstart.clipstart:
            s = clip.getdr(self.posstart.pos, self.posend.pos)
            if len(s) == self.posend.pos - self.posstart.pos:
                seqs.append(s)
        for clip in self.posend.clipend:
            s = clip.getdr(self.posstart.pos, self.posend.pos)
            if len(s) == self.posend.pos - self.posstart.pos:
                seqs.append(s)
        if len(seqs) == 0 :
            return "No sequence to build direct repeat"
        if len(set(map(len,seqs))) != 1:
            raise Exception("Direct repeat region can only compute on same length")
        dr_tmp = []
        for i in range(0,len(seqs[0])):
            allele = set()
            values = []
            for seq in seqs:
                base = seq[i]
                allele.add(base)
                values.append(base)
            toadd = "N"
            for base in allele:
                if len(values) >= 2 and \
                   float(values.count(base))/len(values) >= percent:
                    toadd = base
            dr_tmp.append(toadd)
        return "".join(dr_tmp)
    
    def __str__(self):
        toreturn = str(self.posstart) + "\t" + str(self.posend) + "\t" + self.dr + "\t"
        if self.ir is None:
            toreturn += "No IR"
        else:
            toreturn += str(self.ir)
        toreturn += "\t" + self.cons5prime + "\t" + self.cons3prime
        return  toreturn
        
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
            elif chrom != pos.chrom: ##change chromosome, create couple and reset
                self.__createcouples(chrom, posstart, posend, size)
                posstart = []
                posend = []
                chrom = pos.chrom
            ##add position to list
            if len(pos.clipstart) != 0:
                posstart.append(pos)
            if len(pos.clipend) != 0:
                posend.append(pos)
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
