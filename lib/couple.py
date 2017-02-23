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
        
        
    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return self.posstart + " - " + self.posend


def groupeposition(positions):
    """Fonction to parse the positions and regroup start and end clip positions to a couple"""
    couples = []
    ##TODO
    for pos in positions.nextposition():
        ##TODO
        pass
    return couples
    

if __name__=='__main__':
    import doctest
    doctest.testmod()
