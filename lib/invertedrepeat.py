#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

class InvertRepeat():
    """Class containing inverted repeat at the 5' and the 3'"""
    def __init__(self):
        self.prime5_pos= None ##position beginning ir from clip end consensus (1)
        self.prime5_seq = None ##Sequence ir from clip end consensus
        self.prime3_pos = None ##position - len(seq) -1 beginning ir from clip start consensus (-1)
        self.prime3_seq = None ##sequence ir from clip start consensus (reverse)

    def __str__(self):
        if self.prime5_pos is None:
            return "No IR"
        else :
            return " ".join(self.prime5_pos, self.prime5_seq, self.prime3_pos, self.prime3_seq)


def searchir(prime5seq, prime3seq):
    """Searched ir from 5 to 3 prie consensus and return InvertRepeat class"""
    ##TODO by super Panisa
    ##Only store IR that is close to the 5 and 3 prime (<10bp)
    ##test
    return InvertRepeat()
