#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Read Bam File and store clip read by position"""

from pysam import AlignmentFile
from clipread import ClipRead
from position import Position, Positions

def parse(bamfile, minqual):
    bamhandle =  AlignmentFile(bamfile, 'rb')
    
    positions = Positions()
    
    for read in bamhandle:
        if isclip(read) is False:
            continue
        clip = ClipRead(read)
        pos = positions.getposition(bamhandle.get_reference_name(read.reference_id), \
                                    clip.getclippos())
        pos.addclipread(clip)
        
    bamhandle.close()
    
    return positions


def isclip(read):
    """Look at if the read is a clip read"""
    ##TODO
    return False


if __name__=='__main__':
    import doctest
    doctest.testmod()
