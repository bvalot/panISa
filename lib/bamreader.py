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
import variables

def parse(bamfile, minqual):
    bamhandle =  AlignmentFile(bamfile, 'rb')
    
    positions = Positions()
    
    for read in bamhandle:
        if isclip(read) is False or read.mapping_quality < minqual:
            continue
        clip = ClipRead(read)
        pos = positions.getposition(bamhandle.get_reference_name(read.reference_id), \
                                    clip.getclippos())
        pos.addclipread(clip)
        
    bamhandle.close()    
    return positions


def isclip(read):
    """Look at if the read is a clip read"""
    ##return False for unmapped, or full map
    if read.is_unmapped or read.cigartuples is None or len(read.cigartuples) == 1:
        return False
    ##search cigar value in start and end
    clip_start = read.cigartuples[0][0] in variables.cigarclip
    clip_end = read.cigartuples[-1][0] in variables.cigarclip
    if clip_start :
        if clip_end:
            return False ##clip read at both end not take into account
        else:
            return True
    elif clip_end:
        return True
    else:
        return False

if __name__=='__main__':
    import doctest
    doctest.testmod()
