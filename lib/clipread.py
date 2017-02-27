#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

import variables

class ClipRead():
    """Clip read object"""
    def __init__(self, alignedsegment):
        self.read_seq = alignedsegment.query_sequence
        self.read_name = alignedsegment.query_name
        self.read_start = alignedsegment.query_alignment_start
        self.read_end = alignedsegment.query_alignment_end
        self.read_len = alignedsegment.query_alignment_length
        self.ref_start = alignedsegment.reference_start
        self.ref_end = alignedsegment.reference_end
        self.ref_len = alignedsegment.reference_length
        self.cigar = alignedsegment.cigarstring
        self.cigartuples = alignedsegment.cigartuples
        self.isreverse = alignedsegment.is_reverse
        
    def isstartclip(self):
        """Test if the read is start or end clip, look at """
        if self.cigartuples is None:
            raise Exception("ClipRead must be aligned")
        if self.cigartuples[0][0] in variables.cigarclip:
            return True
        elif self.cigartuples[-1][0] in variables.cigarclip:
            return False
        else:
            raise Exception("ClipRead must contain clip part at start or end")

    def getclippos(self):
        """Return the position of the clip"""
        if self.isstartclip():
            return self.ref_start
        else:
            return self.ref_end

    def getclipseq(self):
        """return clip part of the read, except for hard clip return None"""
        if len(self.read_seq) == self.read_len:
            return None
        if self.isstartclip():
            return self.read_seq[:self.read_start]
        else:
            return self.read_seq[self.read_end:]
        
    def __len__(self):
        return len(self.read_seq)

    def __repr__(self):
    	return self.read_seq
    

if __name__=='__main__':
    import doctest
    doctest.testmod()
