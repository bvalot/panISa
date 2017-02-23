h#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL


class ClipRead():
    """Clip read object"""
    def __init__(self, alignedsegment):
        if type(alignedsegment) is not pysam.AlignedSegment:
            raise Exception("ClipRead must be initialised with a AlignedSegment object")
        self.seq = alignedsegment.seq
        ##TODO store other import variable
        
    def isstartclip(self):
        """Test if the read is start or end clip"""
        ##TODO
        return True

    def getclippos(self):
        """Return the position of the clip"""
        ##TODO
        return -1
        
    def __len__(self):
        return len(self.seq)

    def __repr__(self):
    	return self.seq
    

if __name__=='__main__':
    import doctest
    doctest.testmod()
