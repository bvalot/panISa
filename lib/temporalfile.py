#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL


import tempfile
import traceback


def createfile():
    newtempfile = tempfile.NamedTemporaryFile()
    return newtempfile

        
def closefile(closefile):
    closefile.close()


def writefile(writedata):
    temp = tempfile.NamedTemporaryFile()
    try:
        temp.write(writedata)
        temp.seek(0)    
    except Exception:
        exe = traceback.format_exc()
        print exe
    return temp


def readfile(tempread,separateBy):
    readdata = []
    sep_line = ["\n","line"]
    try:
        tempread.seek(0)
        if separateBy in sep_line:
            readdata = tempread.readlines()
        else:
            readdata = tempread.read().split(separateBy)[1:]
    except Exception:
        exe = traceback.format_exc()
        print exe
    return readdata


if __name__=='__main__':
    import doctest
    doctest.testmod()