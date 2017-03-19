#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

def defineIR(irdata):
    if irdata is None:
        ir = "No IR"
    else:
        ir = irdata
    return str(ir)

def isfoundDR(rawdata):
    if len(rawdata.dr) > 0:##found DR
        templatedata = "\t".join([rawdata.chrom, str(rawdata.posstart.pos+1), rawdata.dr, str(rawdata.posend.pos),\
            defineIR(rawdata.ir), rawdata.cons5prime, rawdata.cons3prime])+"\n"
    else:
        templatedata = "\t".join([rawdata.chrom, str(rawdata.posend.pos), rawdata.dr, str(rawdata.posstart.pos+1),\
            defineIR(rawdata.ir), rawdata.cons5prime, rawdata.cons3prime])+"\n"
    return templatedata

def writeoutput(index,data):
    header = "\t".join(["Chromosome","Start","Direct repeats","End","Inverted repeats","Left sequences","Right sequences"])\
    +"\n"+"-"*103
    outdata = isfoundDR(data)

    #[check index for writing header]
    if index == 0:
        outdata = "\n".join([header,outdata])
    else:
        outdata = outdata
    return outdata


if __name__=='__main__':
    import doctest
    doctest.testmod()