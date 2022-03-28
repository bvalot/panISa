#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##Copyright (c) 2019 Benoit Valot and Charlotte Couchoud
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""automate annotated an ISFinder output with gene informations"""

import argparse
import sys
from lib import gff
from lib import variables

desc = "annotated an ISFinder result with gene informations"
command = argparse.ArgumentParser(prog='annotateISresult.py', \
    description=desc, usage='%(prog)s [options] ISfinder gff')
command.add_argument('-o', '--output', nargs="?", \
    type=argparse.FileType("wt"), default=sys.stdout, \
    help='Return annotated ISFinder result (default:stdout)')
command.add_argument('-f', '--feature', nargs="?", \
    type=str, default="CDS", \
    help='Defined the feature to screen in gff, default:CDS')
command.add_argument('-i', '--id', nargs="?", \
    type=str, default="ID", \
    help='Defined the key in gff containing ID, default:ID')
command.add_argument('-a', '--annote', nargs="?", \
    type=str, default="product", \
    help='Defined the key in gff containing annotation, default:product')
command.add_argument('isfinder', type=argparse.FileType("r"), \
    help='ISFinder search result to annotate')
command.add_argument('gff', type=argparse.FileType("r"), \
    help='Gff file containing annotation of genome')
command.add_argument('-v', '--version', action='version', \
    version='%(prog)s ' + variables.version)

def splitline(line):
    return line.rstrip("\n").split("\t")

def searchInside(chrom, start, end, gffs):
    for g in gffs:
        if g.is_in(chrom, start) or g.is_in(chrom, end):
            return g
    return None

def searchSurround(chrom, start, end, gffs):
    prev_, next_ = None, None
    for g in gffs:
        if g.is_closer_previous(chrom, start, prev_):
            prev_ = g
        if g.is_closer_next(chrom, end, next_):
            next_ = g
    return prev_, next_

def printResult(res, inside, prev_, next_, ident, annot):
    """Add information of ID, start, stop, strand, annotation"""
    res2 = list(res)
    addAnnot(res2, prev_, ident, annot)
    addAnnot(res2, inside, ident, annot)
    addAnnot(res2, next_, ident, annot)
    return res2
        
def addAnnot(res2, gff, ident, annot):
    if gff:
        res2.append(gff.get_annotation(ident))
        res2.append(gff.get_annotation("Start"))
        res2.append(gff.get_annotation("End"))
        res2.append(gff.get_annotation("Strand"))
        res2.append(gff.get_annotation(annot))
    else:
        res2.extend(["", "", "", "", ""])

def headerOut(header):
    headerOut = list(header)
    tmp = ["ID", "Start", "End", "Strand", "Fonction"]
    headerOut.extend(["previous_" + t for t in tmp])
    headerOut.extend(["inside_" + t for t in tmp])
    headerOut.extend(["next_" + t for t in tmp])
    return headerOut

if __name__=='__main__':
    """Performed job on execution script"""
    args = command.parse_args()

    ## load gffs informations
    gffs = gff.read_gff(args.gff, args.feature)
    # print(str(len(gffs)))

    ## Verify  ISFinder input
    ## Sample Chromosome Start_Position Stop_Position Potential_sequence Potential_IS Alignment
    header = splitline(args.isfinder.readline())
    if len(header) != 7 and header[0] != "Sample":
        raise Exception("The file seems to not be an ISFinder result")

    ## start output
    output = args.output
    output.write("\t".join(headerOut(header)) + "\n")
    
    ## read ISFinder result and return annotated one
    for line in args.isfinder.readlines():
        res = splitline(line)
        ## IS are inside gene
        inside = searchInside(res[1], int(res[2]), int(res[3]), gffs)
        if inside:
            # print(res[1] + "|" + res[2]  + " " + str(inside))
            output.write("\t".join(printResult(res, inside, None, None, args.id, args.annote)) + "\n")
        else:
            prev_, next_ = searchSurround(res[1], int(res[2]), int(res[3]), gffs)
            # print(res[1] + "|" + res[2]  + " " + str(prev_) + " " + str(next_))
            output.write("\t".join(printResult(res, None, prev_, next_, args.id, args.annote)) + "\n")
