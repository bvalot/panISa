#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Search integrative element insertion on BAM alignment"""

import sys
import argparse
import os
import tempfile
from lib import bamreader
from lib.couple  import Couples
from lib.writer import writetabular

desc = "Search integrative element (IS) insertion on a genome using BAM alignment"
command = argparse.ArgumentParser(prog='panISa.py', \
    description=desc, usage='%(prog)s [options] bam')
command.add_argument('-o', '--output', nargs="?", \
    type=argparse.FileType("w"), default=sys.stdout, \
    help='Return list of IS insertion by alignment, default=stdout')
command.add_argument('-q', '--quality', nargs="?", \
    type=int, default=20, \
    help='Min alignment quality value to conserve a clipped read, default=20')
command.add_argument('-m', '--minimun', nargs="?", \
    type=int, default=10, \
    help='Min number of clipped reads to look at IS on a position, default=10')
command.add_argument('-s', '--size', nargs="?", \
    type=int, default=20, \
    help='Maximun size of direct repeat region, default=20pb')
command.add_argument('-p', '--percentage', nargs="?", \
    type=float, default=0.8, \
    help='Minimum percentage of same base to create consensus, default=0.8')
command.add_argument('bam', type=argparse.FileType("r"), \
    help='Alignment on BAM/SAM format')
command.add_argument('-v', '--version', action='version', \
    version='%(prog)s 0.1.2')


if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    output = args.output
    ##valided argument
    if args.percentage < 0 or args.percentage > 1:
        raise Exception("--percentage must be comprised between 0 to 1")
    ##Search clip read on bam by position with min quality
    positions = bamreader.parse(args.bam.name, args.quality)

    ##filter positions with not enought clipread
    positions.filterposition(args.minimun)
        
    ##find close position with both start and end clip to make a couple
    couples = Couples()
    couples.groupeposition(positions, args.size)
    couples.filteroverlapcouple()
        
    ##Create consensus and search inverted repeat sequence on it
    for couple in couples:
        couple.createconsensus(args.percentage)
        couple.searchir()
    
    ##return list of possible IS
    writetabular(output, couples)
