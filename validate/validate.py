#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Validate panISa [search (IS) insertion on a genome] by simulation read"""

import sys
import argparse
from shutil import copyfile
from lib import siminsertseq
from lib import simread

desc = "Validate panISa [search (IS) insertion on a genome] by simulation read"
command = argparse.ArgumentParser(prog='validate.py', \
    description=desc, usage='%(prog)s [options] fasta_dir isinfo_dir')
command.add_argument('fasta_dir', action='store', \
    help='Input directory that store fasta file')
command.add_argument('isinfo_dir', action='store', \
    help='Input directory that store info of IS file')
command.add_argument('-s', '--seq_dir', action='store', default='output/seq/',\
    help='Return file of the simulated complete genome sequence with ISs, default=seq')
command.add_argument('-i', '--info_dir', action='store', default='output/info/',\
    help='Return report of information which insert in simulated complete genome, default=info')
command.add_argument('-r', '--ref_file', default='NC_002516.fasta', \
    help='Reference file in fasta format for BAM read simulation, default=NC_002516.fasta')
command.add_argument('-c', '--coverage', nargs="?", \
    type=int, default=50, help='Mean coverage for simulation read, default=50')
command.add_argument('-l', '--length', nargs="?", \
    type=int, default=150, help='Length of the first and second reads, default=150')
command.add_argument('-v', '--version', action='version', \
    version='%(prog)s 0.1.0')


if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    fasta_dir = args.fasta_dir
    isinfo_dir = args.isinfo_dir
    outseq_dir = args.seq_dir
    info_dir = args.info_dir

    ##Simutale the complete genome sequence with ISs
    simfasta_file = siminsertseq.simSequences(isinfo_dir,fasta_dir,outseq_dir,info_dir)


    ##Simulate read from simulated complete genome
    read_dir = 'output/seq/'
    ##Move ref file as same place with simulated complete genome
    copyfile(fasta_dir+args.ref_file,read_dir+args.ref_file)
    args.ref_file = read_dir+args.ref_file
    simread.simread(args.ref_file,simfasta_file,args.length,args.coverage)

