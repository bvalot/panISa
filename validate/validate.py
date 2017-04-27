#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Validate panISa [search (IS) insertion on a genome] by simulation read"""

import argparse
import subprocess
from os import listdir
from util import siminsertseq
from util import simread
from util import managefiledir as man
from util.genreport import GenReport

desc = "Validate panISa [search (IS) insertion on a genome] by simulation read"
command = argparse.ArgumentParser(prog='validate.py', \
    description=desc, usage='%(prog)s [options] fasta_dir isinfo_dir')
command.add_argument('fasta_dir', action='store', \
    help='Input directory that store fasta file')
command.add_argument('isinfo_dir', action='store', \
    help='Input directory that store info of IS file')
command.add_argument('-o', '--output_dir', action='store', default='output/',\
    help='Place all output directories of all subprocess of validation, default=output')
command.add_argument('-s', '--seq_dir', action='store', default='output/seq/',\
    help='Return file of the simulated complete genome sequence with ISs, default=seq')
command.add_argument('-i', '--info_dir', action='store', default='output/info/',\
    help='Return report of information which insert in simulated complete genome, default=info')
command.add_argument('-r', '--report_dir', action='store', default='output/report/',\
    help='Return comparation report between panISa and real IS, default=report')
command.add_argument('-c', '--coverage', nargs="?", \
    type=int, default=60, help='Mean coverage for simulation read, default=60')
command.add_argument('-l', '--length', nargs="?", \
    type=int, default=150, help='Length of the first and second reads, default=150')
command.add_argument('-v', '--version', action='version', \
    version='%(prog)s 0.1.0')


if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    read_dir = '%sread/'%(args.output_dir)
    temp_dir = '%stemp/'%(args.output_dir) ## will be change to use temp file after be sure this scirpt is correct

    ##Prepare directory for all outputs
    man.makedir(args.output_dir)
    man.removedir(args.output_dir)
    ##Prepare directory for simulated read output
    man.makedir(read_dir)
    ##Prepare directory for report output
    man.makedir(temp_dir) ## will be change to use temp file after be sure this scirpt is correct
    man.makedir(args.report_dir)

    for sim in listdir(args.fasta_dir):
        ##Simutale the complete genome sequence with ISs
        simfasta_file = siminsertseq.getSimIS(args.isinfo_dir,args.fasta_dir+sim,args.seq_dir,args.info_dir)

        ##Move ref file and simulated complete genome as the same place 
        man.copyfile(args.fasta_dir+sim,read_dir+sim)
        ref_file = read_dir+sim
        man.copyfile(simfasta_file,read_dir+simfasta_file.split("/")[-1])
        simfasta_file = read_dir+simfasta_file.split("/")[-1]
        ##Simulate read from simulated complete genome
        simread.simread(ref_file,simfasta_file,args.length,args.coverage,read_dir)
        ##Remove intermediate files of read simulation process
        readfile = simfasta_file.split("/")[-1].split(".fasta")[0]+".bam"
        man.rmexceptfile(read_dir,readfile)

        ##Call panISa.py
        cmdIS = "../panISa.py %s -o %s%s"%(read_dir+readfile,temp_dir,sim.split(".")[0])
        subprocess.call(cmdIS, shell=True)

        ##Create comparation report
        name = sim.split(".fasta")[0]
        report = GenReport()
        report.processReport(args.info_dir+name,temp_dir+name,2,args.report_dir+name)


