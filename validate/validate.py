#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""Validate panISa [search (IS) insertion on a genome] by simulation read"""

import sys
import argparse
import tempfile
import subprocess
from os import listdir
from util import siminsertseq
from util import simread
from util import managefiledir as man
from util.genreport import GenReport

desc = "Validate panISa [search (IS) insertion on a genome] by simulation read"
command = argparse.ArgumentParser(prog='validate.py', \
    description=desc, usage='%(prog)s [options] fasta isinfo')
command.add_argument('fasta', nargs="?", type=argparse.FileType("r"),\
    help='Input fasta file')
command.add_argument('isinfo', nargs="?", type=argparse.FileType("r"),\
    help='Input tabular file that store information of ISs')
command.add_argument('-o', '--output', nargs="?", type=argparse.FileType("w"), \
    default=sys.stdout,\
    help='Return comparation report between panISa and simulation, default=stdout')
command.add_argument('-n', '--number_IS', nargs="?", type=int, default=30, \
    help='Number of ISs which were randomly added in simulation, default=30')
command.add_argument('-c', '--coverage', nargs="?", type=int, default=60, \
    help='Mean coverage for simulation read, default=60')
command.add_argument('-l', '--length', nargs="?", type=int, default=150, \
    help='Length of the first and second reads, default=150')
command.add_argument('-s', '--size', nargs="?", type=int, default=15, \
    help='Maximun size of direct repeat region, default=15')
command.add_argument('-q', '--quality', nargs="?", type=int, default=20, \
    help='Min alignment quality value to conserved a clip read, default=20')
command.add_argument('-v', '--version', action='version', \
    version='%(prog)s 0.1.0')


if __name__=='__main__':
    """Performed job on execution script""" 
    args = command.parse_args()
    
    ##[Prepare temp files and temp directory for all outputs]
    temp_dir = tempfile.mkdtemp()+"/" ##[Prepare directory for simulated read output]
    simisinfo_resultfile = tempfile.NamedTemporaryFile()
    panisa_resultfile = tempfile.NamedTemporaryFile()

    ##[Simulate the complete genome sequence with ISs]
    simfasta_file = siminsertseq.getSimIS(args.isinfo.name,args.fasta.name, \
        temp_dir,simisinfo_resultfile.name,args.number_IS)

    ##[Move original fasta file and simulated complete genome as the same place]
    ref_file = temp_dir+args.fasta.name.split("/")[-1]
    man.copyfile(args.fasta.name,ref_file)

    ##[Simulate read from simulated complete genome]
    simread.simread(ref_file,simfasta_file,args.length,args.coverage,temp_dir)

    ##[Remove intermediate files of read simulation process]
    readfile = simfasta_file.split("/")[-1].split(".fasta")[0]+".bam"
    man.rmexceptfile(temp_dir,readfile)

    ##[Call panISa.py]
    cmdIS = "../panISa.py %s -o %s -s %s -q %s"%(temp_dir+readfile,\
        panisa_resultfile.name,args.size,args.quality)
    subprocess.call(cmdIS, shell=True)

    ##[Create comparation report]
    report = GenReport()
    report.processReport(simisinfo_resultfile.name,panisa_resultfile.name,args.output)

    ##[Clear temp files]
    simisinfo_resultfile.close()
    panisa_resultfile.close()
    man.removedir(temp_dir)